spPredict <- function(sp.obj, X.0, coords.0, sub.sample, n.omp.threads = 1, verbose=TRUE, n.report=100, ...){
  
    ####################################################
    ##Check for unused args
    ####################################################
    formal.args <- names(formals(sys.function(sys.parent())))

    elip.args <- list(...)
    for(i in names(elip.args)){
        if(! i %in% formal.args)
            warning("'",i, "' is not an argument")
    }
    
    if(missing(sp.obj)){stop("error: spPredict expects sp.obj\n")}
    if(!class(sp.obj)[1] == "NNGP"){
        stop("error: requires an output object of class NNGP\n")
    }


    ##call
    out <- list()
    out$call <- match.call()
    out$sp.obj.class <- class(sp.obj)
    
    ##conjugate
    if(class(sp.obj)[2] == "conjugate"){

        if(!missing(sub.sample)){
            warning("'sub.sample' is not an argument for prediction using a conjugate model.")
        }
        
        theta.alpha <- as.vector(sp.obj$theta.alpha)
        names(theta.alpha) <- colnames(sp.obj$theta.alpha)
        
        ptm <- proc.time()
        if(length(class(sp.obj)) == 4 & class(sp.obj)[4] == "SLGP"){
            out <- c(out, spConjNNGP(sp.obj$y ~ sp.obj$X-1, coords=sp.obj$coords, knots=sp.obj$knots, sigma.sq.IG=sp.obj$sigma.sq.IG, n.neighbors=sp.obj$n.neighbors,
                                     X.0 = X.0, coords.0=coords.0,
                                     theta.alpha=theta.alpha, cov.model=sp.obj$cov.model, n.omp.threads=n.omp.threads, search.type=sp.obj$search.type))
        }else{
            out <- c(out, spConjNNGP(sp.obj$y ~ sp.obj$X-1, coords=sp.obj$coords, sigma.sq.IG=sp.obj$sigma.sq.IG, n.neighbors=sp.obj$n.neighbors,
                                     X.0 = X.0, coords.0=coords.0,
                                     theta.alpha=theta.alpha, cov.model=sp.obj$cov.model, n.omp.threads=n.omp.threads, search.type=sp.obj$search.type))
        }
        
        out$run.time <- proc.time() - ptm
    }else{ ##sequential and response models
        
        X <- sp.obj$X
        y <- sp.obj$y
        coords <- sp.obj$coords
        family <- sp.obj$family ##family code gaussian=1, binomial=2, ...
        
        family.indx <- 1
        if(class(sp.obj)[3] == "binomial"){
            family.indx <- 2
        }
        
        n <- nrow(X)
        p <- ncol(X)
        
        p.theta.samples <- sp.obj$p.theta.samples
        p.beta.samples <- sp.obj$p.beta.samples
        n.samples <- nrow(p.beta.samples)
        if(class(sp.obj)[2] == "sequential"){
            p.w.samples <- sp.obj$p.w.samples
        }    
        n.neighbors <- sp.obj$n.neighbors
        cov.model.indx <- sp.obj$cov.model.indx
        
        ##subsamples
        if(missing(sub.sample)){
            sub.sample <- list()
        }
     
        start <- ifelse(!"start" %in% names(sub.sample), floor(0.5*n.samples), sub.sample$start)
        end <- ifelse(!"end" %in% names(sub.sample), n.samples, sub.sample$end)
        thin <- ifelse(!"thin" %in% names(sub.sample), 1, sub.sample$thin)   
        if(!is.numeric(start) || start >= n.samples){stop("invalid start")}
        if(!is.numeric(end) || end > n.samples){stop("invalid end")}
        if(!is.numeric(thin) || thin >= n.samples){stop("invalid thin")}
        sub.sample <- list(start=start, end=end, thin=thin)
        s.indx <- seq(as.integer(start), as.integer(end), by=as.integer(thin))
        n.samples <- length(s.indx)
        
        p.theta.samples <- t(p.theta.samples[s.indx,,drop=FALSE])
        p.beta.samples <- t(p.beta.samples[s.indx,,drop=FALSE])
        
        if(class(sp.obj)[2] == "sequential"){
            p.w.samples <- p.w.samples[,s.indx,drop=FALSE]
        }    
        
        ##check X.0 and coords.0
        if(missing(X.0)){stop("error: X.0 must be specified\n")}
        if(!any(is.data.frame(X.0), is.matrix(X.0))){stop("error: X.0 must be a data.frame or matrix\n")}
        if(ncol(X.0) != ncol(X)){ stop(paste("error: X.0 must have ",p," columns\n"))}
        
        if(missing(coords.0)){stop("error: coords.0 must be specified\n")}
        if(!any(is.data.frame(coords.0), is.matrix(coords.0))){stop("error: coords.0 must be a data.frame or matrix\n")}
        if(!ncol(coords.0) == 2){stop("error: coords.0 must have two columns\n")}
        
        q <- nrow(X.0)
        
        ##get nn indx
        nn.indx.0 <- nn2(coords, coords.0, k=n.neighbors)$nn.idx-1 ##obo for cNNGP.cpp indexing
        
        storage.mode(X) <- "double"
        storage.mode(y) <- "double"
        storage.mode(coords) <- "double"
        storage.mode(n) <- "integer"
        storage.mode(p) <- "integer"
        storage.mode(n.neighbors) <- "integer"
        storage.mode(X.0) <- "double"
        storage.mode(coords.0) <- "double"
        storage.mode(q) <- "integer"
        storage.mode(p.beta.samples) <- "double"
        storage.mode(p.theta.samples) <- "double"
        if(class(sp.obj)[2] == "sequential"){
            storage.mode(p.w.samples) <- "double"
        }
        storage.mode(n.samples) <- "integer"
        storage.mode(cov.model.indx) <- "integer"
        storage.mode(nn.indx.0) <- "integer"
        storage.mode(n.omp.threads) <- "integer"
        storage.mode(verbose) <- "integer"
        storage.mode(n.report) <- "integer"
        storage.mode(family.indx) <- "integer"
        
        ptm <- proc.time()
        
        if(class(sp.obj)[2] == "sequential"){
            out <- c(out, .Call("sNNGPPredict", X, y, coords, n, p, n.neighbors, X.0, coords.0, q, nn.indx.0, 
                                p.beta.samples, p.theta.samples, p.w.samples, n.samples, family.indx, cov.model.indx, n.omp.threads, verbose, n.report))
        }else{
            out <- c(out, .Call("rNNGPPredict", X, y, coords, n, p, n.neighbors, X.0, coords.0, q, nn.indx.0, 
                                p.beta.samples, p.theta.samples, n.samples, cov.model.indx, n.omp.threads, verbose, n.report))
        }
        
        out$run.time <- proc.time() - ptm
    }

    class(out) <- "spPredict"
    out
}
