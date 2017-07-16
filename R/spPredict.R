spPredict <- function(sp.obj, X.0, coords.0, start=1, end, thin=1, n.omp.threads = 1, verbose=TRUE, n.report=100, ...){
  
  ####################################################
  ##Check for unused args
  ####################################################
    formal.args <- names(formals(sys.function(sys.parent())))
    elip.args <- names(list(...))
    for(i in elip.args){
        if(! i %in% formal.args)
            warning("'",i, "' is not an argument")
    }
    
    if(missing(sp.obj)){stop("error: spPredict expects sp.obj\n")}
    if(!class(sp.obj) %in% c("rNNGP","sNNGP")){
        stop("error: requires an output object of class rNNGP or sNNGP\n")
    }
    
    obj.class <- class(sp.obj)
    
    X <- sp.obj$X
    y <- sp.obj$y
    n <- nrow(X)
    p <- ncol(X)
    p.theta.samples <- sp.obj$p.theta.samples
    p.beta.samples <- sp.obj$p.beta.samples
    n.samples <- nrow(p.beta.samples)
    if(obj.class == "sNNGP"){
        p.w.samples <- sp.obj$p.w.samples
    }    
    n.neighbors <- sp.obj$n.neighbors
    coords <- sp.obj$coords
    ord <- sp.obj$ord
    cov.model.indx <- sp.obj$cov.model.indx

    ##subsamples
    if(missing(end)){end <- n.samples}
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
       
    s.indx <- seq(start, end, by=as.integer(thin))

    p.theta.samples <- t(p.theta.samples[s.indx,,drop=FALSE])
    p.beta.samples <- t(p.beta.samples[s.indx,,drop=FALSE])
    n.samples <- ncol(p.beta.samples)
    if(obj.class == "sNNGP"){
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
    if(obj.class == "sNNGP"){
        storage.mode(p.w.samples) <- "double"
    }
    storage.mode(n.samples) <- "integer"
    storage.mode(cov.model.indx) <- "integer"
    storage.mode(nn.indx.0) <- "integer"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    
    ptm <- proc.time()
    
    if(obj.class == "sNNGP"){
        out <- .Call("sNNGPPredict", X, y, coords, n, p, n.neighbors, X.0, coords.0, q, nn.indx.0, 
                     p.beta.samples, p.theta.samples, p.w.samples, n.samples, cov.model.indx, n.omp.threads, verbose, n.report)
    }else{
        out <- .Call("rNNGPPredict", X, y, coords, n, p, n.neighbors, X.0, coords.0, q, nn.indx.0, 
                     p.beta.samples, p.theta.samples, n.samples, cov.model.indx, n.omp.threads, verbose, n.report)
    }
    
    out$run.time <- proc.time() - ptm
    
    out
}
