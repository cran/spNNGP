spConjNNGP <- function(formula, data = parent.frame(), coords, knots, n.neighbors = 15,
                       theta.alpha, sigma.sq.IG, cov.model = "exponential",
                       k.fold = 5, score.rule = "crps",
                       X.0, coords.0, 
                       n.omp.threads = 1, search.type = "cb", ord, return.neighbor.info = TRUE,  
                       neighbor.info, fit.rep = FALSE, n.samples, verbose = TRUE, ...){

    
    ####################################################
    ##Check for unused args
    ####################################################
    formal.args <- c(names(formals(sys.function(sys.parent()))), "n.report")
    
    elip.args <- list(...)
    for(i in names(elip.args)){
        if(! i %in% formal.args)
            warning("'",i, "' is not an argument")
    }

    ##call
    cl <- match.call()
    
    ####################################################
    ##Formula
    ####################################################
    if(missing(formula)){stop("error: formula must be specified")}
    
    if(is(formula, "formula")){
        
        holder <- parseFormula(formula, data)
        y <- as.vector(holder[[1]])
        X <- as.matrix(holder[[2]])
        x.names <- holder[[3]]
        
    }else{
        stop("error: formula is misspecified")
    }
    
    p <- ncol(X)
    n <- nrow(X)
    
    ##Coords
    if(missing(coords)){stop("error: coords must be specified")}
    
    if(is.vector(coords)){
        if(is.character(coords)){
            if(all(coords %in% colnames(data))){
                coords <- as.matrix(data[,coords])
            }else{
                stop(paste0("error: coords name ", paste(coords[!(coords %in% colnames(data))], collapse=" and "), " not in data"))
            }
        }else if(all(coords %in% (1:ncol(data)))){
            coords <- as.matrix(data[,coords])
        }else{
            stop(paste0("error: coords column index ", paste(coords[!(coords %in% (1:ncol(data)))], collapse=" and "), " not in data"))
        }
    }else{
        if(!any(is.matrix(coords), is.data.frame(coords))){
            stop("error: coords must n-by-m matrix or dataframe of coordinates or vector indicating the column names or integer indexes in data")
        }
        coords <- as.matrix(coords)
    }
    
    if(nrow(coords) != n || ncol(coords) != 2){
        stop("error: either coords has more than two columns or the number of rows is different than data used in the model formul")
    }

    ## This can be slow if n is large, the onus is on the user to check
    ## if(any(duplicated(coords))){
    ##     stop("error: duplicated coordinates found. Remove duplicates.")
    ## }
    
    ####################################################
    ##Fitted, replicated, and exact samples
    ####################################################
    if(fit.rep){
        if(missing(n.samples)){
            stop("error: if fit.rep=TRUE, specify the number of fitted and replicated samples to generate using the n.samples argument.")
        }
    }

    if(missing(n.samples)){
        n.samples <- 0
    }
    
    ####################################################
    ##Ordering
    ####################################################
    neighbor.info.provided <- FALSE
        
    if(!missing(neighbor.info)){

        if(!all(c("n.neighbors","nn.indx","nn.indx.lu","ord") %in% names(neighbor.info))){stop("The supplied neighbor.info is malformed.")}
        
        nn.indx <- neighbor.info$nn.indx
        nn.indx.lu <- neighbor.info$nn.indx.lu
        ord <- neighbor.info$ord
        n.neighbors <- neighbor.info$n.neighbors
        neighbor.info.provided <- TRUE
        
        warning("Using user defined neighbor.info, no checks are done on the supplied neighbor information.")
        
    }else{
        
        if(missing(ord)){   
            ord <- order(coords[,1])##default order by x column
        }else{
            if(length(ord) != n){stop("error: supplied order vector ord must be of length n")}
            ## if(search.type == "cb"){
            ##     warning("switching to search.type='brute' given user defined ordering ord, this could be slow if n is large")
            ##     search.type <- "brute"
            ## }
        }
    }
    
    coords <- coords[ord,]
    X <- X[ord,,drop=FALSE]
    y <- y[ord]

    ##SLGP knots
    slgp <- FALSE
    
    if(!missing(knots)){
        slgp <- TRUE
        
        if(!is.matrix(knots)){stop("error: knots must r-by-2 matrix of xy-coordinate locations")}
        if(ncol(knots) != 2){
            stop("error: either coords has more than two columns or the number of rows is different than
          data used in the model formula")
        }
        
        r <- nrow(knots)
        storage.mode(r) <- "integer"
        storage.mode(knots) <- "double"
    }
    
    storage.mode(y) <- "double"
    storage.mode(X) <- "double"
    storage.mode(p) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(coords) <- "double"

    ####################################################
    ##Prediction data
    ####################################################
    n.0 <- 0 ##if n.0 is zero in cNNGP.cpp then no prediction occurs
    nn.indx.0 <- 0
    
    if(!missing(X.0)){
        if(!is.matrix(X.0) || ncol(X.0) != p){stop(paste("error: X.0 must n.0-by-",p," matrix", sep=""))}
        n.0 <- nrow(X.0)
        if(missing(coords.0)){stop("error: coords.0 must specified")}
        if(ncol(coords.0) != 2 || nrow(coords.0) != n.0){
            stop("error: either coords.0 has more than two columns or the number of rows is different than
          X.0")
        }

        nn.indx.0 <- nn2(coords, coords.0, k=n.neighbors)$nn.idx-1 ##obo for cNNGP.cpp indexing
        
    }else{
        coords.0 <- 0
        X.0 <- 0
    }

    storage.mode(n.0) <- "integer"
    storage.mode(nn.indx.0) <- "integer"
    storage.mode(coords.0) <- "double"
    storage.mode(X.0) <- "double"
    
    ####################################################
    ##Covariance model
    ####################################################
    if(missing(cov.model)){stop("error: cov.model must be specified")}

    cov.model.names <- c("exponential","spherical","matern","gaussian")##order much match util.cpp spCor
    
    if(!cov.model%in%cov.model.names)
    {stop("error: specified cov.model '",cov.model,"' is not a valid option; choose from ", paste(cov.model.names, collapse=", ", sep="") ,".")}
    
    cov.model.indx <- which(cov.model == cov.model.names)-1##obo for cov model lookup on c side
    storage.mode(cov.model.indx) <- "integer"
    
    ####################################################
    ##Parameters
    ####################################################
    if(any(is.matrix(theta.alpha), is.vector(theta.alpha))){
        if(is.vector(theta.alpha)){
            theta.alpha <- t(as.matrix(theta.alpha))
            k.fold <- NA
            ##message("Only one set of parameters given in theta.alpha, ignoring k.fold and score.rule")
        }else{##is matrix
            if(any(is.na(k.fold), k.fold < 2, k.fold >= n)){
                stop("error: k.fold must be greater than 1 and less than n")
            }
        }   
    }else{
        stop("error: theta.alpha must be a vector or matrix")
    }
    
    if(!"phi"%in%colnames(theta.alpha)){stop("error: column name phi must be specified in theta.alpha")}
    if(!"alpha"%in%colnames(theta.alpha)){stop("error: column name alpha must be specified in theta.alpha")}

    if(cov.model == "matern"){
        if(!"nu"%in%colnames(theta.alpha)){stop("error: column name nu must be specified in theta.alpha")}
        
        ##reorder for cNNGP.cpp
        theta.alpha <- t(theta.alpha[,c("phi","alpha","nu"), drop=FALSE])
    }else{
        ##reorder for cNNGP.cpp
        theta.alpha <- t(theta.alpha[,c("phi","alpha"), drop=FALSE])
    }

    g <- ncol(theta.alpha)

    storage.mode(theta.alpha) <- "double"
    storage.mode(g) <- "integer"
    
    ####################################################
    ##Priors
    ####################################################
    if(missing(sigma.sq.IG)) {stop("error: sigma.sq.IG must be specified")}
    
    if(!is.vector(sigma.sq.IG) || length(sigma.sq.IG) != 2){stop("error: sigma.sq.IG must be a vector of length 2")}
    if(any(sigma.sq.IG <= 0)){stop("error: sigma.sq.IG must be a positive vector of length 2")}

    storage.mode(sigma.sq.IG) <- "double"
    
    ####################################################
    ##nn search 
    ####################################################
    if(!neighbor.info.provided){
        
        search.type.names <- c("brute", "cb")
        
        if(!search.type %in% search.type.names){
            stop("error: specified search.type '",search.type,"' is not a valid option; choose from ", paste(search.type.names, collapse=", ", sep="") ,".")
        }
        
        ##indexes
        if(search.type == "brute"){
            indx <- mkNNIndx(coords, n.neighbors, n.omp.threads)
        }else{
            indx <- mkNNIndxCB(coords, n.neighbors, n.omp.threads)
        }
        
        nn.indx <- indx$nnIndx
        nn.indx.lu <- indx$nnIndxLU
        nn.indx.run.time <- indx$run.time
    }
    
    storage.mode(nn.indx) <- "integer"
    storage.mode(nn.indx.lu) <- "integer"
                      
    ####################################################
    ##Other stuff
    ####################################################
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(n.neighbors) <- "integer"
    storage.mode(verbose) <- "integer"

    ####################################################
    ##Pack it up and off it goes
    ####################################################
    ptm <- proc.time()
    
    if(!is.na(k.fold)){
                
        score.rule.names <- c("rmspe","crps")
        
        if(!score.rule%in%score.rule.names){
            stop("error: specified score.rule '",score.rule,"' is not a valid option; choose from ",
                 paste(score.rule.names, collapse=", ", sep="") ,".")
        }
        
        k.fold <- as.integer(k.fold)
        k.fold.scores <- matrix(0, g, 2)
        colnames(k.fold.scores) <- score.rule.names

        if(verbose){
            cat("----------------------------------------\n")
            cat("\tStarting k-fold\n")
            cat("----------------------------------------\n")
            
            bar <- txtProgressBar(min=0, max=k.fold, initial=0, style=3, char="*", width=30)
        }
        
        for(i in 1:k.fold){

            ho.id <- sort(sample(1:n, floor(length(y)/k.fold)))
            
            X.mod <- X[-ho.id,,drop=FALSE]
            y.mod <- y[-ho.id]
            coords.mod <- coords[-ho.id,,drop=FALSE]
            n.mod <- nrow(X.mod)

            X.ho <- X[ho.id,,drop=FALSE]
            y.ho <- y[ho.id]
            coords.ho <- coords[ho.id,,drop=FALSE]
            n.ho <- nrow(X.ho)
            
            ##indexes
            if(search.type == "brute"){
                indx <- mkNNIndx(coords.mod, n.neighbors, n.omp.threads)
            }else{
                indx <- mkNNIndxCB(coords.mod, n.neighbors, n.omp.threads)
            }

            nn.indx.mod <- indx$nnIndx
            nn.indx.lu.mod <- indx$nnIndxLU

            nn.indx.mod.ho <- nn2(coords.mod, coords.ho, k=n.neighbors)$nn.idx-1 ##obo for cNNGP.cpp indexing

            storage.mode(X.mod) <- "double"
            storage.mode(y.mod) <- "double"
            storage.mode(coords.mod) <- "double"
            storage.mode(n.mod) <- "integer"
            storage.mode(X.ho) <- "double"
            storage.mode(y.ho) <- "double"
            storage.mode(coords.ho) <- "double"
            storage.mode(n.ho) <- "integer"
            storage.mode(nn.indx.mod.ho) <- "integer"
            storage.mode(nn.indx.mod) <- "integer"
            storage.mode(nn.indx.lu.mod) <- "integer"
            
            fold.verbose <- 0
            storage.mode(fold.verbose) <- "integer"

            get.X.str <- 0##only need this if we want to compute fitted and replicated data for SLGP, we don't want fit and rep for k-fold only final run
            storage.mode(get.X.str) <- "integer"

            if(slgp){
                out <- .Call("cSLGP", y.mod, X.mod, p, n.mod, r, coords.mod, knots, theta.alpha, nn.indx.mod, nn.indx.lu.mod, 
                             X.ho, coords.ho, n.ho, nn.indx.mod.ho, g,
                             n.neighbors, sigma.sq.IG, cov.model.indx, n.omp.threads, fold.verbose, get.X.str)
            }else{
                out <- .Call("cNNGP", y.mod, X.mod, p, n.mod, coords.mod, theta.alpha, nn.indx.mod, nn.indx.lu.mod, 
                             X.ho, coords.ho, n.ho, nn.indx.mod.ho, g,
                             n.neighbors, sigma.sq.IG, cov.model.indx, n.omp.threads, fold.verbose)
            }
            
            for(j in 1:g){
                k.fold.scores[j,"crps"] <- k.fold.scores[j,"crps"]+crps(y.ho, out$y.0.hat[,j], out$y.0.hat.var[,j])
                k.fold.scores[j,"rmspe"] <- k.fold.scores[j,"rmspe"]+rmspe(y.ho, out$y.0.hat[,j])
            }

            if(verbose){
                setTxtProgressBar(bar,i)
            }
        }
        if(verbose){
            cat("\n")
        }            
        
        ##calculate lowest mean score and set theta.alpha 
        k.fold.scores <- cbind(t(theta.alpha),k.fold.scores/k.fold)

        theta.alpha <- theta.alpha[,which.min(k.fold.scores[,score.rule])]
        storage.mode(theta.alpha) <- "double"
    }

    ##final run
    g <- 1
    storage.mode(g) <- "integer"
    
    ##get.X.str <- 0##only need this if we want to compute fitted and replicated data for SLGP
    ##if(fit.rep){
    get.X.str <- 1 ##just always return it for subsequent calls to fitted (just note it could be large if n and/or r are large
    ##}
    storage.mode(get.X.str) <- "integer"
    
    if(slgp){
        out <- .Call("cSLGP", y, X, p, n, r, coords, knots, theta.alpha, nn.indx, nn.indx.lu, 
                     X.0, coords.0, n.0, nn.indx.0, g,
                     n.neighbors, sigma.sq.IG, cov.model.indx, n.omp.threads, verbose, get.X.str)
        
        out$X.str <- matrix(out$X.str, nrow=n)
    }else{
        out <- .Call("cNNGP", y, X, p, n, coords, theta.alpha, nn.indx, nn.indx.lu, 
                     X.0, coords.0, n.0, nn.indx.0, g,
                     n.neighbors, sigma.sq.IG, cov.model.indx, n.omp.threads, verbose)
        
    }
    out$run.time <- proc.time() - ptm
 
    out$theta.alpha <- t(theta.alpha)
    
    out$sigma.sq.hat <- out$ab[2,]/(out$ab[1,]-1)
    out$sigma.sq.var <- out$ab[2,]^2/((out$ab[1,]-1)^2*(out$ab[1,]-2))
    out$beta.hat <- t(out$beta.hat)
    if(slgp){
        beta.var <- matrix(out$bB.inv, p+r, p+r)
    }else{
        beta.var <- matrix(out$bB.inv, p, p)
    }
    beta.var[upper.tri(beta.var, FALSE)] <- t(beta.var)[upper.tri(beta.var, FALSE)]
    out$beta.var <- out$ab[2,]*beta.var/(out$ab[1,]-1)
    if(slgp){
        colnames(out$beta.hat) <- c(x.names, paste0("w.str.",1:r))
        rownames(out$beta.var) <- colnames(out$beta.var) <- c(x.names, paste0("w.str.",1:r))
    }else{
        colnames(out$beta.hat) <- x.names
        rownames(out$beta.var) <- colnames(out$beta.var) <- x.names
    }
    
    out$sigma.sq.IG <- sigma.sq.IG
    out$n.neighbors <- n.neighbors
    out$cov.model <- cov.model
    out$cov.model.indx <- cov.model.indx
    out$search.type <- search.type
    out$call <- cl
    
    if(!is.na(k.fold)){
        out$k.fold.scores <- k.fold.scores
    }

    if(return.neighbor.info){
        out$neighbor.info <- list(n.neighbors = n.neighbors, n.indx=mk.n.indx.list(nn.indx, n, n.neighbors),
                                  nn.indx=nn.indx, nn.indx.lu=nn.indx.lu, ord=ord,
                                  nn.indx.run.time=nn.indx.run.time)
    }

    ##do exact sampling and fit and replicated data if requested
    if(n.samples > 0){
        
        b <- out$bb
        B.inv <- matrix(out$bB.inv, length(b), length(b))
        B.inv[upper.tri(B.inv, FALSE)] <- t(B.inv)[upper.tri(B.inv, FALSE)]
        alpha <- out$theta.alpha[,"alpha"]
        out$p.beta.theta.samples <- matrix(0, nrow=n.samples, ncol=length(b)+2)
        
        for(i in 1:n.samples){
            sigma.sq <- rigamma(1, out$ab[1,], out$ab[2,])
            beta <- rmvn(1, B.inv%*%b, sigma.sq*B.inv)
            tau.sq <- alpha*sigma.sq
            out$p.beta.theta.samples[i,] <- c(beta, sigma.sq, tau.sq)
        }
        
        colnames(out$p.beta.theta.samples) <- c(colnames(out$beta.hat), "sigma.sq", "tau.sq")
        
        
        if(fit.rep){
            if(cov.model == "matern"){
                theta <- as.matrix(as.data.frame(lapply(out$theta.alpha[,c("phi","nu")], rep, n.samples)))
            }else{
                theta <- as.matrix(as.data.frame(lapply(out$theta.alpha[,"phi"], rep, n.samples)))
            }
            
            theta <- t(cbind(out$p.beta.theta.samples[,c("sigma.sq","tau.sq")], theta))

            storage.mode(theta) <- "double"
            storage.mode(n.samples) <- "integer"

            n.report <- 100
            if("n.report" %in% names(elip.args)){
                n.report <-  elip.args[["n.report"]]
            }
            
            storage.mode(n.report) <- "integer"
            
            if(slgp){
                q <- p+r
                beta <- out$p.beta.theta.samples[,1:q]

                #out$X.str <- matrix(out$X.str, nrow=n)

                out$y.hat.samples <- (out$X.str%*%t(beta))[order(ord),,drop=FALSE]
                out$y.hat.quants <- t(apply(out$y.hat.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))

                #not yet implemented
                #storage.mode(out$X.str) <- "double"
                #storage.mode(beta) <- "double"
                #storage.mode(q) <- "integer"
                                
                #out$y.rep.samples <- .Call("rNNGPReplicated", out$X.str, q, n, n.neighbors, coords, cov.model.indx, nn.indx, nn.indx.lu, beta, theta,
                #                           n.samples, n.omp.threads, verbose, n.report)$y.rep.samples[order(ord),,drop=FALSE]
                
                ##add this just so it works with other support functions written for NNGP response which, at this point, is only spDiag
                out$sub.sample <- list(start=1, end=n.samples, thin=1)
                out$s.indx <- 1:n.samples
                
            }else{
                beta <- out$p.beta.theta.samples[,1:p]

                out$y.hat.samples <- (X%*%t(beta))[order(ord),,drop=FALSE]
                out$y.hat.quants <- t(apply(out$y.hat.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))
                
                storage.mode(beta) <- "double"
                
                out$y.rep.samples <- .Call("rNNGPReplicated", X, p, n, n.neighbors, coords, cov.model.indx, nn.indx, nn.indx.lu, beta, theta,
                                           n.samples, n.omp.threads, verbose, n.report)$y.rep.samples[order(ord),,drop=FALSE]

                out$y.rep.quants <- t(apply(out$y.rep.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))

                ##add this just so it works with other support functions written for NNGP response which, at this point, is only spDiag
                out$sub.sample <- list(start=1, end=n.samples, thin=1)
                out$s.indx <- 1:n.samples
            }

        }
    }
    
    
    ##put everthing back in the original order
    out$coords <- coords[order(ord),]
    out$y <- y[order(ord)]
    out$X <- X[order(ord),,drop=FALSE]

    if("X.str" %in% names(out)){
        out$X.str <- out$X.str[order(ord),]
    }

    out$type <- c("conjugate", "gaussian")
    
    class(out) <- "NNGP"
    
    if(slgp){
        out$knots <- knots
        out$type <- c(out$type, "SLGP")
    }
    
    out
}
