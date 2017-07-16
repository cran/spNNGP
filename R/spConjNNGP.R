spConjNNGP <- function(formula, data = parent.frame(), coords, n.neighbors = 15,
                       theta.alpha, sigma.sq.IG, cov.model = "exponential",
                       k.fold, score.rule,
                       X.0, coords.0, 
                       n.omp.threads = 1, search.type = "tree",
                       return.neighbors = FALSE, verbose=TRUE, ...){
    
    ####################################################
    ##Check for unused args
    ####################################################
    formal.args <- names(formals(sys.function(sys.parent())))
    elip.args <- names(list(...))
    for(i in elip.args){
        if(! i %in% formal.args)
            warning("'",i, "' is not an argument")
    }

    ####################################################
    ##Formula
    ####################################################
    if(missing(formula)){stop("error: formula must be specified")}
    
    if(class(formula) == "formula"){
        
        holder <- parseFormula(formula, data)
        y <- holder[[1]]
        X <- as.matrix(holder[[2]])
        x.names <- holder[[3]]
        
    }else{
        stop("error: formula is misspecified")
    }
    
    p <- ncol(X)
    n <- nrow(X)

    ##Coords and ordering
    if(!is.matrix(coords)){stop("error: coords must n-by-2 matrix of xy-coordinate locations")}
    if(ncol(coords) != 2 || nrow(coords) != n){
        stop("error: either coords has more than two columns or the number of rows is different than
          data used in the model formula")
    }

    ##order by x
    ord <- order(coords[,1])
    coords <- coords[ord,]
    X <- X[ord,,drop=FALSE]
    y <- y[ord]

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
        if(!is.matrix(X.0) || ncol(X.0) != p){stop(paste("error: coords must n.0-by-",p," matrix"), sep="")}
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
    if(is.vector(theta.alpha)){
        theta.alpha <- t(as.matrix(theta.alpha))
    }

    if(!is.matrix(theta.alpha)){stop("error: theta.alpha must matrix of parameter values")}
   
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
    search.type.names <- c("brute", "tree")
    
    if(!search.type %in% search.type.names){
        stop("error: specified search.type '",search.type,"' is not a valid option; choose from ", paste(search.type.names, collapse=", ", sep="") ,".")
    }
    
    search.type.indx <- which(search.type == search.type.names)-1##obo for cov model lookup on c side
    storage.mode(search.type.indx) <- "integer"
    storage.mode(return.neighbors) <- "integer"
    
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
    
    if(!missing(k.fold)){
        
        if(k.fold < 2 || k.fold >= n){
            stop("error: k.fold must be greater than 1 and less than n")
        }
        
        if(missing(score.rule)){
            warning("score.rule not specified, defaulting to CRPS")
            score.rule <- "crps"
        }

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

            fold.verbose <- 0
            storage.mode(fold.verbose) <- "integer"

            fold.return.neighbors <- 0
            storage.mode(fold.return.neighbors) <- "integer"
            
            out <- .Call("cNNGP", y.mod, X.mod, p, n.mod, coords.mod, theta.alpha,
                         X.ho, coords.ho, n.ho, nn.indx.mod.ho, g,
                         n.neighbors, sigma.sq.IG, cov.model.indx, n.omp.threads, search.type.indx, fold.return.neighbors, fold.verbose)
            
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
        g <- 1
        
        storage.mode(theta.alpha) <- "double"
        storage.mode(g) <- "integer"
    }


    ##final run
    out <- .Call("cNNGP", y, X, p, n, coords, theta.alpha,
                 X.0, coords.0, n.0, nn.indx.0, g,
                 n.neighbors, sigma.sq.IG, cov.model.indx, n.omp.threads, search.type.indx, return.neighbors, verbose)
        
    out$run.time <- proc.time() - ptm

    sigma.sq.hat <- out$ab[2,]/(out$ab[1,]-1)
    out$theta.alpha.sigmaSq <- cbind(t(theta.alpha), sigma.sq.hat)

    out$beta.hat <- t(out$beta)
    colnames(out$beta.hat) <- x.names

    out$y <- y
    out$X <- X
    out$n.neighbors <- n.neighbors
    out$coords <- coords
    out$ord <- ord
    out$cov.model <- cov.model

    if(!missing(k.fold)){
        out$k.fold.scores <- k.fold.scores
    }

    if(return.neighbors){
        out$n.indx <- mk.n.indx.list(out$n.indx, n, n.neighbors)
    }

    class(out) <- "cNNGP"
        
    out
}
