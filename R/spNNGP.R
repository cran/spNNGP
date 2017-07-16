spNNGP <- function(formula, data = parent.frame(), coords, method = "response", n.neighbors = 15, 
                   starting, tuning, priors, cov.model = "exponential",
                   n.samples, n.omp.threads = 1, search.type = "tree", return.neighbors = FALSE, verbose = TRUE, n.report=100, ...){
    
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
        stop("error: either the coords have more than two columns or then number of rows is different than
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
    ##NNGP method
    ####################################################
    method.names <- c("response","sequential")
    
    if(!method%in%method.names)
    {stop("error: specified method '",method,"' is not a valid option; choose from ", paste(method.names, collapse=", ", sep="") ,".")}
        
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
    ##Priors
    ####################################################
    sigma.sq.IG <- 0
    tau.sq.IG <- 0
    nu.Unif <- 0
    phi.Unif <- 0
    
    if(missing(priors)) {stop("error: prior list for the parameters must be specified")}
    
    names(priors) <- tolower(names(priors))

    if(!"sigma.sq.ig" %in% names(priors)){stop("error: sigma.sq.IG must be specified")}
    sigma.sq.IG <- priors[["sigma.sq.ig"]]
    
    if(!is.vector(sigma.sq.IG) || length(sigma.sq.IG) != 2){stop("error: sigma.sq.IG must be a vector of length 2")}
    if(any(sigma.sq.IG <= 0)){stop("error: sigma.sq.IG must be a positive vector of length 2")}

    if(!"tau.sq.ig" %in% names(priors)){stop("error: tau.sq.IG must be specified")}
    tau.sq.IG <- priors[["tau.sq.ig"]]
    
    if(!is.vector(tau.sq.IG) || length(tau.sq.IG) != 2){stop("error: tau.sq.IG must be a vector of length 2")}
    if(any(tau.sq.IG <= 0)){stop("error: tau.sq.IG must be a positive vector of length 2")}

    if(!"phi.unif" %in% names(priors)){stop("error: phi.Unif must be specified")}
    phi.Unif <- priors[["phi.unif"]]
    
    if(!is.vector(phi.Unif) || length(phi.Unif) != 2){stop("error: phi.Unif must be a vector of length 2")}
    if(any(phi.Unif <= 0, phi.Unif[1] >= phi.Unif[2])){stop("error: phi.Unif must be a positive vector of length 2 with element 1 < element 2")}
    
    if(cov.model == "matern"){
        
        if(!"nu.unif" %in% names(priors)){stop("error: nu.Unif must be specified")}
        nu.Unif <- priors[["nu.unif"]]
        
        if(!is.vector(nu.Unif) || length(nu.Unif) != 2){stop("error: nu.Unif must be a vector of length 2")}
        if(any(nu.Unif <= 0, nu.Unif[1] >= nu.Unif[2])){stop("error: nu.Unif must be a positive vector of length 2 with element 1 < element 2")}
    }
    
    storage.mode(sigma.sq.IG) <- "double"
    storage.mode(tau.sq.IG) <- "double"
    storage.mode(phi.Unif) <- "double"
    storage.mode(nu.Unif) <- "double"
   
    ####################################################
    ##Starting values
    ####################################################
    beta.starting <- 0
    sigma.sq.starting <- 0
    tau.sq.starting <- 0
    phi.starting <- 0
    nu.starting <- 0
    
    if(missing(starting)){stop("error: starting value list for the parameters must be specified")}
    
    names(starting) <- tolower(names(starting))   
    
    if(!"beta" %in% names(starting)){
        beta.starting <- as.vector(coefficients(lm(y~X-1)))
    }else{
        beta.starting <- starting[["beta"]]
        if(length(beta.starting) != p){stop("error: supplied beta.starting is the wrong length")}
    }
    
    if(!"sigma.sq" %in% names(starting)){stop("error: sigma.sq must be specified in starting value list")}
    sigma.sq.starting <- starting[["sigma.sq"]][1]
    
    if(!"tau.sq" %in% names(starting)){stop("error: tau.sq must be specified in starting value list")}
    tau.sq.starting <- starting[["tau.sq"]][1]
    
    if(!"phi" %in% names(starting)){stop("error: phi must be specified in starting value list")}
    phi.starting <- starting[["phi"]][1]
    
    if(cov.model == "matern"){
        if(!"nu" %in% names(starting)){stop("error: nu must be specified in starting value list")}
        nu.starting <- starting[["nu"]][1]
    }
    
    storage.mode(beta.starting) <- "double"
    storage.mode(sigma.sq.starting) <- "double"
    storage.mode(tau.sq.starting) <- "double"
    storage.mode(phi.starting) <- "double"
    storage.mode(nu.starting) <- "double"
    
    ####################################################
    ##Tuning values
    ####################################################
    sigma.sq.tuning <- 0
    tau.sq.tuning <- 0
    phi.tuning <- 0
    nu.tuning <- 0
    
    if(missing(tuning)){stop("error: tuning value vector for the spatial parameters must be specified")}
    
    names(tuning) <- tolower(names(tuning))
    
    if(!"sigma.sq" %in% names(tuning) & method == "response"){stop("error: sigma.sq must be specified in tuning value list")}
    sigma.sq.tuning <- tuning[["sigma.sq"]][1]
    
    if(!"tau.sq" %in% names(tuning) & method == "response"){stop("error: tau.sq must be specified in tuning value list")}
    tau.sq.tuning <- tuning[["tau.sq"]][1]
    
    if(!"phi" %in% names(tuning)){stop("error: phi must be specified in tuning value list")}
    phi.tuning <- tuning[["phi"]][1]
    
    if(cov.model == "matern"){
        if(!"nu" %in% names(tuning)){stop("error: nu must be specified in tuning value list")}
        nu.tuning <- tuning[["nu"]][1]
    }    

    storage.mode(sigma.sq.tuning) <- "double"
    storage.mode(tau.sq.tuning) <- "double"
    storage.mode(phi.tuning) <- "double"
    storage.mode(nu.tuning) <- "double"

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
    storage.mode(n.samples) <- "integer"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(n.neighbors) <- "integer"
    storage.mode(n.report) <- "integer"
    storage.mode(verbose) <- "integer"

    ####################################################
    ##Pack it up and off it goes
    ####################################################

    ptm <- proc.time()

    nngp <- paste(substr(method,1,1),"NNGP",collapse="",sep="")
    
    out <- .Call(nngp, y, X, p, n, n.neighbors, coords, cov.model.indx,
                 sigma.sq.IG, tau.sq.IG, phi.Unif, nu.Unif, 
                 beta.starting, sigma.sq.starting, tau.sq.starting, phi.starting, nu.starting,
                 sigma.sq.tuning, tau.sq.tuning, phi.tuning, nu.tuning,
                 n.samples, search.type.indx, return.neighbors, n.omp.threads, verbose, n.report)
        
    out$run.time <- proc.time() - ptm

    out$p.theta.samples <- mcmc(t(out$p.theta.samples))
    out$p.beta.samples <- mcmc(t(out$p.beta.samples))

    col.names <- c("sigma.sq", "tau.sq", "phi")
    
    if(cov.model == "matern"){
        col.names <- c("sigma.sq", "tau.sq", "phi", "nu")
    }
    
    colnames(out$p.theta.samples) <- col.names
    colnames(out$p.beta.samples) <- x.names

    out$y <- y
    out$X <- X
    out$n.neighbors <- n.neighbors
    out$coords <- coords
    out$ord <- ord
    out$cov.model <- cov.model
    out$cov.model.indx <- cov.model.indx

    if(return.neighbors){
        out$n.indx <- mk.n.indx.list(out$n.indx, n, n.neighbors)
    }

    class(out) <- nngp
    
    out
    
}
