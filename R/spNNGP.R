spNNGP <- function(formula, data = parent.frame(), coords, method = "response", family="gaussian", weights, n.neighbors = 15, 
                   starting, tuning, priors, cov.model = "exponential",
                   n.samples, n.omp.threads = 1, search.type = "cb", ord, 
                   return.neighbor.info = FALSE, neighbor.info,
                   fit.rep = FALSE, sub.sample, verbose = TRUE, n.report=100, ...){
    
    ####################################################
    ##Check for unused args
    ####################################################
    formal.args <- c(names(formals(sys.function(sys.parent()))), "u.search.type")
    
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
    ##Family
    #################################################### 
    family.names <- c("gaussian","binomial")
    family <- tolower(family)
    
    if(!family%in%family.names){stop("error: specified family '",family,"' is not a valid option; choose from ", paste(family.names, collapse=", ", sep="") ,".")}    
    if(method == "response" & family != "gaussian"){stop("error: only the latent method can be used with non-Gaussain family.")}
    
    if(family == "binomial"){
        if(missing(weights)){
            weights <- rep(1, n)
        }else if(length(weights) == 1){
            weights <- rep(weights, n)
        }else{
            if(length(weights) != n){
                stop("error: for non-Gaussian models weights must be of length n or 1 (in which case the specified value is repeted n times)")
            }
        }
    }else{
        weights <- NA
    }
    
    ####################################################
    ##Neighbors and ordering
    ####################################################
    neighbor.info.provided <- FALSE
    
    u.search.type <- 2 ##2 is very fast, 1 is slightly faster than 0, and 0 is the orginal slow one (2 and 0 should match, 1 is also corrected just different opposite sorting among the children)
    if("u.search.type" %in% names(elip.args)){
        u.search.type <- elip.args[["u.search.type"]]
    }
    
    if(!missing(neighbor.info)){

        warning("Using user defined neighbor.info. No checks are done on the supplied neighbor information.")
        
        if(!all(c("n.neighbors","nn.indx","nn.indx.lu","ord") %in% names(neighbor.info))){stop("The supplied neighbor.info is malformed.")}
        
        nn.indx <- neighbor.info$nn.indx
        nn.indx.lu <- neighbor.info$nn.indx.lu
        ord <- neighbor.info$ord
        n.neighbors <- neighbor.info$n.neighbors
        nn.indx.run.time <- neighbor.info$nn.indx.run.time
        neighbor.info.provided <- TRUE

        if(method == "latent"){
            
            if(!all(c("u.indx", "u.indx.lu", "ui.indx") %in% names(neighbor.info))){
                ##must be neighbor.info from a response or conjugate model
                neighbor.info <- c(neighbor.info, mkUIndx(n, n.neighbors, nn.indx, nn.indx.lu, u.search.type))
                print("Computing additional index needed for method = 'latent' using user defined neighbor.info. This might take a while if n is large.")
            }

            u.indx <- neighbor.info$u.indx
            u.indx.lu <- neighbor.info$u.indx.lu
            ui.indx <- neighbor.info$ui.indx
            u.indx.run.time <- neighbor.info$u.indx.run.time
        }
        
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

    if(family == "binomial"){
        weights <- weights[ord]
    }
    
    storage.mode(y) <- "double"
    storage.mode(X) <- "double"
    storage.mode(p) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(coords) <- "double"
    storage.mode(weights) <- "integer"
    
    ####################################################
    ##NNGP method
    ####################################################
    method.names <- c("response","latent")
    method <- tolower(method)
    
    if(!method%in%method.names){stop("error: specified method '",method,"' is not a valid option; choose from ", paste(method.names, collapse=", ", sep="") ,".")}

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

    if(family == "gaussian"){
        if(!"tau.sq.ig" %in% names(priors)){stop("error: tau.sq.IG must be specified")}
        tau.sq.IG <- priors[["tau.sq.ig"]]
   
        if(!is.vector(tau.sq.IG) || length(tau.sq.IG) != 2){stop("error: tau.sq.IG must be a vector of length 2")}
        if(any(tau.sq.IG <= 0)){stop("error: tau.sq.IG must be a positive vector of length 2")}
    }
    
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

    if(family == "gaussian"){
        if(!"tau.sq" %in% names(starting)){stop("error: tau.sq must be specified in starting value list")}
        tau.sq.starting <- starting[["tau.sq"]][1]
    }
    
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

    if(family == "gaussian"){
        if(!"tau.sq" %in% names(tuning) & method == "response"){stop("error: tau.sq must be specified in tuning value list")}
        tau.sq.tuning <- tuning[["tau.sq"]][1]
    }
        
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
    if(!neighbor.info.provided){

        if(verbose){
            cat("----------------------------------------\n");
            cat("\tBuilding the neighbor list\n");
            cat("----------------------------------------\n");
        }
        
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
        
        storage.mode(nn.indx) <- "integer"
        storage.mode(nn.indx.lu) <- "integer"
        
        if(method == "latent"){

            if(verbose){
                cat("----------------------------------------\n");
                cat("Building the neighbors of neighbors list\n");
                cat("----------------------------------------\n");
            }
            
            indx <- mkUIndx(n, n.neighbors, nn.indx, nn.indx.lu, u.search.type)
            
            u.indx <- indx$u.indx
            u.indx.lu <- indx$u.indx.lu
            ui.indx <- indx$ui.indx
            u.indx.run.time <- indx$run.time
            
            storage.mode(u.indx) <- "integer"
            storage.mode(u.indx.lu) <- "integer"
            storage.mode(ui.indx) <- "integer"
        }
    }  

    ####################################################
    ##fitted and replicated y 
    ####################################################
    n.rep <- 0
    rep.indx <- 0
    
    if(fit.rep){
        
        if(missing(sub.sample)){
            sub.sample <- list()
        }

        start <- ifelse(!"start" %in% names(sub.sample), floor(0.5*n.samples), sub.sample$start)
        end <- ifelse(!"end" %in% names(sub.sample), n.samples, sub.sample$end)
        thin <- ifelse(!"thin" %in% names(sub.sample), 1, sub.sample$thin)   
        if(!is.numeric(start) || start >= n.samples){stop("error: in sub.sample, invalid start")}
        if(!is.numeric(end) || end > n.samples){stop("error: in sub.sample, invalid end")}
        if(!is.numeric(thin) || thin >= n.samples){stop("error: in sub.sample, invalid thin")}
        
        s.indx <- seq(as.integer(start), as.integer(end), by=as.integer(thin))
        sub.sample <- list(start=start, end=end, thin=thin)
        n.rep <- length(s.indx)
        rep.indx <- rep(0, n.samples)
        rep.indx[1:n.samples %in% s.indx] <- 1
    }
    
    storage.mode(n.rep) <- "integer"
    storage.mode(rep.indx) <- "integer"
    
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

    if(family == "gaussian"){

        if(method == "response"){
            
            out <- .Call("rNNGP", y, X, p, n, n.neighbors, coords, cov.model.indx, nn.indx, nn.indx.lu, 
                         sigma.sq.IG, tau.sq.IG, phi.Unif, nu.Unif, 
                         beta.starting, sigma.sq.starting, tau.sq.starting, phi.starting, nu.starting,
                         sigma.sq.tuning, tau.sq.tuning, phi.tuning, nu.tuning,
                         n.samples, n.omp.threads, verbose, n.report, n.rep, rep.indx)
            
        }else{##sequential
            
            out <- .Call("sNNGP", y, X, p, n, n.neighbors, coords, cov.model.indx, nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx,
                         sigma.sq.IG, tau.sq.IG, phi.Unif, nu.Unif, 
                         beta.starting, sigma.sq.starting, tau.sq.starting, phi.starting, nu.starting,
                         sigma.sq.tuning, tau.sq.tuning, phi.tuning, nu.tuning,
                         n.samples, n.omp.threads, verbose, n.report)
        }
            
        col.names <- c("sigma.sq", "tau.sq", "phi")
        
        if(cov.model == "matern"){
            col.names <- c("sigma.sq", "tau.sq", "phi", "nu")
        }
        
    }else{
        
        out <- .Call("sNNGPLogit", y, X, p, n, n.neighbors, coords, weights, cov.model.indx, nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx,
                     sigma.sq.IG, phi.Unif, nu.Unif, 
                     beta.starting, sigma.sq.starting, phi.starting, nu.starting,
                     sigma.sq.tuning, phi.tuning, nu.tuning,
                     n.samples, n.omp.threads, verbose, n.report)
        
        col.names <- c("sigma.sq", "phi")
        
        if(cov.model == "matern"){
            col.names <- c("sigma.sq", "phi", "nu")
        }        
    }

    out$run.time <- proc.time() - ptm

    out$p.theta.samples <- mcmc(t(out$p.theta.samples))
    out$p.beta.samples <- mcmc(t(out$p.beta.samples))
        
    colnames(out$p.theta.samples) <- col.names
    colnames(out$p.beta.samples) <- x.names
    
    if(return.neighbor.info){
        if(method == "response"){
            out$neighbor.info <- list(n.neighbors = n.neighbors, n.indx=mk.n.indx.list(nn.indx, n, n.neighbors),
                                      nn.indx=nn.indx, nn.indx.lu=nn.indx.lu, ord=ord,
                                      nn.indx.run.time=nn.indx.run.time)
        }else{
            out$neighbor.info <- list(n.neighbors = n.neighbors, n.indx=mk.n.indx.list(nn.indx, n, n.neighbors),
                                      nn.indx=nn.indx, nn.indx.lu=nn.indx.lu, u.indx=u.indx, u.indx.lu=u.indx.lu, ui.indx=ui.indx, ord=ord,
                                      nn.indx.run.time=nn.indx.run.time, u.indx.run.time=u.indx.run.time)
            
        }
    }

    ##put everthing back in the original order
    out$coords <- coords[order(ord),]
    out$y <- y[order(ord)]
    out$X <- X[order(ord),,drop=FALSE]
    out$weights <- weights[order(ord)]

    if(method == "latent"){
        out$p.w.samples <- out$p.w.samples[order(ord),,drop=FALSE]
    }

    ##finish off the replicate and fitted values
    if(fit.rep){
        if(method == "response"){
            out$y.hat.samples <- out$X%*%t(as.matrix(out$p.beta.samples)[s.indx,,drop=FALSE])
            out$y.hat.quants <- t(apply(out$y.hat.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))           
            out$y.rep.samples <- out$y.rep.samples[order(ord),,drop=FALSE]
            out$y.rep.quants <- t(apply(out$y.rep.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))
        }else{##sequential
            out$y.hat.samples <- out$X%*%t(as.matrix(out$p.beta.samples)[s.indx,,drop=FALSE]) + out$p.w.samples[,s.indx,drop=FALSE]
            out$y.hat.quants <- t(apply(out$y.hat.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))

            if(family == "gaussian"){
                out$y.rep.samples <- out$X%*%t(as.matrix(out$p.beta.samples)[s.indx,,drop=FALSE]) + out$p.w.samples[,s.indx,drop=FALSE] + sapply(out$p.theta.samples[s.indx,"tau.sq"], function(x) sqrt(x)*rnorm(n))
                out$y.rep.quants <- t(apply(out$y.rep.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))
            }else{##binomial
                out$y.rep.samples <- apply(1/(1+exp(-out$y.hat.samples)), 2, function(x) rbinom(n, out$weights, x))
                out$y.rep.quants <- t(apply(out$y.rep.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))
            }
            
        }
        out$sub.sample <- sub.sample
        out$s.indx <- s.indx
    }
    
    out$n.neighbors <- n.neighbors
    out$cov.model <- cov.model
    out$cov.model.indx <- cov.model.indx
    out$call <- cl
    out$starting <- starting
    out$priors <- priors
    out$tuning <- tuning
    out$type <- c(method, family)
    class(out) <- "NNGP"

    out
    
}
