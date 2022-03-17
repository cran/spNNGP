PGLogit <- function(formula, weights = 1, data = parent.frame(), n.samples, n.omp.threads = 1, fit.rep = FALSE, sub.sample, verbose = TRUE, ...){
    
    ####################################################
    ##Check for unused args
    ####################################################
    formal.args <- names(formals(sys.function(sys.parent())))
    elip.args <- names(list(...))
    for(i in elip.args){
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
    
    if(length(weights) == 1){
        weights <- rep(weights, n)
    }else{
        if(length(weights) != n){
            stop("error: for non-Gaussian models weights must be of length n or 1 (in which case the specified value is repeted n times)")
        }
    }

    beta.starting <- coefficients(glm((y/weights)~X-1, weights=weights, family="binomial"))
    
    storage.mode(y) <- "double"
    storage.mode(X) <- "double"
    storage.mode(p) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(weights) <- "integer"
    storage.mode(beta.starting) <- "double"

    ####################################################
    ##fitted and replicated y 
    ####################################################
 
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

        sub.sample <- list(start=start, end=end, thin=thin)
        s.indx <- seq(as.integer(start), as.integer(end), by=as.integer(thin))
    }

    ####################################################
    ##Other stuff
    ####################################################
    storage.mode(n.samples) <- "integer"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"

    ptm <- proc.time()
    
    out <- .Call("PGLogit", y, X, p, n, weights, beta.starting, n.samples, n.omp.threads, verbose)

    out$run.time <- proc.time() - ptm

    out$p.beta.samples <- mcmc(t(out$p.beta.samples))
    colnames(out$p.beta.samples) <- x.names
    out$X <- X
    out$y <- y
    out$weights <- weights
    out$call <- cl

    ##finish off the replicate and fitted values
    if(fit.rep){

        ##fitted
        out$y.hat.samples <- out$X%*%t(as.matrix(out$p.beta.samples)[s.indx,,drop=FALSE])
        out$y.hat.quantiles <- t(apply(out$y.hat.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))

        ##replicates
        prob.samples <- 1/(1+exp(-out$y.hat.samples))
        out$y.rep.samples <- t(apply(prob.samples, 1, function(x) rbinom(n, out$weights, x)))
        out$y.rep.quants <- t(apply(out$y.rep.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))
        out$sub.sample <- sub.sample
        out$s.indx <- s.indx
    }
    
    class(out) <- "PGLogit"
    
    out
}
