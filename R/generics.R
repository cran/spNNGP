fitted.NNGP <- function(object, sub.sample, ...){

    elip.args <- list(...)

    n.omp.threads <- 1
    if("n.omp.threads" %in% names(elip.args)){
        n.omp.threads <-  elip.args[["n.omp.threads"]]
    }
    
    verbose <- 1
    if("verbose" %in% names(elip.args)){
        verbose <-  elip.args[["verbose"]]
    }
    
    n.report <- 100
    if("n.report" %in% names(elip.args)){
        n.report <-  elip.args[["n.report"]]
    }
    
    if(class(object)[2] == "conjugate"){

        if(!missing(sub.sample)){
            if(!(is.atomic(sub.sample) && length(sub.sample) == 1L && sub.sample > 0)){
                stop("error: for an object of class conjugate, the sub.sample argument should be in integer value indicating the number of fitted and replicated samples to collect.")
            }
        }
                       
        if("SLGP" %in% class(object)){

            if("y.hat.samples" %in% names(object)){
                if(missing(sub.sample) || (!missing(sub.sample) && sub.sample == ncol(object$y.hat.samples))){
                    return(list(y.hat.quants=object$y.hat.quants, y.hat.samples=object$y.hat.samples))
                }
            }

            ##y.rep.samples is missing or requesting a different number of samples
            if(missing(sub.sample)){
                stop("error: for an object of class conjugate created with fit.rep=FALSE, the sub.sample argument should be provided to indicate the number of fitted and replicated samples to collect.")
            }
            
            out <- list()
            
            ##get samples
            b <- object$bb
            p <- length(b)
            B.inv <- matrix(object$bB.inv, length(b), length(b))
            B.inv[upper.tri(B.inv, FALSE)] <- t(B.inv)[upper.tri(B.inv, FALSE)]
            alpha <- object$theta.alpha[,"alpha"]
            out$p.beta.theta.samples <- matrix(0, nrow=sub.sample, ncol=length(b)+2)
            
            for(i in 1:sub.sample){
                sigma.sq <- rigamma(1, object$ab[1,], object$ab[2,])
                beta <- rmvn(1, B.inv%*%b, sigma.sq*B.inv)
                tau.sq <- alpha*sigma.sq
                out$p.beta.theta.samples[i,] <- c(beta, sigma.sq, tau.sq)
            }
            
            colnames(out$p.beta.theta.samples) <- c(colnames(object$beta.hat), "sigma.sq", "tau.sq")
 
            beta <- out$p.beta.theta.samples[,1:ncol(object$X.str)]

            out$y.hat.samples <- object$X.str%*%t(beta)
            out$y.hat.quants <- t(apply(out$y.hat.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))

            ##add this just so it works with other support functions written for NNGP response which, at this point, is only spDiag
            out$sub.sample <- list(start=1, end=sub.sample, thin=1)
            out$s.indx <- 1:sub.sample
            
        }else{##NNGP
            
            if("y.rep.samples" %in% names(object)){
                if(missing(sub.sample) || (!missing(sub.sample) && sub.sample == ncol(object$y.rep.samples))){
                    return(list(y.hat.quants=object$y.hat.quants, y.hat.samples=object$y.hat.samples,
                                y.rep.quants=object$y.rep.quants, y.rep.samples=object$y.rep.samples))
                }
            }
            
            ##y.rep.samples is missing or requesting a different number of samples
            if(missing(sub.sample)){
                stop("error: for an object of class conjugate created with fit.rep=FALSE, the sub.sample argument should be provided to indicate the number of fitted and replicated samples to collect.")
            }

            if(!"neighbor.info" %in% names(object)){
                stop("error: the neighbor.info is not in the model object, call spConjNNGP again with return.neighbor.info = TRUE.")
            }
            
            out <- list()
            
            ##get samples
            b <- object$bb
            p <- length(b)
            B.inv <- matrix(object$bB.inv, length(b), length(b))
            B.inv[upper.tri(B.inv, FALSE)] <- t(B.inv)[upper.tri(B.inv, FALSE)]
            alpha <- object$theta.alpha[,"alpha"]
            out$p.beta.theta.samples <- matrix(0, nrow=sub.sample, ncol=length(b)+2)
            
            for(i in 1:sub.sample){
                sigma.sq <- rigamma(1, object$ab[1,], object$ab[2,])
                beta <- rmvn(1, B.inv%*%b, sigma.sq*B.inv)
                tau.sq <- alpha*sigma.sq
                out$p.beta.theta.samples[i,] <- c(beta, sigma.sq, tau.sq)
            }
            
            colnames(out$p.beta.theta.samples) <- c(colnames(object$beta.hat), "sigma.sq", "tau.sq")
            
            if(object$cov.model == "matern"){
                theta <- as.matrix(as.data.frame(lapply(object$theta.alpha[,c("phi","nu")], rep, sub.sample)))
            }else{
                theta <- as.matrix(as.data.frame(lapply(object$theta.alpha[,"phi"], rep, sub.sample)))
            }
            
            theta <- t(cbind(out$p.beta.theta.samples[,c("sigma.sq","tau.sq")], theta))                   
            beta <- out$p.beta.theta.samples[,1:p]
            
            ##get fitted and replicated
            X <- object$X
            out$y.hat.samples <- X%*%t(beta)
            out$y.hat.quants <- t(apply(out$y.hat.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))
            
            n <- nrow(X)
            n.neighbors <- object$n.neighbors
            cov.model.indx <- object$cov.model.indx
            nn.indx <- object$neighbor.info$nn.indx
            nn.indx.lu <- object$neighbor.info$nn.indx.lu
            ord <- object$neighbor.info$ord
            
            coords <- object$coords[ord,]
            X <- X[ord,,drop=FALSE]
            
            storage.mode(X) <- "double"
            storage.mode(p) <- "integer"
            storage.mode(n) <- "integer"
            storage.mode(n.neighbors) <- "integer"
            storage.mode(coords) <- "double"
            storage.mode(cov.model.indx) <- "integer"
            storage.mode(nn.indx) <- "integer"
            storage.mode(nn.indx.lu) <- "integer"
            storage.mode(beta) <- "double"
            storage.mode(theta) <- "double"
            storage.mode(sub.sample) <- "integer" 
            storage.mode(n.omp.threads) <- "integer"
            storage.mode(verbose) <- "integer"
            storage.mode(n.report) <- "integer"
            
            out$y.rep.samples <- .Call("rNNGPReplicated", X, p, n, n.neighbors, coords, cov.model.indx, nn.indx, nn.indx.lu, beta, theta,
                                       sub.sample, n.omp.threads, verbose, n.report)$y.rep.samples[order(ord),,drop=FALSE]
            
            out$y.rep.quants <- t(apply(out$y.rep.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))

            ##add this just so it works with other support functions written for NNGP response which, at this point, is only spDiag
            out$sub.sample <- list(start=1, end=sub.sample, thin=1)
            out$s.indx <- 1:sub.sample
                        
        }
        
        return(out)
        
    }else{##not a conjugate class
        
        update.fit.rep <- FALSE
        
        if(!missing(sub.sample)){
            update.fit.rep <- TRUE
        }
        
        if(missing(sub.sample) & !"s.indx" %in% names(object)){
            update.fit.rep <- TRUE
            sub.sample <- list() #uses default sub.sample start, end, thin
        }        
        
        if(update.fit.rep){
            n.samples <- nrow(object$p.beta.samples)
            start <- ifelse(!"start" %in% names(sub.sample), floor(0.5*n.samples), sub.sample$start)
            end <- ifelse(!"end" %in% names(sub.sample), n.samples, sub.sample$end)
            thin <- ifelse(!"thin" %in% names(sub.sample), 1, sub.sample$thin)
            if(!is.numeric(start) || start >= n.samples){stop("invalid start")}
            if(!is.numeric(end) || end > n.samples){stop("invalid end")}
            if(!is.numeric(thin) || thin >= n.samples){stop("invalid thin")}
            sub.sample <- list(start=start, end=end, thin=thin)
            s.indx <- seq(as.integer(start), as.integer(end), by=as.integer(thin))
            n.samples <- length(s.indx)
        }
        
        if(class(object)[2] == "latent"){
            
            if(update.fit.rep){
                beta <- as.matrix(object$p.beta.samples)[s.indx,,drop=FALSE]
                w <- object$p.w.samples[,s.indx,drop=FALSE]
                
                X <- object$X
                y.hat.samples <- X%*%t(beta) + w
                y.hat.quants <- t(apply(y.hat.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))
                
                if(class(object)[2] == "gaussian"){
                    tau.sq <- object$p.theta.samples[s.indx,"tau.sq"]
                    y.rep.samples <- y.hat.samples + sapply(tau.sq, function(x) sqrt(x)*rnorm(nrow(X)))
                }else{
                    y.rep.samples <- apply(1/(1+exp(-y.hat.samples)), 2, function(x) rbinom(nrow(X), object$weights, x))
                }
                
                y.rep.quants <- t(apply(y.rep.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))

                return(list(y.hat.quants=y.hat.quants, y.hat.samples=y.hat.samples,
                            y.rep.quants=y.rep.quants, y.rep.samples=y.rep.samples, sub.sample=sub.sample, s.indx=s.indx))     
            }else{                
                return(list(y.hat.quants=object$y.hat.quants, y.hat.samples=object$y.hat.samples,
                            y.rep.quants=object$y.rep.quants, y.rep.samples=object$y.rep.samples, sub.sample=object$sub.sample, s.indx=object$s.indx))
            }
            
            
        }else{##response
            
            if(update.fit.rep){

                if(!"neighbor.info" %in% names(object)){
                    stop("error: the neighbor.info is not in the model object, call spConjNNGP again with return.neighbor.info = TRUE.")
                }
                
                beta <- as.matrix(object$p.beta.samples)[s.indx,,drop=FALSE]
                
                X <- object$X
                y.hat.samples <- X%*%t(beta)
                y.hat.quants <- t(apply(y.hat.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))
                
                theta <- t(as.matrix(object$p.theta.samples)[s.indx,])
                n.samples <- ncol(theta)
                
                p <- ncol(X)
                n <- nrow(X)
                
                n.neighbors <- object$n.neighbors
                cov.model.indx <- object$cov.model.indx
                nn.indx <- object$neighbor.info$nn.indx
                nn.indx.lu <- object$neighbor.info$nn.indx.lu
                ord <- object$neighbor.info$ord
                
                coords <- object$coords[ord,]
                X <- X[ord,,drop=FALSE]
                                       
                storage.mode(X) <- "double"
                storage.mode(p) <- "integer"
                storage.mode(n) <- "integer"
                storage.mode(n.neighbors) <- "integer"
                storage.mode(coords) <- "double"
                storage.mode(cov.model.indx) <- "integer"
                storage.mode(nn.indx) <- "integer"
                storage.mode(nn.indx.lu) <- "integer"
                storage.mode(beta) <- "double"
                storage.mode(theta) <- "double"
                storage.mode(n.samples) <- "integer" 
                storage.mode(n.omp.threads) <- "integer"
                storage.mode(verbose) <- "integer"
                storage.mode(n.report) <- "integer"
                
                y.rep.samples <- .Call("rNNGPReplicated", X, p, n, n.neighbors, coords, cov.model.indx, nn.indx, nn.indx.lu, beta, theta,
                                       n.samples, n.omp.threads, verbose, n.report)$y.rep.samples[order(ord),,drop=FALSE]
                
                y.rep.quants <- t(apply(y.rep.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))
                
                return(list(y.hat.quants=y.hat.quants, y.hat.samples=y.hat.samples,
                            y.rep.quants=y.rep.quants, y.rep.samples=y.rep.samples, sub.sample=sub.sample, s.indx=s.indx))
            }else{
                
                return(list(y.hat.quants=object$y.hat.quants, y.hat.samples=object$y.hat.samples,
                            y.rep.quants=object$y.rep.quants, y.rep.samples=object$y.rep.samples, sub.sample=object$sub.sample, s.indx=object$s.indx))
            }
            
        }
    }
}

fitted.PGLogit <- function(object, sub.sample, ...){

    update.fit.rep <- FALSE
    
    if(!missing(sub.sample)){
        update.fit.rep <- TRUE
    }

    if(missing(sub.sample) & !"s.indx" %in% names(object)){
        update.fit.rep <- TRUE
        sub.sample <- list() #uses default sub.sample start, end, thin
    }
        
    if(update.fit.rep){
        n.samples <- nrow(object$p.beta.samples)
        start <- ifelse(!"start" %in% names(sub.sample), floor(0.5*n.samples), sub.sample$start)
        end <- ifelse(!"end" %in% names(sub.sample), n.samples, sub.sample$end)
        thin <- ifelse(!"thin" %in% names(sub.sample), 1, sub.sample$thin)
        if(!is.numeric(start) || start >= n.samples){stop("invalid start")}
        if(!is.numeric(end) || end > n.samples){stop("invalid end")}
        if(!is.numeric(thin) || thin >= n.samples){stop("invalid thin")}
        sub.sample <- list(start=start, end=end, thin=thin)
        s.indx <- seq(as.integer(start), as.integer(end), by=as.integer(thin))
        n.samples <- length(s.indx)
        
        ##fitted
        y.hat.samples <- object$X%*%t(as.matrix(object$p.beta.samples)[s.indx,,drop=FALSE])
        y.hat.quants <- t(apply(y.hat.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))
        
        ##replicates
        y.rep.samples <- apply(1/(1+exp(-y.hat.samples)), 2, function(x) rbinom(nrow(object$X), object$weights, x))
        y.rep.quants <- t(apply(y.rep.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))

        return(list(y.hat.quants=y.hat.quants, y.hat.samples=y.hat.samples, y.rep.quants=y.rep.quants, y.rep.samples=y.rep.samples, sub.sample=sub.sample, s.indx=s.indx))
    }else{
        return(list(y.hat.quants=object$y.hat.quants, y.hat.samples=object$y.hat.samples,
                    y.rep.quants=object$y.rep.quants, y.rep.samples=object$y.rep.samples, sub.sample=object$sub.sample, s.indx=object$s.indx))
    }
}  

residuals.NNGP <- function(object, sub.sample, ...){

    out <- fitted(object, sub.sample)   
    y <- object$y
        
    if(class(object)[3] == "binomial"){
        y <- y/object$weights
        y.hat.samples <- 1/(1+exp(-out$y.hat.samples))
    }
   
    residuals.samples <- y-out$y.hat.samples
    residuals.quants <- t(apply(residuals.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))
    
    return(list(residuals.samples=residuals.samples, residuals.quantiles=residuals.quants, sub.sample=out$sub.sample))
}

residuals.PGLogit <- function(object, sub.sample, ...){

    out <- fitted(object, sub.sample)
    y <- object$y
    
    y <- y/object$weights
    y.hat.samples <- 1/(1+exp(-out$y.hat.samples))
    residuals.samples <- y-out$y.hat.samples
    residuals.quants <- t(apply(residuals.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))
    
    return(list(residuals.samples=residuals.samples, residuals.quantiles=residuals.quants, sub.sample=out$sub.sample))    
}

print.PGLogit <- function(x, ...){
    cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), "", sep = "\n")
}

summary.PGLogit <- function(object, sub.sample, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), digits = max(3L, getOption("digits") - 3L), ...){
    
    print(object)
    
    n.samples <- nrow(object$p.beta.samples)

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

    cat("Chain sub.sample:\n")
    cat(paste("start = ",start,"\n", sep=""))
    cat(paste("end = ",end,"\n", sep=""))
    cat(paste("thin = ",thin,"\n", sep=""))
    cat(paste("samples size = ",length(s.indx),"\n", sep=""))
    
    ##print(summary(mcmc(object$p.beta.samples[s.indx,]), quantiles=quantiles), digits=digits)
    print(noquote(apply(t(apply(object$p.beta.samples[s.indx,], 2, function(x) quantile(x, prob=quantiles))), 2, function(x) formatC(x, format = "f", digits = digits))))
    
}   

print.NNGP <- function(x, ...){

    cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), "", sep = "\n")
    
    cat("Model class is",class(x)[1],class(x)[2],class(x)[3],"family.\n")
    
    if(class(x)[2] != "conjugate"){
        cat("\nModel object contains",nrow(x$p.theta.samples), "MCMC samples.\n")
    }
}

print.spPredict <- function(x, ...){

    cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), "", sep = "\n")
    
    cat("Predictions for a model of class",x$sp.obj.class[1],x$sp.obj.class[2],x$sp.obj.class[3],"family.\n")
    
    if(x$sp.obj.class[2] != "conjugate"){
        cat("\nspPredict object contains",ncol(x$p.y.0), "samples from the posterior predictive distribution.\n")
    }
}

summary.NNGP <- function(object, sub.sample, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), digits = max(3L, getOption("digits") - 3L), ...){

    elip.args <- list(...)
    
    print(object)
    cat("\n")
    
    if(class(object)[2] == "conjugate"){
        
        slgp <- FALSE
        w.str <- FALSE
        if("SLGP" %in% class(object)){
            slgp <- TRUE
            if("w.str" %in% names(elip.args)){
                w.str <- elip.args[["w.str"]]
            }
        }

        ##point estimate
        p <- ncol(object$X)
        
        if(slgp && w.str){
            tab <- data.frame("Estimate" = c(as.vector(object$beta.hat),object$sigma.sq.hat,object$theta.alpha),
                              "Variance" = c(diag(object$beta.var),object$sigma.sq.var,rep(0,length(object$theta.alpha))))
            rownames(tab) <- c(colnames(object$beta.hat),"sigma.sq",colnames(object$theta.alpha))
        }else{
            tab <- data.frame("Estimate" = c(as.vector(object$beta.hat)[1:p],object$sigma.sq.hat,object$theta.alpha),
                              "Variance" = c(diag(object$beta.var)[1:p],object$sigma.sq.var,rep(0,length(object$theta.alpha))))
            rownames(tab) <- c(colnames(object$beta.hat)[1:p],"sigma.sq",colnames(object$theta.alpha))
        }
        
        print(noquote(apply(tab, 2, function(x) formatC(x, format = "f", digits = digits))))
        cat("\n")
        
        ##mcmc summary
        update.samples <- FALSE
        
        if(!"p.beta.theta.samples" %in% names(object)){
            if(missing(sub.sample)){
                message("If posterior summaries are desired, then either rerun spConjNNGP with fit.rep=TRUE, or specify the summary argument sub.sample to indicate the number of fitted and replicated samples to collect.")
            }else{
                update.samples <- TRUE
            }
        }else{
            if(missing(sub.sample) || (!missing(sub.sample) && sub.sample == nrow(object$p.beta.theta.samples))){
                cat(paste("samples size = ",nrow(object$p.beta.theta.samples),"\n", sep=""))
                if(slgp && w.str){
                    print(noquote(apply(t(apply(cbind(object$p.beta.theta.samples,object$p.beta.theta.samples[,c("sigma.sq","tau.sq")]), 2, function(x) quantile(x, prob=quantiles))), 2, function(x) formatC(x, format = "f", digits = digits))))
                }else{##nngp and slgp without printing w.str
                    print(noquote(apply(t(apply(cbind(object$p.beta.theta.samples[,1:p,drop=FALSE],object$p.beta.theta.samples[,c("sigma.sq","tau.sq")]), 2, function(x) quantile(x, prob=quantiles))), 2, function(x) formatC(x, format = "f", digits = digits))))
                }
            }else{
                update.samples <- TRUE
            }
        }
        
        if(update.samples){
            out <- list()
            
            ##get samples
            b <- object$bb
            B.inv <- matrix(object$bB.inv, length(b), length(b))
            B.inv[upper.tri(B.inv, FALSE)] <- t(B.inv)[upper.tri(B.inv, FALSE)]
            alpha <- object$theta.alpha[,"alpha"]
            out$p.beta.theta.samples <- matrix(0, nrow=sub.sample, ncol=length(b)+2)
            
            for(i in 1:sub.sample){
                sigma.sq <- rigamma(1, object$ab[1,], object$ab[2,])
                beta <- rmvn(1, B.inv%*%b, sigma.sq*B.inv)
                tau.sq <- alpha*sigma.sq
                out$p.beta.theta.samples[i,] <- c(beta, sigma.sq, tau.sq)
            }
            
            colnames(out$p.beta.theta.samples) <- c(colnames(object$X), "sigma.sq", "tau.sq")
            
            cat(paste("samples size = ",sub.sample,"\n", sep=""))
            if(slgp && w.str){
                print(noquote(apply(t(apply(cbind(out$p.beta.theta.samples,out$p.beta.theta.samples[,c("sigma.sq","tau.sq")]), 2, function(x) quantile(x, prob=quantiles))), 2, function(x) formatC(x, format = "f", digits = digits))))
            }else{##nngp and slgp and without printing w.str
                print(noquote(apply(t(apply(cbind(out$p.beta.theta.samples[,1:p,drop=FALSE],out$p.beta.theta.samples[,c("sigma.sq","tau.sq")]), 2, function(x) quantile(x, prob=quantiles))), 2, function(x) formatC(x, format = "f", digits = digits))))
            }
            
            return(out)
        }
        
    }else{##not conjugate
        n.samples <- nrow(object$p.theta.samples)
        
        if(missing(sub.sample)){
            sub.sample <- list()
        }
        
        start <- ifelse(!"start" %in% names(sub.sample), floor(0.5*n.samples), sub.sample$start)
        end <- ifelse(!"end" %in% names(sub.sample), n.samples, sub.sample$end)
        thin <- ifelse(!"thin" %in% names(sub.sample), 1, sub.sample$thin)
        if(!is.numeric(start) || start >= n.samples){stop("invalid start")}
        if(!is.numeric(end) || end > n.samples){stop("invalid end")}
        if(!is.numeric(thin) || thin >= n.samples){stop("invalid thin")}
        
        s.indx <- seq(as.integer(start), as.integer(end), by=as.integer(thin))
        
        cat("Chain sub.sample:\n")
        cat(paste("start = ",start,"\n", sep=""))
        cat(paste("end = ",end,"\n", sep=""))
        cat(paste("thin = ",thin,"\n", sep=""))
        cat(paste("samples size = ",length(s.indx),"\n", sep=""))
        
        ##print(summary(mcmc(cbind(object$p.beta.samples, object$p.theta.samples)[s.indx,]), quantiles=quantiles), digits=digits)
        print(noquote(apply(t(apply(cbind(object$p.beta.samples, object$p.theta.samples)[s.indx,], 2, function(x) quantile(x, prob=quantiles))), 2, function(x) formatC(x, format = "f", digits = digits))))
    }
}

print.spDiag <- function(x, ...){
    
    cat("Chain sub.sample:\n")
    cat(paste("start = ",x$sub.sample[[1]],"\n", sep=""))
    cat(paste("end = ",x$sub.sample[[2]],"\n", sep=""))
    cat(paste("thin = ",x$sub.sample[[3]],"\n", sep=""))
    cat(paste("samples size = ",length(x$s.indx),"\n", sep=""))
    
    if("DIC" %in% names(x)){
        print(x$DIC)
        cat("\n")
    }
    if("WAIC" %in% names(x)){
        print(x$WAIC)
        cat("\n")
    }
    if("GP" %in% names(x)){
        print(x$GP)
        cat("\n")
    }
    if("GRS" %in% names(x)){
        print(x$GRS)
    }

}
