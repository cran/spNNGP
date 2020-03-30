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
    
    ##if(class(object)[2] == "conjugate"){
    if(object$type[1] == "conjugate"){
        
        if(!missing(sub.sample)){
            if(!(is.atomic(sub.sample) && length(sub.sample) == 1L && sub.sample > 0)){
                stop("error: for a conjugate model, the sub.sample argument should be in integer value indicating the number of fitted and replicated samples to collect.")
            }
        }
                       
        ##if("SLGP" %in% class(object)){
        if("SLGP" %in% object$type){

            if("y.hat.samples" %in% names(object)){
                if(missing(sub.sample) || (!missing(sub.sample) && sub.sample == ncol(object$y.hat.samples))){
                    return(list(y.hat.quants=object$y.hat.quants, y.hat.samples=object$y.hat.samples))
                }
            }

            ##y.rep.samples is missing or requesting a different number of samples
            if(missing(sub.sample)){
                stop("error: for a conjugate model created with fit.rep=FALSE, the sub.sample argument should be provided to indicate the number of fitted and replicated samples to collect.")
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
                stop("error: for an conjugate model created with fit.rep=FALSE, the sub.sample argument should be provided to indicate the number of fitted and replicated samples to collect.")
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
        
        ##if(class(object)[2] == "latent"){
        if(object$type[1] == "latent"){
            
            if(update.fit.rep){
                beta <- as.matrix(object$p.beta.samples)[s.indx,,drop=FALSE]
                w <- object$p.w.samples[,s.indx,drop=FALSE]
                
                X <- object$X
                y.hat.samples <- X%*%t(beta) + w
                y.hat.quants <- t(apply(y.hat.samples, 1, function(x) quantile(x, prob=c(0.5, 0.05, 0.975))))
                
                ##if(class(object)[2] == "gaussian"){##this was a bug
                if(object$type[2] == "gaussian"){
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
        
    ##if(class(object)[3] == "binomial"){
    if(object$type[2] == "binomial"){
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
    
    ##cat("Model class is",class(x)[1],class(x)[2],class(x)[3],"family.\n")
    cat("Model class is ",class(x)[1],", method ", x$type[1],", family ",x$type[2],".\n", sep="")
    
    ##if(class(x)[2] != "conjugate"){
    if(x$type[1] != "conjugate"){
        cat("\nModel object contains",nrow(x$p.theta.samples), "MCMC samples.\n")
    }
}

summary.NNGP <- function(object, sub.sample, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), digits = max(3L, getOption("digits") - 3L), ...){

    elip.args <- list(...)
    
    print(object)
    cat("\n")
    
    ##if(class(object)[2] == "conjugate"){
    if(object$type[1] == "conjugate"){
        slgp <- FALSE
        w.str <- FALSE
        ##if("SLGP" %in% class(object)){
        if("SLGP" %in% object$type){
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

predict.NNGP <- function(object, X.0, coords.0, sub.sample, n.omp.threads = 1, verbose=TRUE, n.report=100, ...){
  
    ####################################################
    ##Check for unused args
    ####################################################
    formal.args <- names(formals(sys.function(sys.parent())))

    elip.args <- list(...)
    for(i in names(elip.args)){
        if(! i %in% formal.args)
            warning("'",i, "' is not an argument")
    }
    
    if(missing(object)){stop("error: predict expects object\n")}
    if(!class(object)[1] == "NNGP"){
        stop("error: requires an output object of class NNGP\n")
    }


    ##call
    out <- list()
    out$call <- match.call()
    out$object.class <- class(object)
    out$type <- object$type
    
    ##conjugate
    ##if(class(object)[2] == "conjugate"){
    if(object$type[1] == "conjugate"){

        if(!missing(sub.sample)){
            warning("'sub.sample' is not an argument for prediction using a conjugate model.")
        }
        
        theta.alpha <- as.vector(object$theta.alpha)
        names(theta.alpha) <- colnames(object$theta.alpha)
        
        ptm <- proc.time()
        ##if(length(class(object)) == 4 & class(object)[4] == "SLGP"){
        if("SLGP" %in% out$type){
            out <- c(out, spConjNNGP(object$y ~ object$X-1, coords=object$coords, knots=object$knots, sigma.sq.IG=object$sigma.sq.IG, n.neighbors=object$n.neighbors,
                                     X.0 = X.0, coords.0=coords.0,
                                     theta.alpha=theta.alpha, cov.model=object$cov.model, n.omp.threads=n.omp.threads, search.type=object$search.type, verbose=verbose))
        }else{
            out <- c(out, spConjNNGP(object$y ~ object$X-1, coords=object$coords, sigma.sq.IG=object$sigma.sq.IG, n.neighbors=object$n.neighbors,
                                     X.0 = X.0, coords.0=coords.0,
                                     theta.alpha=theta.alpha, cov.model=object$cov.model, n.omp.threads=n.omp.threads, search.type=object$search.type, verbose=verbose))
        }
        
        out$run.time <- proc.time() - ptm
    }else{ ##sequential and response models
        
        X <- object$X
        y <- object$y
        coords <- object$coords
        family <- object$family ##family code gaussian=1, binomial=2, ...
        
        family.indx <- 1
        ##if(class(object)[3] == "binomial"){
        if(out$type[2] == "binomial"){
            family.indx <- 2
        }
        
        n <- nrow(X)
        p <- ncol(X)
        
        p.theta.samples <- object$p.theta.samples
        p.beta.samples <- object$p.beta.samples
        n.samples <- nrow(p.beta.samples)
        ##if(class(object)[2] == "latent"){
        if(out$type[1] == "latent"){
            p.w.samples <- object$p.w.samples
        }    
        n.neighbors <- object$n.neighbors
        cov.model.indx <- object$cov.model.indx
        
        ##subsamples
        if(missing(sub.sample)){
            sub.sample <- list()
        }
     
        start <- ifelse(!"start" %in% names(sub.sample), 1, sub.sample$start)
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
        
        ##if(class(object)[2] == "latent"){
        if(out$type[1] == "latent"){
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
        ##if(class(object)[2] == "latent"){
        if(out$type[1] == "latent"){
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
        
        ##if(class(object)[2] == "latent"){
        if(out$type[1] == "latent"){
            out <- c(out, .Call("sNNGPPredict", X, y, coords, n, p, n.neighbors, X.0, coords.0, q, nn.indx.0, 
                                p.beta.samples, p.theta.samples, p.w.samples, n.samples, family.indx, cov.model.indx, n.omp.threads, verbose, n.report))
        }else{
            out <- c(out, .Call("rNNGPPredict", X, y, coords, n, p, n.neighbors, X.0, coords.0, q, nn.indx.0, 
                                p.beta.samples, p.theta.samples, n.samples, cov.model.indx, n.omp.threads, verbose, n.report))
        }
        
        out$run.time <- proc.time() - ptm
    }

    class(out) <- "predict.NNGP"
    out
}


print.predict.NNGP <- function(x, ...){

    cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), "", sep = "\n")
    
    ##cat("Predictions for a model of class",x$object.class[1],x$object.class[2],x$object.class[3],"family.\n")
    cat("Predictions for a model of class ",x$object.class[1],", method ", x$type[1],", family ",x$type[2],".\n", sep="")
    
    ##if(x$object.class[1] != "conjugate"){
    if(x$type[1] != "conjugate"){
        cat("\npredict.NNGP object contains",ncol(x$p.y.0), "samples from the posterior predictive distribution.\n")
    }
}
