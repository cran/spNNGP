rigamma <- function(n, a, b){
    1/rgamma(n = n, shape = a, rate = b)
}

rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

GP <- function(Y.rep, y){
    mu.rep <- apply(Y.rep, 1, mean)
    var.rep <- apply(Y.rep, 1, var)
    G <- sum((y-mu.rep)^2)
    P <- sum(var.rep)
    D <- G+P
    
    GPD <- matrix(0,3,1)
    rownames(GPD) <- c("G","P","D")
    colnames(GPD) <- "value"
    GPD[1,1] <- G
    GPD[2,1] <- P
    GPD[3,1] <- D
    as.data.frame(GPD)
}

GR <- function(Y.rep, y){
    mu.rep <- apply(Y.rep, 1, mean)
    var.rep <- apply(Y.rep, 1, var)
    GRS <- as.data.frame(-sum(((y-mu.rep)/sqrt(var.rep))^2) - sum(log(var.rep)))

    rownames(GRS) <- "GRS"
    colnames(GRS) <- "value"
    GRS
}

PGLogitDiag <- function(object, sub.sample, ...){

    out <- list()
    
    ##for DIC and WAIC
    if(missing(sub.sample)){
        if("s.indx" %in% names(object)){
            out$y.hat.samples <- object$y.hat.samples
            out$s.indx <- object$s.indx
            out$sub.sample <- object$sub.sample
        }else{
            out <- c(out, fitted(object))#uses default sub.sample start, end, thin in fitted
        }
    }else{
        out <- c(out, fitted(object, sub.sample=sub.sample))
    }

    n.samples <- length(out$s.indx)

    y <- object$y
    X <- object$X
    n <- nrow(X)

    ##############
    ##DIC
    ##############
    beta.mu <- apply(object$p.beta.samples[out$s.indx,,drop=FALSE], 2, mean)
    
    prob <- 1/(1+exp(-(X%*%beta.mu)))
    L <- sum(dbinom(y, size=object$weights, prob=prob, log=TRUE));
    
    llSum <- 0
    for(i in 1:n.samples){
        prob <- 1/(1+exp(-(out$y.hat.samples[,i])))
        llSum <- llSum+sum(dbinom(y, size=object$weights, prob=prob, log=TRUE))
    }
    
    P <- 2*(L-(1/n.samples*llSum))
    DIC <- data.frame(value=c(-2*(L-P), P, L))
    rownames(DIC) <- c("DIC", "pD", "D")
    out$DIC <- DIC
    
    ##############
    ##WAIC
    ##############
    LPPD <- 0
    P1 <- 0
    P2 <- 0
    
    for(i in 1:n){
        
        prob <- 1/(1+exp(-(out$y.hat.samples[i,])))
        L <- dbinom(y[i], size=object$weights[i], prob=prob, log=FALSE)
        
        LPPD <- LPPD+log(mean(L))
        
        a <- log(mean(L))
        b <- mean(log(L))
        
        P1 <- P1 + 2*(a-b)
        P2 <- P2 + var(log(L))
    }
    
    WAIC.1 <- -2*(LPPD-P1)
    WAIC.2 <- -2*(LPPD-P2)
    
    out$WAIC <- data.frame(value=c(WAIC.1, WAIC.2, P1, P2, LPPD))
    row.names(out$WAIC) <- c("WAIC.1","WAIC.2","P.1","P.2","LPPD")
    
    ## ##############
    ## ##GPD and GRS
    ## ##############
    ## out$GP <- GP(out$y.rep.samples, y)
    ## out$GRS <- GR(out$y.rep.samples, y)

    class(out) <- "spDiag"
    out
   
}

spDiag <- function(object, sub.sample, ...){

    if(class(object)[1] == "PGLogit"){
        return(PGLogitDiag(object, sub.sample))
    }
    
    out <- list()

    if(missing(sub.sample)){
        if("s.indx" %in% names(object)){
            out$y.hat.samples <- object$y.hat.samples
            out$y.rep.samples <- object$y.rep.samples
            out$s.indx <- object$s.indx
            out$sub.sample <- object$sub.sample
        }else{
            out <- c(out, fitted(object))#uses default sub.sample start, end, thin in fitted
        }
    }else{
        out <- c(out, fitted(object, sub.sample=sub.sample))
    }

    n.samples <- length(out$s.indx)
    
    ##if(class(object)[2] == "latent"){
    if(object$type[1] == "latent"){                 
        n.samples <- length(out$s.indx)
        y <- object$y
        X <- object$X
        n <- nrow(X)
        
        ##############
        ##DIC
        ##############
        beta <- object$p.beta.samples[out$s.indx,,drop=FALSE]
        w <- object$p.w.samples[,out$s.indx,drop=FALSE]
        
        beta.mu <- apply(beta, 2, mean)
        w.mu <- apply(w, 1, mean)
            
        ##if(class(object)[3] == "gaussian"){
        if(object$type[2] == "gaussian"){
            tau.sq <- object$p.theta.samples[out$s.indx,"tau.sq"]
            tau.sq.mu <- mean(tau.sq) 
            
            L <- sum(dnorm(y, X%*%beta.mu+w.mu, sqrt(tau.sq.mu), log=TRUE))
            
            llSum <- 0
            for(i in 1:n.samples){
                llSum <- llSum+sum(dnorm(y, X%*%t(beta[i,,drop=FALSE])+w[,i], sqrt(tau.sq[i]), log=TRUE))
            }
            
        }else{##binomial
            prob <- 1/(1+exp(-(X%*%beta.mu+w.mu)))
            L <- sum(dbinom(y, size=object$weights, prob=prob, log=TRUE));
            
            llSum <- 0
            for(i in 1:n.samples){
                prob <- 1/(1+exp(-(X%*%t(beta[i,,drop=FALSE])+w[,i])))
                llSum <- llSum+sum(dbinom(y, size=object$weights, prob=prob, log=TRUE))
            }
        }           
        
        P <- 2*(L-(1/n.samples*llSum))
        DIC <- data.frame(value=c(-2*(L-P), P, L))
        rownames(DIC) <- c("DIC", "pD", "L")
        out$DIC <- DIC

        ##############
        ##WAIC
        ##############
        LPPD <- 0
        P1 <- 0
        P2 <- 0

        for(i in 1:n){
            w <- object$p.w.samples[i,out$s.indx,drop=FALSE]
                  
            ##if(class(object)[3] == "gaussian"){
            if(object$type[2] == "gaussian"){
                tau.sq <- object$p.theta.samples[out$s.indx,"tau.sq"]
                L <- dnorm(y[i], as.vector(X[i,]%*%t(beta)+w), sqrt(tau.sq), log=FALSE)
            }else{##binomial
                prob <- as.vector(1/(1+exp(-(X[i,]%*%t(beta)+w))))
                L <- dbinom(y[i], size=object$weights[i], prob=prob, log=FALSE)
            }           
            
            LPPD <- LPPD+log(mean(L))##log pointwise predictive density
            
            a <- log(mean(L))
            b <- mean(log(L))
            
            P1 <- P1 + 2*(a-b)
            P2 <- P2 + var(log(L))
        }
        
        WAIC.1 <- -2*(LPPD-P1)
        WAIC.2 <- -2*(LPPD-P2)
        
        out$WAIC <- data.frame(value=c(WAIC.1, WAIC.2, P1, P2, LPPD))
        row.names(out$WAIC) <- c("WAIC.1","WAIC.2","P.1","P.2","LPPD")

        ##############
        ##GPD and GRS
        ##############
              
        ##if(class(object)[3] == "gaussian"){
        if(object$type[2] == "gaussian"){
            out$GP <- GP(out$y.rep.samples, y)
            out$GRS <- GR(out$y.rep.samples, y)
        }
               
    }else{##response and conj
        
        ##############
        ##GPD and GRS
        ##############
 
        ##if(class(object)[3] == "gaussian"){
        if(object$type[2] == "gaussian"){
            out$GP <- GP(out$y.rep.samples, object$y)
            out$GRS <- GR(out$y.rep.samples, object$y)
        }
    }

    class(out) <- "spDiag"
    out
}

crps <- function(y, y.hat, y.var){
    sd <- sqrt(y.var)
    y.std <- (y-y.hat)/sd
    mean(-sd*(1/sqrt(pi) - 2*dnorm(y.std) - y.std*(2*pnorm(y.std) - 1)))
}

rmspe <- function(y, y.hat){
    sqrt(mean((y-y.hat)^2))
}

get.n.indx <- function(i, m){
    i <- i-1
    if(i == 0){
        return(NA)
    }else if(i < m){
        n.indx.i <- i/2*(i-1)
        m.i <- i
        return((n.indx.i+1):((n.indx.i+1)+i-1))
    }else{
        n.indx.i <- m/2*(m-1)+(i-m)*m
        m.i <- m
        return((n.indx.i+1):((n.indx.i+1)+m-1))
    }
}

mk.n.indx.list <- function(n.indx, n, m){
    n.indx.list <- vector("list", n)
    n.indx.list[1] <- NA
    for(i in 2:n){
        n.indx.list[[i]] <- n.indx[get.n.indx(i, m)]+1
    }
    n.indx.list
}
