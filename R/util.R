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
