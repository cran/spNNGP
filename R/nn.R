mkNNIndx <- function(coords, m, n.omp.threads=1){

    n <- nrow(coords)
    nIndx <- (1+m)/2*m+(n-m-1)*m
    nnIndx <- rep(0, nIndx)
    nnDist <- rep(0, nIndx)
    nnIndxLU <- matrix(0, n, 2)

    n <- as.integer(n)
    m <- as.integer(m)
    coords <- as.double(coords)
    nnIndx <- as.integer(nnIndx)
    nnDist <- as.double(nnDist)
    nnIndxLU <- as.integer(nnIndxLU)
    n.omp.threads <- as.integer(n.omp.threads)
    
    ptm <- proc.time()
    
    out <- .Call("mkNNIndx", n, m, coords, nnIndx, nnDist, nnIndxLU, n.omp.threads)

    run.time <- proc.time() - ptm
    
    list("run.time"=run.time, "nnIndx"=as.integer(nnIndx), "nnDist"=as.double(nnDist), "nnIndxLU"=nnIndxLU)

}

mkNNIndxCB <- function(coords, m, n.omp.threads=1){
    
    n <- nrow(coords)
    nIndx <- (1+m)/2*m+(n-m-1)*m
    nnIndx <- rep(0, nIndx)
    nnDist <- rep(0, nIndx)
    nnIndxLU <- matrix(0, n, 2)
    
    n <- as.integer(n)
    m <- as.integer(m)
    coords <- as.double(coords)
    nnIndx <- as.integer(nnIndx)
    nnDist <- as.double(nnDist)
    nnIndxLU <- as.integer(nnIndxLU)
    n.omp.threads <- as.integer(n.omp.threads)

    ptm <- proc.time()
    
    out <- .Call("mkNNIndxCB", n, m, coords, nnIndx, nnDist, nnIndxLU, n.omp.threads)

    run.time <- proc.time() - ptm
    
    list("run.time"=run.time, "nnIndx"=as.integer(nnIndx), "nnDist"=as.double(nnDist), "nnIndxLU"=nnIndxLU)
}
