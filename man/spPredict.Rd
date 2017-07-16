\name{spPredict}
\alias{spPredict}
\title{Function for prediction at new locations using \code{spNNGP} models.}

\description{
  The function \code{spPredict} collects posterior predictive samples
  for a set of new locations given a \code{\link{spNNGP}} object.
}

\usage{
spPredict(sp.obj, X.0, coords.0, start=1, end, thin=1, n.omp.threads = 1,
          verbose=TRUE, n.report=100, ...)
}

\arguments{
  \item{sp.obj}{an object returned by \code{\link{spNNGP}}. }
  \item{X.0}{the design matrix for prediction locations. An
    intercept should be provided in the first column if one is specified
    in \code{sp.obj} model.}
  \item{coords.0}{the spatial coordinates corresponding to
    \code{X.0}. }
  \item{start}{specifies the first sample included in the composition sampling.}
  \item{end}{specifies the last sample included in the composition.
    The default is to use all posterior samples in \code{sp.obj}. }
  \item{thin}{a sample thinning factor.  The default of 1 considers all
    samples between \code{start} and \code{end}.  For example, if \code{thin = 10} then 1 in 10 samples are considered between \code{start} and
    \code{end}. }
  \item{n.omp.threads}{a positive integer indicating
    the number of threads to use for SMP parallel processing. The package must
    be compiled for OpenMP support. For most Intel-based machines, we recommend setting
    \code{n.omp.threads} up to the number of hyperthreaded cores.}
  \item{verbose}{if \code{TRUE}, model specification and progress of the
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.  }
  \item{n.report}{the interval to report sampling progress. }  
  \item{...}{currently no additional arguments.  }
}

\value{
  A list comprising:
  \item{p.y.0}{a matrix that holds the response variable posterior
    predictive samples where rows are locations corresponding to
  \code{coords.0} and columns are samples. }
  \item{p.w.0}{a matrix that holds the random effect posterior
    predictive samples where rows are locations corresponding to
  \code{coords.0} and columns are samples. This is only returned if \code{\link{spNNGP}} \code{method = "sequential".}}
\item{run.time}{execution time reported using \code{proc.time()}.}
}

\references{
  Datta, A., S. Banerjee, A.O. Finley, and A.E. Gelfand. (2016)
  Hierarchical Nearest-Neighbor Gaussian process models for large
  geostatistical datasets. Journal of the American Statistical
  Association, 111:800-812.
  
  Finley, A.O., A. Datta, B.C. Cook, D.C. Morton, H.E. Andersen, and
  S. Banerjee (2017) Applying Nearest Neighbor Gaussian Processes to
  massive spatial data sets: Forest canopy height prediction across
  Tanana Valley Alaska, \url{https://arxiv.org/abs/1702.00434v2}.
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Abhirup Datta \email{abhidatta@jhu.edu}, \cr
  Sudipto Banerjee \email{sudipto@ucla.edu}
}

\examples{

rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)\%*\%D + rep(mu,rep(n,p)))
}

##Make some data
set.seed(1)
n <- 100
coords <- cbind(runif(n,0,1), runif(n,0,1))

x <- cbind(1, rnorm(n))

B <- as.matrix(c(1,5))

sigma.sq <- 5
tau.sq <- 1
phi <- 3/0.5

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, x\%*\%B + w, sqrt(tau.sq))

ho <- sample(1:n, 50)

y.ho <- y[ho]
x.ho <- x[ho,,drop=FALSE]
w.ho <- w[ho]
coords.ho <- coords[ho,]

y <- y[-ho]
x <- x[-ho,,drop=FALSE]
w <- w[-ho,,drop=FALSE]
coords <- coords[-ho,]

##Fit a Response and Sequential NNGP model
n.samples <- 500

starting <- list("phi"=phi, "sigma.sq"=5, "tau.sq"=1)

tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)

priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1))

cov.model <- "exponential"

n.report <- 500

##Predict for holdout set using both models
m.s <- spNNGP(y~x-1, coords=coords, starting=starting, method="sequential", n.neighbors=10,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, n.omp.threads=2, n.report=n.report)

m.r <- spNNGP(y~x-1, coords=coords, starting=starting, method="response", n.neighbors=10,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, n.omp.threads=2, n.report=n.report)

##Prediction for holdout data
p.s <- spPredict(m.s, X.0 = x.ho, coords.0 = coords.ho, n.omp.threads=2)

plot(apply(p.s$p.w.0, 1, mean), w.ho)
plot(apply(p.s$p.y.0, 1, mean), y.ho)

p.r <- spPredict(m.r, X.0 = x.ho, coords.0 = coords.ho, n.omp.threads=2)

plot(apply(p.r$p.y.0, 1, mean), y.ho)

}

\keyword{model}