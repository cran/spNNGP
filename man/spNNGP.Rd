\name{spNNGP}
\alias{spNNGP}
\title{Function for fitting univariate Bayesian spatial regression
  models}

\description{
  The function \code{spNNGP} fits Gaussian univariate Bayesian spatial
  regression models using Nearest Neighbor Gaussian Processes (NNGP).
}

\usage{
spNNGP(formula, data = parent.frame(), coords, method = "response", n.neighbors = 15,
      starting, tuning, priors, cov.model = "exponential",
      n.samples, n.omp.threads = 1, search.type = "tree",
      return.neighbors = FALSE, verbose = TRUE, n.report = 100, ...)
}

\arguments{
  \item{formula}{a symbolic description of the regression model to be
    fit. See example below. }
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which \code{spNNGP} is called.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix of the observation coordinates
    in \eqn{R^2}{R^2} (e.g., easting and northing). }
  \item{method}{a quoted keyword that specifies the NNGP sampling algorithm. Supported method keywords are: \code{"response"} and
  \code{"sequential"}. See below for details.}
  \item{n.neighbors}{number of neighbors used in the NNGP.}
  \item{starting}{a list with each tag corresponding to a
    parameter name. Valid tags are \code{beta}, \code{sigma.sq},
    \code{tau.sq}, \code{phi}, and \code{nu}. \code{nu} is only specified if \code{cov.model="matern"}. The value portion of each
  tag is the parameter's starting value. }
  \item{tuning}{a list with each tag corresponding to a
    parameter name. Valid tags are \code{sigma.sq},
    \code{tau.sq}, \code{phi}, and \code{nu}. If \code{method="sequential"}
  then only \code{phi} and \code{nu} need to be specified. The value portion of each
  tag defines the variance of the Metropolis sampler Normal proposal
  distribution. }
  \item{priors}{a list with each tag corresponding to a
    parameter name. Valid tags are \code{sigma.sq.ig},
    \code{tau.sq.ig}, \code{phi.unif}, and \code{nu.unif}. Variance parameters, \code{simga.sq} and
    \code{tau.sq}, are assumed to follow an
    inverse-Gamma distribution, whereas the spatial decay \code{phi}
    and smoothness \code{nu} parameters are assumed to follow Uniform distributions. The
    hyperparameters of the inverse-Gamma are
    passed as a vector of length two, with the first and second elements corresponding
    to the \emph{shape} and \emph{scale}, respectively. The hyperparameters
    of the Uniform are also passed as a vector of length two with the first
    and second elements corresponding to the lower and upper support,
    respectively.
 }
  \item{cov.model}{a quoted keyword that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.  See below for details.}
  \item{n.samples}{the number of posterior samples to collect.}
  \item{n.omp.threads}{a positive integer indicating
the number of threads to use for SMP parallel processing. The package must
be compiled for OpenMP support. For most Intel-based machines, we
  recommend setting \code{n.omp.threads} up to the number of hyperthreaded cores.}
  \item{verbose}{if \code{TRUE}, model specification and progress of the
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.}
  \item{search.type}{a quoted keyword that specifies type of nearest
    neighbor search algorithm. Supported method key words are: \code{"tree"} and
    \code{"brute"} both will yield the same solution but \code{"tree"}
    should be much faster.}
  \item{return.neighbors}{if \code{TRUE}, a list containing the indices
    for each locations' nearest neighbors will be returned. See
    \code{n.indx} below for more details. }
  \item{n.report}{the interval to report Metropolis sampler acceptance and MCMC progress.}  
  \item{...}{currently no additional arguments.  }
}

\value{
  An object of class \code{rNNGP} or \code{sNNGP} depending on the \code{method}, which is a list comprising:
  \item{p.beta.samples}{a \code{coda} object of posterior samples for
    the regression coefficients. }
  \item{p.theta.samples}{a \code{coda} object of posterior samples for
    covariance parameters. }
  \item{n.indx}{if \code{return.neighbors=TRUE} then \code{n.indx} will
    be a list of length \eqn{n}{n}. The \emph{i}-th element in the
    list corresponds to the \emph{i}-th row in the ordered
    \code{coords} matrix and the elements are the nearest neighbor indices
    for the given location. }
  \item{coords}{the input \code{coords} matrix ordered by the first column. }
  \item{run.time}{execution time for building the nearest neighbor index
    and MCMC sampler reported using \code{proc.time()}.}
  The return object will include additional objects used for subsequent
  prediction and/or model fit evaluation.
}

\details{
Model parameters can be fixed at their \code{starting} values by setting their
\code{tuning} values to zero.

The \emph{no nugget} model is specified by setting \code{tau.sq} to zero
in the \code{starting} and \code{tuning} lists.
}

\references{
  Datta, A., S. Banerjee, A.O. Finley, and A.E. Gelfand. (2016)
  Hierarchical Nearest-Neighbor Gaussian process models for large
  geostatistical datasets. Journal of the American Statistical
  Association, 111:800-812

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

##Fit a Response and Sequential NNGP model
n.samples <- 500

starting <- list("phi"=phi, "sigma.sq"=5, "tau.sq"=1)

tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)

priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1))

cov.model <- "exponential"

m.s <- spNNGP(y~x-1, coords=coords, starting=starting, method="sequential", n.neighbors=10,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, n.omp.threads=2)

round(summary(m.s$p.beta.samples)$quantiles[,c(3,1,5)],2)
round(summary(m.s$p.theta.samples)$quantiles[,c(3,1,5)],2)

m.r <- spNNGP(y~x-1, coords=coords, starting=starting, method="response", n.neighbors=10,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, n.omp.threads=2)

round(summary(m.r$p.beta.samples)$quantiles[,c(3,1,5)],2)
round(summary(m.r$p.theta.samples)$quantiles[,c(3,1,5)],2)

}

\keyword{model}