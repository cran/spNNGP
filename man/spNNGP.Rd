\name{spNNGP}
\alias{spNNGP}
\title{Function for Fitting Univariate Bayesian Spatial Regression
  Models}

\description{
  The function \code{spNNGP} fits Gaussian and non-Gaussian univariate Bayesian spatial
  regression models using Nearest Neighbor Gaussian Processes (NNGP).
}

\usage{
spNNGP(formula, data = parent.frame(), coords, method = "response",
      family="gaussian", weights, n.neighbors = 15, 
      starting, tuning, priors, cov.model = "exponential",
      n.samples, n.omp.threads = 1, search.type = "cb", ord,
      return.neighbor.info = FALSE, neighbor.info,
      fit.rep = FALSE, sub.sample, verbose = TRUE, n.report = 100, ...)
}

\arguments{
  \item{formula}{a symbolic description of the regression model to be
    fit. See example below. }
  
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which \code{spNNGP} is called.  }
  
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix of the observation
    coordinates in \eqn{R^2}{R^2} (e.g., easting and northing), or if
    \code{data} is a data frame then \code{coords} can be a vector of
    length two comprising coordinate column names or indices. There can
    be no duplicate locations. }
  
  \item{method}{a quoted keyword that specifies the NNGP sampling
    algorithm. Supported method keywords are: \code{"response"} and
    \code{"latent"}. When \eqn{n}{n} is large, the \code{"response"}
    algorithm should be faster and provide finer control over Metropolis
    acceptance rate for covariance parameters. In general,
    unless estimates of spatial random effects are needed, the
    \code{"response"} algorithm should be used. }
  
  \item{family}{a quoted keyword that specifies the data
    likelihood. Choices are "gaussian" for continuous outcome and
    "binomial" for discrete outcome which assumes a logistic link modeled
    using Polya-Gamma latent variables. }
  
  \item{weights}{specifies the number of trials for each observation when
    \code{family="binomial"}. The default is 1 trial for each
    observation. Valid arguments are a scalar value that specifies the number
    of trials used if all observations have the same number of trials, and a
    vector of length \eqn{n}{n} that specifies the number of trials for
    each observation used there are differences in the number of trials. }
  
  \item{n.neighbors}{number of neighbors used in the NNGP.}
  
  \item{starting}{a list with each tag corresponding to a parameter name. Valid tags are \code{beta}, \code{sigma.sq},
    \code{tau.sq}, \code{phi}, and \code{nu}. \code{nu} is only
    specified if \code{cov.model="matern"}. The value portion of each
    tag is the parameter's startingvalue. }
  
  \item{tuning}{a list with each tag corresponding to a parameter
    name. Valid tags are \code{sigma.sq}, \code{tau.sq}, \code{phi}, and
    \code{nu}. If \code{method="latent"} then only \code{phi} and
    \code{nu} need to be specified. The value portion of each tag defines
    the variance of the Metropolis sampler Normal proposal distribution. }
  
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
    respectively.}
  
  \item{cov.model}{a quoted keyword that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.  See below for details.}
  
  \item{n.samples}{the number of posterior samples to collect.}
  
  \item{n.omp.threads}{a positive integer indicating
    the number of threads to use for SMP parallel processing. The package must
    be compiled for OpenMP support. For most Intel-based machines, we
    recommend setting \code{n.omp.threads} up to the number of
    hyperthreaded cores. Note, \code{n.omp.threads} > 1 might not
    work on some systems.}
  
  \item{fit.rep}{if \code{TRUE}, regression fitted and replicate data will be
    returned. The argument \code{sub.sample} controls which MCMC samples
    are used to generate the fitted and replicated data. }
  
  \item{sub.sample}{an optional list that specifies the samples used
    for \code{fit.rep}. Valid tags are \code{start},
    \code{end}, and \code{thin}. Given the value associated with the tags,
    the sample subset is selected using \code{seq(as.integer(start),
      as.integer(end), by=as.integer(thin))}. The default values are
    \code{start=floor(0.5*n.samples)}, \code{end=n.samples} and \code{thin=1}.}
  
  \item{verbose}{if \code{TRUE}, model specification and progress of the
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.}
  
    \item{search.type}{a quoted keyword that specifies type of nearest
      neighbor search algorithm. Supported method key words are: \code{"cb"} and
      \code{"brute"}. The \code{"cb"} should generally be much
      faster. If locations do not have identical coordinate values on
      the axis used for the nearest neighbor ordering (see
      \code{ord} argument) then \code{"cb"} and \code{"brute"} should
      produce identical neighbor sets. However, if there are identical
      coordinate values on the axis used for nearest neighbor ordering,
      then \code{"cb"} and \code{"brute"} might produce different, but
      equally valid, neighbor sets, e.g., if data are on a grid. }
    
  \item{ord}{an index vector of length \eqn{n}{n} used for the nearest neighbor
    search. Internally, this vector is used to order \code{coords},
    i.e., \code{coords[ord,]}, and associated data. Nearest neighbor
    candidates for the \emph{i}-th row in the ordered \code{coords} are
    rows \code{1:(i-1)}, with the \code{n.neighbors} nearest neighbors being
    those with the minimum euclidean distance to the location defined by
    ordered \code{coords[i,]}. The default used when \code{ord} is not
    specified is x-axis ordering, i.e.,
    \code{order(coords[,1])}. This argument should typically be left
    blank. This argument will be ignored if the \code{neighbor.info}
    argument is used. }
  
  \item{return.neighbor.info}{if \code{TRUE}, a list called
    \code{neighbor.info} containing several data structures used for
    fitting the NNGP model is returned. If there is no change in input
    data or \code{n.neighbors}, this list can be passed to subsequent \code{spNNGP} calls via the \code{neighbor.info} argument to avoid
    the neighbor search, which can be time consuming if \eqn{n}{n} is
    large. In addition to the several cryptic data structures in
    \code{neighbor.info} there is a list called \code{n.indx} that is of
    length \eqn{n}{n}. The \emph{i}-th element in \code{n.indx}
    corresponds to the \emph{i}-th row in \code{coords[ord,]} and holds
    the vector of that location's nearest neighbor indices. This list can be useful for plotting the
    neighbor graph if desired. }
  
    \item{neighbor.info}{see the \code{return.neighbor.info} argument
      description above.}    

    \item{n.report}{the interval to report Metropolis sampler acceptance
      and MCMC progress.}  
  
  \item{...}{currently no additional arguments.  }
}

\value{
  An object of class \code{NNGP} with additional class designations for
  \code{method} and \code{family}. The return object is a list comprising:
  \item{p.beta.samples}{a \code{coda} object of posterior samples for
    the regression coefficients. }
  
  \item{p.theta.samples}{a \code{coda} object of posterior samples for
    covariance parameters. }
  
  \item{p.w.samples}{is a matrix of posterior samples for the spatial
    random effects, where rows correspond to locations in \code{coords} and columns hold the \code{n.samples} posterior samples. This is only returned if \code{method="latent"}.}

  \item{y.hat.samples}{if \code{fit.rep=TRUE}, regression fitted values from posterior samples
    specified using \code{sub.sample}. See additional details below. }

  \item{y.hat.quants}{if \code{fit.rep=TRUE}, 0.5, 0.025, and 0.975 quantiles of the \code{y.hat.samples}. }

  \item{y.rep.samples}{if \code{fit.rep=TRUE}, replicated outcome from posterior samples
    specified using \code{sub.sample}. See additional details below. }
  
  \item{y.rep.quants}{if \code{fit.rep=TRUE}, 0.5, 0.025, and 0.975
    quantiles of the \code{y.rep.samples}. }

  \item{s.indx}{if \code{fit.rep=TRUE}, the subset index specified with \code{sub.sample}.}
  
  \item{neighbor.info}{returned if \code{return.neighbor.info=TRUE} see the \code{return.neighbor.info} argument
    description above.}
    
  \item{run.time}{execution time for parameter estimation reported using
    \code{proc.time()}. This time does not include nearest neighbor
    search time for building the neighbor set.}
  
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
  Finley, A.O., A. Datta, S. Banerjee (2022) spNNGP R Package for
  Nearest Neighbor Gaussian Process Models. \emph{Journal of Statistical
  Software}, \doi{10.18637/jss.v103.i05}.
  
  Datta, A., S. Banerjee, A.O. Finley, and A.E. Gelfand. (2016)
  Hierarchical Nearest-Neighbor Gaussian process models for large
  geostatistical datasets. \emph{Journal of the American Statistical
    Association}, \doi{10.1080/01621459.2015.1044091}.
  
  Finley, A.O., A. Datta, B.D. Cook, D.C. Morton, H.E. Andersen, and
  S. Banerjee. (2019) Efficient algorithms for Bayesian Nearest Neighbor
  Gaussian Processes. \emph{Journal of Computational and Graphical
    Statistics}, \doi{10.1080/10618600.2018.1537924}.
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

##Fit a Response and Latent NNGP model
n.samples <- 500

starting <- list("phi"=phi, "sigma.sq"=5, "tau.sq"=1)

tuning <- list("phi"=0.5, "sigma.sq"=0.5, "tau.sq"=0.5)

priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1))

cov.model <- "exponential"

m.s <- spNNGP(y~x-1, coords=coords, starting=starting, method="latent", n.neighbors=10,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, n.omp.threads=1)

summary(m.s)
plot(apply(m.s$p.w.samples, 1, median), w)

m.r <- spNNGP(y~x-1, coords=coords, starting=starting, method="response", n.neighbors=10,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, n.omp.threads=1)

summary(m.r)

##Fit with some user defined neighbor ordering

##ord <- order(coords[,2]) ##y-axis 
ord <- order(coords[,1]+coords[,2]) ##x+y-axis
##ord <- sample(1:n, n) ##random

m.r.xy <- spNNGP(y~x-1, coords=coords, starting=starting, method="response", n.neighbors=10,
              tuning=tuning, priors=priors, cov.model=cov.model,
              ord=ord, return.neighbor.info=TRUE,
              n.samples=n.samples, n.omp.threads=1)

summary(m.r.xy)

\dontrun{
##Visualize the neighbor sets and ordering constraint
n.indx <- m.r.xy$neighbor.info$n.indx
ord <- m.r.xy$neighbor.info$ord

##This is how the data are ordered internally for model fitting
coords.ord <- coords[ord,]

for(i in 1:n){

    plot(coords.ord, cex=1, xlab="Easting", ylab="Northing")
    points(coords.ord[i,,drop=FALSE], col="blue", pch=19, cex=1)
    points(coords.ord[n.indx[[i]],,drop=FALSE], col="red", pch=19, cex=1)

    readline(prompt = "Pause. Press <Enter> to continue...")
}
}


}

\keyword{model}
