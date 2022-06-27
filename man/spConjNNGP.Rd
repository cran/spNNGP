\name{spConjNNGP}
\alias{spConjNNGP}
\title{Function for Fitting Univariate Bayesian Conjugate Spatial
  Regression Models}

\description{
  The function \code{spConjNNGP} fits Gaussian univariate Bayesian conjugate spatial
  regression models using Nearest Neighbor Gaussian Processes (NNGP).
}

\usage{
    spConjNNGP(formula, data = parent.frame(), coords, knots, n.neighbors = 15,
               theta.alpha, sigma.sq.IG, cov.model = "exponential",
               k.fold = 5, score.rule = "crps",
               X.0, coords.0, n.omp.threads = 1, search.type = "cb",
               ord, return.neighbor.info = TRUE, 
               neighbor.info, fit.rep = FALSE, n.samples, verbose = TRUE, ...)
}

\arguments{
  \item{formula}{a symbolic description of the regression model to be
    fit. See example below. }
  
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which
    \code{spConjNNGP} is called.  }
  
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix of the observation
    coordinates in \eqn{R^2}{R^2} (e.g., easting and northing), or if
    \code{data} is a data frame then \code{coords} can be a vector of
    length two comprising coordinate column names or indices. There can
    be no duplicate locations. }

  \item{knots}{an \eqn{r \times 2}{r x 2} matrix of the observation coordinates
    in \eqn{R^2}{R^2} (e.g., easting and northing). Adding the
    \code{knots} argument invokes SLGP, see Shin et al. (2019) below.}
  
  \item{n.neighbors}{number of neighbors used in the NNGP.}
  
  \item{theta.alpha}{a vector or matrix of parameter values for
    \code{phi}, \code{nu}, and \code{alpha}, where
    \eqn{\alpha=\tau^2/\sigma^2}{alpha=tau^2/sigma^2} and \code{nu} is
    only required if \code{cov.model="matern"}. A vector is passed to
    run the model using one set of parameters. The vector elements must
    be named and hold values for \code{phi}, \code{nu}, and
    \code{alpha}. If a matrix is passed, columns must be named and hold
    values for \code{phi}, \code{nu}, and \code{alpha}. Each row in the matrix defines a set of parameters for which the model will be run. }

  \item{sigma.sq.IG}{a vector of length two that holds the
    hyperparameters, \emph{shape} and \emph{scale} respectively, for the
    inverse-Gamma prior on \eqn{\sigma^2}{sigma^2}. }
  
  \item{cov.model}{a quoted keyword that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.  See below for details.}
  
  \item{k.fold}{specifies the number of
    \emph{k} folds for cross-validation. If \code{theta.alpha} is a
    vector then cross-validation is not performed and
    \code{k-fold} and \code{score.rule} are ignored. In \emph{k}-fold
    cross-validation, the data specified in \code{model} is randomly
    partitioned into \emph{k} equal sized subsamples. Of the \emph{k}
    subsamples, \emph{k}-1 subsamples are used to fit the model and the
    remaining \emph{k} samples are used for prediction. The
    cross-validation process is repeated \emph{k} times (the folds). Root
    mean squared prediction error (RMSPE) and continuous ranked
    probability score (CRPS; Gneiting and Raftery, 2007) rules are
    averaged over the \emph{k} fold prediction results and reported for
    the parameter sets defined by \code{theta.alpha}. The parameter set
    that yields the \emph{best} performance based on the scoring rule
    defined by \code{score.rule} is used to fit the final model that uses all
    the data and make predictions if \code{X.0} and \code{coords.0} are
    supplied.  Results from the \emph{k}-fold cross-validation are returned in
    the \code{k.fold.scores} matrix. }
  
  \item{score.rule}{a quoted keyword \code{"rmspe"} or \code{"crps"}
    that specifies the scoring rule used to select the \emph{best}
    parameter set, see argument definition for \code{k.fold} for more
      details.}
    
    \item{X.0}{the design matrix for prediction locations. An
      intercept should be provided in the first column if one is specified
      in \code{model}.}
    
    \item{coords.0}{the spatial coordinates corresponding to
      \code{X.0}. }
    
    \item{n.omp.threads}{a positive integer indicating
      the number of threads to use for SMP parallel processing. The package must
      be compiled for OpenMP support. For most Intel-based machines, we recommend setting
      \code{n.omp.threads} up to the number of hyperthreaded cores. Note, \code{n.omp.threads} > 1 might not
    work on some systems.}

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
      \code{order(coords[,1])}. This argument should typically be
      left blank. This argument will be ignored if the \code{neighbor.info}
      argument is used. }
    
    \item{return.neighbor.info}{if \code{TRUE}, a list called
      \code{neighbor.info} containing several data structures used for
      fitting the NNGP model is returned. If there is no change in input
      data or \code{n.neighbors}, this list can be passed to subsequent
      \code{spNNGP} calls via the \code{neighbor.info} argument to avoid
      the neighbor search, which can be time consuming if \eqn{n}{n} is
      large. In addition to the several cryptic data structures in
      \code{neighbor.info} there is a list called \code{n.indx} that is of
      length \eqn{n}{n}. The \emph{i}-th element in \code{n.indx}
      corresponds to the \emph{i}-th row in \code{coords[ord,]} and holds
      the vector of that location's nearest neighbor indices. This list can be useful for plotting the
      neighbor graph if desired. }
    
    \item{neighbor.info}{see the \code{return.neighbor.info} argument
      description above.}

    \item{fit.rep}{if \code{TRUE}, regression fitted and replicate data will be
      returned. The argument \code{n.samples} controls the number of
      fitted and replicated data samples. }

    \item{n.samples}{gives the number of posterior samples
      returned. Note, point and associated variance estimates for model
      parameters are not based on posterior samples. Only specify
      \code{n.samples} if you wish to generate samples from parameters'
      posteriors (this is an exact sampling algorithm). If
      \code{fit.rep} is \code{TRUE}, then \code{n.samples} also controls
      the number of fitted and replicated data samples.}
    
  \item{verbose}{if \code{TRUE}, model specification and progress is printed to the screen. Otherwise, nothing is printed to
    the screen.}
  
  \item{...}{currently no additional arguments.  }
}

\value{
  An object of class \code{NNGP} and \code{conjugate}, and, if \code{knots}
  is provided, \code{SLGP}. Among other elements, the object contains:
  
  \item{theta.alpha}{the input \code{theta.alpha} vector, or \emph{best} (according to the selected scoring
    rule) set of parameters in the \code{theta.alpha} matrix. All
    subsequent parameter estimates are based on this parameter set.}
  
  \item{beta.hat}{a matrix of regression coefficient estimates
    corresponding to the returned \code{theta.alpha}.}
  
  \item{beta.var}{\code{beta.hat} variance-covariance matrix.}
  
  \item{sigma.sq.hat}{estimate of \eqn{\sigma^2}{sigma^2} corresponding to the returned \code{theta.alpha}.}
  
  \item{sigma.sq.var}{\code{sigma.sq.hat} variance.}
  
  \item{k.fold.scores}{results from the k-fold cross-validation if
    \code{theta.alpha} is a matrix.}
  
  \item{y.0.hat}{prediction if \code{X.0} and \code{coords.0} are
    specified.}
  
  \item{y.0.var.hat}{\code{y.0.hat} variance.}
  
  \item{n.neighbors}{number of neighbors used in the NNGP.}
  
  \item{neighbor.info}{returned if \code{return.neighbor.info=TRUE} see the \code{return.neighbor.info} argument
      description above.}
      
  \item{run.time}{execution time for parameter estimation reported using
    \code{proc.time()}. This time does not include nearest neighbor
    search time for building the neighbor set.}
}

\references{
  Datta, A., S. Banerjee, A.O. Finley, and A.E. Gelfand. (2016)
  Hierarchical Nearest-Neighbor Gaussian process models for large
  geostatistical datasets. \emph{Journal of the American Statistical
  Association}, \doi{10.1080/01621459.2015.1044091}.

  Finley, A.O., A. Datta, S. Banerjee (2022) spNNGP R Package for
  Nearest Neighbor Gaussian Process Models. \emph{Journal of Statistical
  Software}, \doi{10.18637/jss.v103.i05}.

  Finley, A.O., A. Datta, B.D. Cook, D.C. Morton, H.E. Andersen, and
  S. Banerjee. (2019) Efficient algorithms for Bayesian Nearest Neighbor
  Gaussian Processes. \emph{Journal of Computational and Graphical
  Statistics}, \doi{10.1080/10618600.2018.1537924}.
  
  Gneiting, T and A.E. Raftery. (2007) Strictly proper scoring rules,
  prediction, and estimation. \emph{Journal of the American Statistical
    Association}, \doi{10.1198/016214506000001437}.

  Shirota, S., A.O. Finley, B.D. Cook, and S. Banerjee (2019) Conjugate Nearest Neighbor Gaussian Process models for efficient statistical interpolation of large spatial data. \url{https://arxiv.org/abs/1907.10109}.
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
n <- 2000
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

ho <- sample(1:n, 1000)

y.ho <- y[ho]
x.ho <- x[ho,,drop=FALSE]
w.ho <- w[ho]
coords.ho <- coords[ho,]

y <- y[-ho]
x <- x[-ho,,drop=FALSE]
w <- w[-ho,,drop=FALSE]
coords <- coords[-ho,]

##Fit a Conjugate NNGP model and predict for the holdout
sigma.sq.IG <- c(2, sigma.sq)

cov.model <- "exponential"

g <- 10
theta.alpha <- cbind(seq(phi,30,length.out=g), seq(tau.sq/sigma.sq,5,length.out=g))

colnames(theta.alpha) <- c("phi", "alpha")

m.c <- spConjNNGP(y~x-1, coords=coords, n.neighbors = 10,
                  X.0 = x.ho, coords.0 = coords.ho,
                  k.fold = 5, score.rule = "crps",
                  n.omp.threads = 1,
                  theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, cov.model = cov.model)

m.c


}


\keyword{model}
