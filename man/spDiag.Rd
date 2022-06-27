\name{spDiag}
\alias{spDiag}          
\title{Model fit diagnostics}

\description{
  The function \code{spDiag} calculates measurements of model fit for
  objects of class \code{NNGP} and \code{PGLogit}.
}

\usage{
  spDiag(object, sub.sample, ...)
}

\arguments{
  \item{object}{an object of class \code{NNGP} or \code{PGLogit}.}
  \item{sub.sample}{an optional list that specifies the samples to included in
  the computations. Valid tags are \code{start},
  \code{end}, and \code{thin}. Given the value associated with the tags,
  the sample subset is selected using \code{seq(as.integer(start),
    as.integer(end), by=as.integer(thin))}. The default values are
  \code{start=floor(0.5*n.samples)}, \code{end=n.samples} and
  \code{thin=1}. If \code{sub.samples} is not specified, then it is
  taken from \code{object}, or, if not aviable in \code{object} the
  default values of \code{start}, \code{end}, and \code{thin} are
  used. Note, if the \code{object} is a \code{NNGP} \code{response} model
  and \code{n} is large, then computing the replicated data needed for
  \code{GPD} and \code{GRS} can take a long time. }
  \item{...}{currently no additional arguments.  }
}

\value{
  A list with the following tags:
  \item{DIC}{a data frame holding Deviance information criterion (DIC) and associated values. Values in \code{DIC} include
    \code{DIC} the criterion (lower is better), \code{D} a goodness of fit, and \code{pD} the effective
    number of parameters, see Spiegelhalter et al. (2002) for details. }
  \item{GPD}{a data frame holding D=G+P and associated values. Values in
    \code{GPD} include \code{G} a goodness of fit, \code{P} a penalty term, and \code{D} the
    criterion (lower is better), see Gelfand
    and Ghosh (1998) for details. }
  \item{GRS}{a scoring rule, see Equation 27 in Gneiting and Raftery
    (2007) for details.}
  \item{WAIC}{a data frame hold Watanabe-Akaike information criteria (WAIC) and associated values. Values in
    \code{WAIC} include \code{LPPD} log pointwise predictive density,
    \code{P.1} penalty term defined in unnumbered equation above Equation
    (11) in Gelman et al. (2014), \code{P.2} an alternative penalty
    term defined in Equation (11), and the criteria \code{WAIC.1} and
    \code{WAIC.2} (lower is better) computed using \code{P.1} and
    \code{P.2}, respectively.}
  \item{y.rep.samples}{if \code{y.rep.samples} in \code{object} were not
    used (or not available), then the newly computed \code{y.rep.samples} is returned.}
  \item{y.fit.samples}{if \code{y.fit.samples} in \code{object} were not
    used (or not available), then the newly computed \code{y.fit.samples} is returned.}
  \item{s.indx}{the index of samples used for the computations.}
}


\references{
  Finley, A.O., A. Datta, S. Banerjee (2022) spNNGP R Package for
  Nearest Neighbor Gaussian Process Models. \emph{Journal of Statistical
  Software}, \doi{10.18637/jss.v103.i05}.
  
  Gelfand A.E. and Ghosh, S.K. (1998). Model choice: a minimum posterior
  predictive loss approach. \emph{Biometrika}, 85:1-11.
  
  Gelman, A., Hwang, J., and  Vehtari, A. (2014). Understanding predictive
  information criteria for Bayesian models. \emph{Statistics and
    Computing}, 24:997-1016.

  Gneiting, T. and Raftery, A.E. (2007). Strictly proper scoring rules, prediction, and estimation.
  \emph{Journal of the American Statistical Association}, 102:359-378.

  Spiegelhalter, D.J., Best, N.G., Carlin, B.P., van der Linde, A. (2002). Bayesian measures of model complexity and fit (with discussion). \emph{Journal of the Royal Statistical Society, Series B.}, 64:583-639.
    
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{sudipto@ucla.edu}
}


\keyword{utilities}