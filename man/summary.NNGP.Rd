\name{summary.NNGP}

\alias{summary.NNGP}
\alias{fitted.NNGP}
\alias{residuals.NNGP}
\alias{print.NNGP}
\alias{print.predict.NNGP}

\title{Methods for NNGP and Derived Objects}

\description{
  Methods for extracting information from fitted NNGP model of class
  \code{NNGP} and \code{predict.NNGP} objects from \code{predict}.
}

\usage{
\method{summary}{NNGP}(object, sub.sample, quantiles = c(0.025, 0.25,
0.5, 0.75, 0.975), digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{NNGP}(x, \dots)
\method{print}{predict.NNGP}(x, \dots)
}

\arguments{
  \item{object, x}{object of class \code{NNGP} or \code{predict.NNGP}.}
  \item{sub.sample}{an optional list that specifies the samples to included in
    the summary or composition sampling. Valid tags are \code{start},
    \code{end}, and \code{thin}. Given the value associated with the tags,
    the sample subset is selected using \code{seq(as.integer(start),
      as.integer(end), by=as.integer(thin))}. The default values are
    \code{start=floor(0.5*n.samples)}, \code{end=n.samples} and \code{thin=1}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report in summary.}
  \item{\dots}{currently not used.}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{NNGP} and prediction object of class \code{predict.NNGP}, including methods to the generic functions \code{\link[base]{print}} and \code{\link[base]{summary}}. 
}

\keyword{model}