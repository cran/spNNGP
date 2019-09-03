\name{summary.PGLogit}

\alias{summary.PGLogit}
\alias{fitted.PGLogit}
\alias{residuals.PGLogit}
\alias{print.PGLogit}

\title{Methods for PGLogit Object}

\description{
  Methods for extracting information from fitted \code{PGLogit} model.
}

\usage{
\method{summary}{PGLogit}(object, sub.sample, quantiles
=c(0.025, 0.25, 0.5, 0.75, 0.975), digits = max(3L, getOption("digits")
- 3L), \dots)
\method{print}{PGLogit}(x, \dots)
}

\arguments{
  \item{object, x}{object of class \code{PGLogit}.}
  \item{sub.sample}{an optional list that specifies the samples to included in
  the summary or composition sampling. Valid tags are \code{start},
  \code{end}, and \code{thin}. Given the value associated with the tags,
  the sample subset is selected using \code{seq(as.integer(start),
    as.integer(end), by=as.integer(thin))}. The default values are
  \code{start=floor(0.5*n.samples)}, \code{end=n.samples} and \code{thin=1}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{\dots}{currently not used.}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{PGLogit}, including methods to the generic functions \code{\link[base]{print}} and \code{\link[base]{summary}}. 
}

\keyword{model}
