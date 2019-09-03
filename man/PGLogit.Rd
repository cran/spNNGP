\name{PGLogit}
\alias{PGLogit}
\title{Function for Fitting Logistic Models using Polya-Gamma Latent Variables}

\description{
  The function \code{PGLogit} fits logistic models to binomial data using Polya-Gamma latent variables.
}

\usage{
PGLogit(formula, weights = 1, data = parent.frame(), n.samples,
        n.omp.threads = 1, fit.rep = FALSE, sub.sample, verbose = TRUE, ...)
}

\arguments{
  \item{formula}{a symbolic description of the regression model to be
    fit. See example below. }
  
  \item{weights}{specifies the number of trials for each observation. The
    default is 1 trial for each observation. Valid arguments are a
    scalar value that specifies the number of trials used if all
    observations have the same number of trials, and a vector of length \eqn{n}{n} that specifies the number of trials for
    each observation when there are differences in the number of trials. }
  
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which \code{PGLogit} is called.  }
  
  \item{n.samples}{the number of posterior samples to collect.}
  
  \item{n.omp.threads}{a positive integer indicating
    the number of threads to use for SMP parallel processing. The package must
    be compiled for OpenMP support. For most Intel-based machines, we
    recommend setting \code{n.omp.threads} up to the number of
    hyperthreaded cores. Note, \code{n.omp.threads} > 1 might not
    work on all systems.}

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
  
  \item{...}{currently no additional arguments.  }
}

\value{
  An object of class \code{PGLogit} which is a list comprising:
  
  \item{p.beta.samples}{a \code{coda} object of posterior samples for
    the regression coefficients. }

  \item{y.hat.samples}{if \code{fit.rep=TRUE}, regression fitted values from posterior samples
    specified using \code{sub.sample}. }

  \item{y.hat.quants}{if \code{fit.rep=TRUE}, 0.5, 0.025, and 0.975 quantiles of the \code{y.hat.samples}. }

  \item{y.rep.samples}{if \code{fit.rep=TRUE}, replicated outcome from posterior samples
    specified using \code{sub.sample}. }
  
  \item{y.rep.quants}{if \code{fit.rep=TRUE}, 0.5, 0.025, and 0.975
    quantiles of the \code{y.rep.samples}. }

  \item{s.indx}{if \code{fit.rep=TRUE}, the subset index specified with \code{sub.sample}.}
  
  \item{run.time}{MCMC sampler execution time reported using \code{proc.time()}.}
  
  The return object will include additional objects used for subsequent
  prediction and/or model fit evaluation.
}

\references{
  Polson, N.G., J.G. Scott, and J. Windle. (2013) Bayesian Inference for
  Logistic Models Using Polya-Gamma Latent Variables. Journal of the
  American Statistical Association, 108:1339-1349.

  
}

\note{
 Some of the underlying code used for generating random number from the
 Polya-Gamma distribution is taken from the \pkg{pgdraw} package written
 by Daniel F. Schmidt and Enes Makalic. Their code implements Algorithm
 6 in PhD thesis of Jesse Bennett Windle (2013) \url{https://repositories.lib.utexas.edu/handle/2152/21842}.
}


\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Abhirup Datta \email{abhidatta@jhu.edu}, \cr
  Sudipto Banerjee \email{sudipto@ucla.edu}
}

\examples{

##Generate binary data
set.seed(1)
n <- 100

x <- cbind(1, rnorm(n), runif(n,0,1))
beta <- c(0.1,-5, 5)
p <- 1/(1+exp(-(x\%*\%beta)))

##Assume 5 trials per outcome
weights <- rep(5, n)

y <- rbinom(n, size=weights, prob=p)

m <- PGLogit(y~x-1, weights = rep(5, n), n.samples = 1000)

summary(m)

}

\keyword{model}