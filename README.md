# spNNGP

`spNNGP` fits univariate Bayesian spatial regression models for large
datasets using nearest neighbor Gaussian processes.

The package supports response and latent NNGP models for Gaussian
outcomes, latent NNGP models for binomial outcomes, and conjugate NNGP
models.

## Installation

Install the development version from GitHub after the repository is
published:

```r
remotes::install_github("finleya/spNNGP")
```

## Example

```r
library(spNNGP)

set.seed(1)
n <- 100
coords <- cbind(runif(n), runif(n))
x <- rnorm(n)
y <- 1 + 0.5 * x + rnorm(n)

fit <- spNNGP(
  y ~ x,
  coords = coords,
  method = "latent",
  n.neighbors = 10,
  starting = list(beta = c(0, 0), sigma.sq = 1, tau.sq = 1, phi = 3),
  tuning = list(phi = 0.1),
  priors = list(
    sigma.sq.IG = c(2, 1),
    tau.sq.IG = c(2, 1),
    phi.Unif = c(0.1, 30)
  ),
  n.samples = 100,
  n.omp.threads = 1,
  verbose = FALSE
)

summary(fit)
```

## References

Finley, A. O., Datta, A., and Banerjee, S. (2022). spNNGP R Package for
Nearest Neighbor Gaussian Process Models. *Journal of Statistical
Software*, 103(5). doi:10.18637/jss.v103.i05.

Finley, A. O., Datta, A., Cook, B. D., Morton, D. C., Andersen, H. E.,
and Banerjee, S. (2019). Efficient algorithms for Bayesian nearest
neighbor Gaussian processes. *Journal of Computational and Graphical
Statistics*, 28(2), 401-414. doi:10.1080/10618600.2018.1537924.

Datta, A., Banerjee, S., Finley, A. O., and Gelfand, A. E. (2016).
Hierarchical nearest-neighbor Gaussian process models for large
geostatistical datasets. *Journal of the American Statistical
Association*, 111(514), 800-812. doi:10.1080/01621459.2015.1044091.
