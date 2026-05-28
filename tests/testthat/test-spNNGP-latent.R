make_gaussian_data <- function(n = 40) {
  coords <- cbind(runif(n), runif(n))
  x <- rnorm(n)
  y <- 1 + 0.5 * x + rnorm(n)

  list(
    y = y,
    x = x,
    coords = coords,
    starting = list(beta = c(0, 0), sigma.sq = 1, tau.sq = 1, phi = 3),
    tuning = list(phi = 0.1),
    priors = list(
      sigma.sq.IG = c(2, 1),
      tau.sq.IG = c(2, 1),
      phi.Unif = c(0.1, 30)
    )
  )
}

make_binomial_data <- function(n = 40) {
  coords <- cbind(runif(n), runif(n))
  x <- rnorm(n)
  y <- rbinom(n, 1, plogis(-0.2 + 0.8 * x))

  list(
    y = y,
    x = x,
    coords = coords,
    starting = list(beta = c(0, 0), sigma.sq = 1, phi = 3),
    tuning = list(phi = 0.1),
    priors = list(
      sigma.sq.IG = c(2, 1),
      phi.Unif = c(0.1, 30)
    )
  )
}

test_that("spNNGP validates neighbor and thread counts before native calls", {
  set.seed(1)
  dat <- make_gaussian_data(n = 10)

  expect_error(
    spNNGP(
      dat$y ~ dat$x,
      coords = dat$coords,
      method = "latent",
      n.neighbors = 15,
      starting = dat$starting,
      tuning = dat$tuning,
      priors = dat$priors,
      n.samples = 5,
      n.omp.threads = 1,
      verbose = FALSE
    ),
    "n.neighbors must be between 1 and n-1"
  )

  expect_error(
    spNNGP(
      dat$y ~ dat$x,
      coords = dat$coords,
      method = "latent",
      n.neighbors = 5,
      starting = dat$starting,
      tuning = dat$tuning,
      priors = dat$priors,
      n.samples = 5,
      n.omp.threads = 0,
      verbose = FALSE
    ),
    "n.omp.threads must be a positive integer"
  )
})

test_that("spNNGP uses the default covariance model", {
  set.seed(2)
  dat <- make_gaussian_data(n = 30)

  fit <- spNNGP(
    dat$y ~ dat$x,
    coords = dat$coords,
    method = "latent",
    n.neighbors = 5,
    starting = dat$starting,
    tuning = dat$tuning,
    priors = dat$priors,
    n.samples = 5,
    n.omp.threads = 1,
    verbose = FALSE
  )

  expect_equal(fit$cov.model, "exponential")
  expect_equal(dim(fit$p.w.samples), c(30L, 5L))
})

test_that("latent Gaussian and binomial samplers produce expected sample shapes", {
  set.seed(3)
  gdat <- make_gaussian_data(n = 35)
  gfit <- spNNGP(
    gdat$y ~ gdat$x,
    coords = gdat$coords,
    method = "latent",
    n.neighbors = 6,
    starting = gdat$starting,
    tuning = gdat$tuning,
    priors = gdat$priors,
    n.samples = 6,
    n.omp.threads = 2,
    verbose = FALSE
  )

  expect_equal(dim(gfit$p.w.samples), c(35L, 6L))
  expect_equal(colnames(gfit$p.theta.samples), c("sigma.sq", "tau.sq", "phi"))

  set.seed(4)
  bdat <- make_binomial_data(n = 35)
  bfit <- spNNGP(
    bdat$y ~ bdat$x,
    family = "binomial",
    coords = bdat$coords,
    method = "latent",
    n.neighbors = 6,
    starting = bdat$starting,
    tuning = bdat$tuning,
    priors = bdat$priors,
    n.samples = 6,
    n.omp.threads = 2,
    verbose = FALSE
  )

  expect_equal(dim(bfit$p.w.samples), c(35L, 6L))
  expect_equal(colnames(bfit$p.theta.samples), c("sigma.sq", "phi"))
})
