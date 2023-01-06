context("test-R2.R ")
set.seed(56468541)

test_that("r2 outputs correct dimension for pibblefit", {
  sim <- pibble_sim(true_priors=TRUE)
  fit <- pibble(sim$Y, sim$X, sim$upsilon, sim$Theta, sim$Gamma, sim$Xi)
  rsquared <- r2(fit)
  expect_equal(dim(fit$Eta)[3], length(rsquared))
})

test_that("r2 outputs correct dimension for bassetfit", {
  sim <- pibble_sim()
  Gamma <- function(X) SE(X)
  Theta <- function(X) matrix(0, nrow(sim$Y)-1, ncol(X))
  fit <- basset(sim$Y, sim$X, Gamma = Gamma, Theta = Theta)
  rsquared <- r2(fit)
  expect_equal(dim(fit$Eta)[3], length(rsquared))
})

