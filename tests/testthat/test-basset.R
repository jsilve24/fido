context("test-basset")

set.seed(569)

test_that("basset and predict.bassetfit run", {
  sim <- pibble_sim()
  Gamma <- function(X) SE(X)
  Theta <- function(X) matrix(0, nrow(sim$Y)-1, ncol(X))
  fit <- basset(sim$Y, sim$X, Gamma = Gamma, Theta = Theta)
  foo <- predict(fit, matrix(c(1,2)))
  expect_true(TRUE)
})

test_that("basset matches old if list is used as input",{
  sim <- pibble_sim()
  Gamma <- function(X) SE(X)
  Theta <- function(X) matrix(0, nrow(sim$Y)-1, ncol(X))
  fit <- basset(sim$Y, sim$X, Gamma = Gamma, Theta = Theta, n_samples = 10000)
  fit.new <- basset(sim$Y, sim$X, Gamma = list(Gamma), Theta = list(Theta), n_samples = 10000)
  
  ## Estimates match
  fit.lam <- apply(fit$Lambda, c(1,2), mean)
  fit.new.lam <- apply(fit.new$Lambda[[1]], c(1,2), mean)
  
  expect_equal(fit.lam[1,], fit.new.lam[1,], tolerance = .1)
  expect_equal(fit.lam[2,], fit.new.lam[2,], tolerance = .1)
  expect_equal(fit.lam[3,], fit.new.lam[3,], tolerance = .1)
})

test_that("basset matches fido when only linear terms are used", {
  sim <- pibble_sim()

  Theta <- matrix(0, nrow = sim$D-1, ncol = sim$Q)
  Gamma <- diag(sim$Q)
  
  ##Running basset
  mod <- basset(sim$Y, sim$X, sim$upsilon, list(Theta), list(Gamma), sim$Xi)
  
  mod.pib <- pibble(sim$Y, sim$X, sim$upsilon, Theta, Gamma, sim$Xi)
  
  mod.Eta <- apply(mod$Eta, c(1,2), mean)
  mod.pib.Eta <- apply(mod.pib$Eta, c(1,2), mean)
  
  expect_equal(mod.Eta[1,], mod.pib.Eta[1,], tol = 1e-1)
  expect_equal(mod.Eta[2,], mod.pib.Eta[2,], tol = 1e-1)
  expect_equal(mod.Eta[3,], mod.pib.Eta[3,], tol = 1e-1)
  
  mod.Lam <- apply(mod$Lambda[[1]], c(1,2), mean)
  mod.pib.Lam <- apply(mod.pib$Lambda, c(1,2), mean)
  
  expect_equal(mod.Lam[,1], unname(mod.pib.Lam[,1]), tol = 5e-1)
  expect_equal(mod.Lam[,2], unname(mod.pib.Lam[,2]), tol = 5e-1)
  
})

test_that("testing that predict works",{
  sim <- pibble_sim()
  
  Theta <- function(X) matrix(0, nrow = sim$D-1, ncol = ncol(X))
  Gamma <- function(X) diag(ncol(X))
  
  ##Running basset
  mod <- basset(sim$Y, sim$X, sim$upsilon, Theta,Gamma, sim$Xi)
  foo <- predict(mod)
  
  ##Running basset
  mod <- basset(sim$Y, sim$X, sim$upsilon, list(Theta), list(Gamma), sim$Xi)
  foo <- predict(mod)
  
  expect_true(TRUE)
  
})
