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
  # adding seed internally because pibble_sim causing tests to fail
  set.seed(1234)
  sim <- pibble_sim()
  Gamma <- function(X) SE(X)
  Theta <- function(X) matrix(0, nrow(sim$Y)-1, ncol(X))
  fit <- basset(sim$Y, sim$X, Gamma = Gamma, Theta = Theta, n_samples = 10000, seed = 123)
  fit.new <- basset(sim$Y, sim$X, Gamma = list(Gamma), Theta = list(Theta), n_samples = 10000, seed = 123)
  
  ## Estimates match
  fit.lam <- apply(fit$Lambda, c(1,2), mean)
  fit.new.lam <- apply(fit.new$Lambda[[1]], c(1,2), mean)
  
  expect_true(mean(abs(fit.lam-fit.new.lam)) < .1)
})

test_that("basset can have multiple GP components",{
  # adding seed internally because pibble_sim causing tests to fail
  set.seed(1234)
  sim <- pibble_sim()
  Gamma <- function(X) SE(X)
  Theta <- function(X) matrix(0, nrow(sim$Y)-1, ncol(X))
  fit.new <- basset(sim$Y, sim$X, Gamma = list(Gamma,Gamma), Theta = list(Theta,Theta), n_samples = 10000)
  expect_true(TRUE)
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
  
  expect_true(mean(abs(mod.Eta - mod.pib.Eta)) < 0.05)
  
  mod.Lam <- apply(mod$Lambda[[1]], c(1,2), mean)
  mod.pib.Lam <- apply(mod.pib$Lambda, c(1,2), mean)
  
  expect_true(mean(abs(mod.Lam - mod.pib.Lam)) < 0.05)
  
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

test_that("testing that r2 works",{
  # adding seed internally because pibble_sim causing tests to fail
  set.seed(1234)
  sim <- pibble_sim()
  
  Theta <- list(matrix(0, nrow = sim$D-1, ncol = sim$Q), function(X) matrix(0, nrow = sim$D-1, ncol = ncol(X)))
  Gamma <- list(diag(sim$Q), function(X) SE(X, sigma=5, rho=10))
  
  ##Running basset
  mod <- basset(sim$Y, sim$X, sim$upsilon, Theta,Gamma, sim$Xi)
  expect_equal(mean(r2.bassetfit(mod)),  mean(r2.bassetfit(mod, components = c(1,2))), tolerance = 1e-3)
  expect_equal(mean(r2.bassetfit(mod, components = 1)), mean(r2.bassetfit(mod, components = 1, covariates = 1)) + mean(r2.bassetfit(mod, components = 1, covariates=2)), tolerance = 1e-3)
  #expect_equal(mean(r2.bassetfit(mod)), mean(r2.bassetfit(mod, components = 1)) + mean(r2.bassetfit(mod, components = 2)), tolerance = 1e-3)
  
})


test_that("transforms work with new basset", {
  # adding seed internally because pibble_sim causing tests to fail
  set.seed(123)
  sim <- pibble_sim()
  Gamma <- function(X) SE(X)
  Theta <- function(X) matrix(0, nrow(sim$Y)-1, ncol(X))
  fit.new <- basset(sim$Y, sim$X, Gamma = list(Gamma), Theta = list(Theta))
  clr.fit <- to_clr(fit.new)
  preds <- predict(clr.fit)
  alr.fit <- to_alr(fit.new, 1)
  preds <- predict(alr.fit)
  ilr.fit <- to_ilr(fit.new)
  preds <- predict(ilr.fit)
  expect_true(TRUE)
})

#' @param eta as an (D-1 x N x iter) array
uncollapse <- function(eta, X, upsilon, Theta, Xi, Gamma, GammaComb, Sigma){
  d <- dim(eta)
  iter <- as.integer(d[3])
  N <- as.integer(d[2])
  D <- as.integer(d[1] + 1)
  Q <- as.integer(nrow(Gamma))
  Lambda <- array(0, c(D-1, Q, iter))
  
  GammaInv <- solve(Gamma)
  GammaCombInv <- solve(GammaComb)
  GammaN <- solve(X %*% GammaCombInv %*% t(X) + GammaInv)
  for (i in 1:iter){
    LambdaN <- (eta[,,i] %*% GammaCombInv %*% t(X) + Theta %*% GammaInv) %*% GammaN
    Z <- matrix(rnorm((D-1)*Q), D-1, Q)
    Lambda[,,i] <- LambdaN + t(chol(Sigma[,,i]))%*%Z%*%chol(GammaN)
  }
  return(list(Lambda=Lambda, Sigma=Sigma))
}



test_that("basset c++ matches R implementation", {
  # adding seed internally because pibble_sim causing tests to fail
  set.seed(1234)
  sim <- pibble_sim()
  fit <- pibble(sim$Y, sim$X, Gamma = sim$Gamma, Theta = sim$Theta, Xi = sim$Xi, upsilon = sim$upsilon, seed = 1, n_samples = 5000)
  
  # Now check uncollapsing for Lambda with C++
  fit.test <- fido:::uncollapsePibble_sigmaKnown(fit$Eta, sim$X, sim$Theta, sim$Gamma,
                                                 GammaComb = diag(sim$N), fit$Xi, sigma = fit$Sigma,
                                                 upsilon= sim$upsilon, seed = 12345)
  
  r.fit <- uncollapse(fit$Eta, sim$X, sim$upsilon, sim$Theta, sim$Xi, 
                      sim$Gamma, diag(sim$N), fit$Sigma)
  
  Lambda.fit <- apply(fit.test$Lambda, c(1,2), mean)
  Lambda.r <- apply(r.fit$Lambda, c(1,2), mean)
  
  expect_true(mean(abs(Lambda.fit - Lambda.r)) < 0.05)
  
  ##Now with arbitrary GammaComb
  # Now check uncollapsing for Lambda with C++
  GammaComb <- matrix(.5, nrow = sim$N, ncol = sim$N)
  diag(GammaComb) <- 1
  fit.test <- fido:::uncollapsePibble_sigmaKnown(fit$Eta, sim$X, sim$Theta, sim$Gamma,
                                                 GammaComb = GammaComb, fit$Xi, sigma = fit$Sigma,
                                                 upsilon= sim$upsilon, seed = 12345)
  
  r.fit <- uncollapse(fit$Eta, sim$X, sim$upsilon, sim$Theta, sim$Xi, 
                      sim$Gamma, GammaComb, fit$Sigma)
  
  Lambda.fit <- apply(fit.test$Lambda, c(1,2), mean)
  Lambda.r <- apply(r.fit$Lambda, c(1,2), mean)
  
  expect_true(mean(abs(Lambda.fit - Lambda.r)) < 0.05)
})


test_that("new collapse sampler matches old when inputs match",{
  sim <- pibble_sim()
  fit <- pibble(sim$Y, sim$X, Gamma = sim$Gamma, Theta = sim$Theta, Xi = sim$Xi, upsilon = sim$upsilon, seed = 1, n_samples = 2000)
  Lambda.fit <- apply(fit$Lambda, c(1,2), mean)
  fit.test <- fido:::uncollapsePibble_sigmaKnown(fit$Eta, sim$X, sim$Theta, sim$Gamma,
                                     GammaComb = diag(sim$N), fit$Xi, sigma = fit$Sigma,
                                     upsilon= sim$upsilon, seed = 1)
  Lambda.fit.test <- apply(fit.test$Lambda, c(1,2), mean)
  
  expect_true(mean(abs(Lambda.fit - Lambda.fit.test)) < 0.05)
})



