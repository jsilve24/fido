context("degrees of freedom in shared model paths")


df_shared_logdet <- function(x) {
  as.numeric(determinant(x, logarithm = TRUE)$modulus)
}


df_shared_orthus_kernel <- function(eta, problem) {
  eta <- matrix(eta, problem$D - 1L, problem$N)
  latent <- rbind(eta, problem$Z)
  residual <- latent - problem$Theta %*% problem$X
  A <- diag(problem$N) +
    crossprod(problem$X, problem$Gamma %*% problem$X)
  XiN <- problem$Xi + residual %*% solve(A, t(residual))

  logits_kernel <- sum(problem$Y[seq_len(problem$D - 1L), ] * eta) -
    sum(colSums(problem$Y) * log(1 + colSums(exp(eta))))
  logits_kernel -
    (problem$upsilon + problem$N) * df_shared_logdet(XiN) / 2
}


test_that("orthus conditional df optimizes the full NIW-collapsed joint", {
  skip_if_not_installed("numDeriv")

  problem <- list(
    D = 3L,
    P = 2L,
    N = 4L,
    X = rbind(1, c(-0.8, -0.1, 0.4, 0.9)),
    Y = matrix(c(
      18, 10, 7,
      12, 15, 8,
      9, 13, 16,
      14, 8, 17
    ), nrow = 3L),
    Z = matrix(c(
      -0.5, 0.3,
      -0.1, -0.4,
      0.4, 0.2,
      0.8, 0.6
    ), nrow = 2L),
    Theta = matrix(c(
      0.10, -0.06,
      -0.08, 0.04,
      0.15, 0.03,
      -0.12, 0.07
    ), nrow = 4L, byrow = TRUE),
    Gamma = matrix(c(0.55, 0.08, 0.08, 0.35), 2L, 2L),
    Xi = matrix(c(
      1.00, 0.14, 0.18, -0.08,
      0.14, 0.90, -0.10, 0.16,
      0.18, -0.10, 1.20, 0.20,
      -0.08, 0.16, 0.20, 1.10
    ), 4L, 4L),
    upsilon = 9
  )
  init <- log(problem$Y[seq_len(problem$D - 1L), , drop = FALSE] /
                problem$Y[problem$D, ])

  fit <- orthus(
    problem$Y, problem$Z, problem$X,
    upsilon = problem$upsilon, Theta = problem$Theta,
    Gamma = problem$Gamma, Xi = problem$Xi, init = init,
    n_samples = 0, calcGradHess = FALSE, ret_mean = TRUE,
    max_iter = 5000, eps_f = 1e-13, eps_g = 1e-8,
    ncores = 1, seed = 1401, use_names = FALSE
  )
  map <- fit$Eta[, , 1]
  full_gradient <- numDeriv::grad(
    function(x) df_shared_orthus_kernel(x, problem), c(map)
  )

  # For the full row dimension D-1+P, IW df first becomes
  # upsilon-(D-1+P)+1, then conditioning on P rows adds P. Consequently
  # the optimizer's determinant exponent is (upsilon+N)/2.
  expect_lt(max(abs(full_gradient)), 2e-5)
  expect_equal(fit$upsilon, problem$upsilon)
})


test_that("maltipoo converts public IW df before its collapsed optimizer", {
  D <- 3L
  N <- 5L
  Q <- 2L
  X <- rbind(1, seq(-0.9, 0.8, length.out = N))
  Y <- matrix(c(
    16, 11, 8,
    13, 14, 9,
    10, 16, 12,
    9, 12, 18,
    14, 9, 15
  ), nrow = D)
  Theta <- matrix(c(0.08, -0.04, -0.06, 0.05), D - 1L, Q,
                  byrow = TRUE)
  Xi <- matrix(c(1.0, 0.16, 0.16, 0.85), D - 1L, D - 1L)
  U <- rbind(
    diag(Q),
    matrix(c(0.65, 0.09, 0.09, 0.45), Q, Q)
  )
  ellinit <- log(c(0.35, 0.20))
  init <- log(Y[seq_len(D - 1L), , drop = FALSE] / Y[D, ])
  upsilon <- 9
  matrix_t_df <- upsilon - (D - 1L) + 1

  fit <- fido:::maltipoo(
    Y, X, upsilon = upsilon, Theta = Theta, U = U, Xi = Xi,
    init = init, ellinit = ellinit, n_samples = 0,
    calcGradHess = FALSE, ret_mean = TRUE, max_iter = 2500,
    eps_f = 1e-6, eps_g = 1e-3, seed = 1501, use_names = FALSE
  )
  low <- fido:::optimMaltipooCollapsed(
    Y, matrix_t_df, Theta, X, solve(Xi), U, init, ellinit,
    n_samples = 0, calcGradHess = FALSE, max_iter = 2500,
    eps_f = 1e-6, eps_g = 1e-3
  )

  expect_equal(fit$Eta[, , 1], low$Pars, tolerance = 1e-12)
  expect_equal(fit$VCScale, low$VCScale, tolerance = 1e-12)
  expect_equal(fit$upsilon, upsilon)
})


test_that("a serialized pibble fit retains IW df when refitted", {
  D <- 3L
  N <- 4L
  X <- rbind(1, c(-0.7, -0.1, 0.5, 0.9))
  Y <- matrix(c(
    15, 10, 8,
    12, 14, 9,
    9, 15, 13,
    11, 8, 17
  ), nrow = D)
  Theta <- matrix(c(0.06, -0.03, -0.05, 0.04), D - 1L, 2L,
                  byrow = TRUE)
  Gamma <- matrix(c(0.5, 0.07, 0.07, 0.3), 2L, 2L)
  Xi <- matrix(c(1.0, 0.12, 0.12, 0.9), D - 1L, D - 1L)
  init <- log(Y[seq_len(D - 1L), , drop = FALSE] / Y[D, ])
  upsilon <- 8
  common <- list(
    n_samples = 0, calcGradHess = FALSE, ret_mean = TRUE,
    max_iter = 2500, eps_f = 1e-12, eps_g = 1e-7,
    ncores = 1, seed = 1601, use_names = FALSE
  )
  original <- do.call(pibble, c(list(
    Y = Y, X = X, upsilon = upsilon, Theta = Theta,
    Gamma = Gamma, Xi = Xi, init = init
  ), common))
  path <- tempfile(fileext = ".rds")
  saveRDS(original, path)
  on.exit(unlink(path), add = TRUE)
  restored <- readRDS(path)

  refitted <- do.call(refit, c(list(
    m = restored, init = init
  ), common))
  fresh <- do.call(pibble, c(list(
    Y = Y, X = X, upsilon = upsilon, Theta = Theta,
    Gamma = Gamma, Xi = Xi, init = init
  ), common))

  expect_equal(refitted$Eta, fresh$Eta, tolerance = 1e-12)
  expect_equal(refitted$Lambda, fresh$Lambda, tolerance = 1e-12)
  expect_equal(refitted$Sigma, fresh$Sigma, tolerance = 1e-12)
  expect_equal(refitted$upsilon, upsilon)
})
