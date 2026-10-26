context("inverse-Wishart and matrix-t degrees of freedom")


# These helpers deliberately derive the collapsed density from the NIW model,
# rather than reproducing the determinant-only matrix-t code under test.
log_multigamma_iw <- function(a, p) {
  p * (p - 1) * log(pi) / 4 +
    sum(lgamma(a + (1 - seq_len(p)) / 2))
}

logdet_iw <- function(x) {
  as.numeric(determinant(x, logarithm = TRUE)$modulus)
}

iw_problem <- function(P, N) {
  D <- P + 1L
  X <- rbind(1, seq(-0.7, 0.8, length.out = N))
  Theta <- outer(seq_len(P), c(0.08, -0.05))
  Gamma <- matrix(c(0.7, 0.12, 0.12, 0.45), 2, 2)
  Xi <- diag(seq(0.8, 1.2, length.out = P), nrow = P, ncol = P)
  if (P > 1L) Xi[row(Xi) != col(Xi)] <- 0.08
  Y <- outer(seq_len(D), seq_len(N), function(d, n) d + 2L * n)
  eta <- outer(seq_len(P), seq_len(N), function(d, n) 0.09 * d - 0.06 * n)
  list(
    P = P, D = D, N = N, X = X, Theta = Theta, Gamma = Gamma,
    Xi = Xi, Y = Y, eta = eta, v = P + 6.75
  )
}

log_collapsed_niw <- function(eta, z) {
  eta <- matrix(eta, z$P, z$N)
  A <- diag(z$N) + crossprod(z$X, z$Gamma %*% z$X)
  E <- eta - z$Theta %*% z$X
  XiN <- z$Xi + E %*% solve(A, t(E))

  log_prior <- log_multigamma_iw((z$v + z$N) / 2, z$P) -
    log_multigamma_iw(z$v / 2, z$P) +
    z$v * logdet_iw(z$Xi) / 2 -
    z$P * z$N * log(pi) / 2 -
    z$P * logdet_iw(A) / 2 -
    (z$v + z$N) * logdet_iw(XiN) / 2

  log_mult <- sum(vapply(seq_len(z$N), function(j) {
    logits <- c(eta[, j], 0)
    log_norm <- max(logits) + log(sum(exp(logits - max(logits))))
    lgamma(sum(z$Y[, j]) + 1) - sum(lgamma(z$Y[, j] + 1)) +
      sum(z$Y[, j] * (logits - log_norm))
  }, numeric(1)))

  log_prior + log_mult
}

lowlevel_log_kernel <- function(eta, z, matrix_t_df, sylv = FALSE) {
  AInv <- solve(diag(z$N) + crossprod(z$X, z$Gamma %*% z$X))
  loglikPibbleCollapsed(
    z$Y, matrix_t_df, z$Theta %*% z$X, solve(z$Xi), AInv,
    eta, sylv = sylv
  )
}

posterior_moments_iw <- function(eta, X, Theta, Gamma, Xi, v) {
  GammaInv <- solve(Gamma)
  GammaN <- solve(GammaInv + tcrossprod(X))
  LambdaN <- (eta %*% t(X) + Theta %*% GammaInv) %*% GammaN
  residual <- eta - LambdaN %*% X
  delta <- LambdaN - Theta
  XiN <- Xi + tcrossprod(residual) + delta %*% GammaInv %*% t(delta)
  P <- nrow(eta)
  list(
    Lambda = LambdaN,
    Sigma = XiN / (v + ncol(X) - P - 1),
    precision = (v + ncol(X)) * solve(XiN),
    Gamma = GammaN,
    Xi = XiN
  )
}

posterior_mean_arrays_iw <- function(eta, X, Theta, Gamma, Xi, v) {
  dimensions <- dim(eta)
  lambda <- array(NA_real_, c(dimensions[1], nrow(X), dimensions[3]))
  sigma <- array(NA_real_, c(dimensions[1], dimensions[1], dimensions[3]))
  for (i in seq_len(dimensions[3])) {
    moments <- posterior_moments_iw(eta[, , i], X, Theta, Gamma, Xi, v)
    lambda[, , i] <- moments$Lambda
    sigma[, , i] <- moments$Sigma
  }
  list(Lambda = lambda, Sigma = sigma)
}


test_that("low-level pibble keeps matrix-t df for P = 1", {
  z <- iw_problem(P = 1, N = 3)
  nu <- z$v - z$P + 1
  eta2 <- z$eta + matrix(c(0.04, -0.03, 0.02), z$P, z$N)

  observed <- lowlevel_log_kernel(eta2, z, nu) -
    lowlevel_log_kernel(z$eta, z, nu)
  expected <- log_collapsed_niw(eta2, z) - log_collapsed_niw(z$eta, z)

  expect_equal(observed, expected, tolerance = 1e-11)
})


test_that("low-level matrix-t df and the Sylvester branch agree with NIW", {
  z <- iw_problem(P = 4, N = 2)
  nu <- z$v - z$P + 1
  eta2 <- z$eta + matrix(seq(-0.035, 0.04, length.out = z$P * z$N),
                         z$P, z$N)
  AInv <- solve(diag(z$N) + crossprod(z$X, z$Gamma %*% z$X))
  E <- eta2 - z$Theta %*% z$X
  S <- diag(z$P) + solve(z$Xi) %*% E %*% AInv %*% t(E)
  logits_kernel <- sum(z$Y[seq_len(z$P), ] * eta2) -
    sum(colSums(z$Y) * log(1 + colSums(exp(eta2))))
  manual_kernel <- logits_kernel -
    (nu + z$N + z$P - 1) * logdet_iw(S) / 2

  # This pins the low-level API to its longstanding matrix-t convention.
  expect_equal(lowlevel_log_kernel(eta2, z, nu, sylv = TRUE),
               manual_kernel, tolerance = 1e-11)

  observed <- lowlevel_log_kernel(eta2, z, nu, sylv = TRUE) -
    lowlevel_log_kernel(z$eta, z, nu, sylv = TRUE)
  expected <- log_collapsed_niw(eta2, z) - log_collapsed_niw(z$eta, z)
  expect_equal(observed, expected, tolerance = 1e-10)
})


test_that("collapsed derivatives match the normalized NIW joint for P = 1 and P > N", {
  skip_if_not_installed("numDeriv")
  for (dimensions in list(c(P = 1, N = 3), c(P = 3, N = 2))) {
    z <- iw_problem(P = dimensions[["P"]], N = dimensions[["N"]])
    nu <- z$v - z$P + 1
    theta_x <- z$Theta %*% z$X
    AInv <- solve(diag(z$N) + crossprod(z$X, z$Gamma %*% z$X))
    target <- function(x) log_collapsed_niw(x, z)

    expected_gradient <- numDeriv::grad(target, c(z$eta))
    expected_hessian <- numDeriv::hessian(target, c(z$eta))
    observed_gradient <- gradPibbleCollapsed(
      z$Y, nu, theta_x, solve(z$Xi), AInv, z$eta, sylv = TRUE
    )
    observed_hessian <- hessPibbleCollapsed(
      z$Y, nu, theta_x, solve(z$Xi), AInv, z$eta, sylv = TRUE
    )

    expect_equal(as.numeric(observed_gradient), expected_gradient,
                 tolerance = 2e-6, info = paste(dimensions, collapse = " x "))
    expect_equal(observed_hessian, expected_hessian, tolerance = 2e-5,
                 info = paste(dimensions, collapse = " x "))
  }
})


test_that("optimizer normalizes Laplace likelihood for P = 1 and P > N", {
  skip_if_not_installed("numDeriv")
  for (dimensions in list(c(P = 1, N = 3), c(P = 3, N = 2))) {
    z <- iw_problem(P = dimensions[["P"]], N = dimensions[["N"]])
    nu <- z$v - z$P + 1
    AInv <- solve(diag(z$N) + crossprod(z$X, z$Gamma %*% z$X))
    fit <- optimPibbleCollapsed(
      z$Y, nu, z$Theta %*% z$X, solve(z$Xi), AInv, z$eta,
      n_samples = 1, calcGradHess = TRUE, max_iter = 2000,
      eps_f = 1e-12, eps_g = 1e-8, useSylv = TRUE,
      ncores = 1, seed = 801
    )
    dimension <- z$P * z$N
    numeric_hessian <- numDeriv::hessian(
      function(x) log_collapsed_niw(x, z), c(fit$Pars)
    )
    independent_log_inverse_det <- -logdet_iw(-numeric_hessian)
    independent_laplace <- log_collapsed_niw(fit$Pars, z) +
      dimension * log(2 * pi) / 2 + independent_log_inverse_det / 2

    expect_equal(fit$logInvNegHessDet, independent_log_inverse_det,
                 tolerance = 2e-5,
                 info = paste(dimensions, collapse = " x "))
    expect_equal(fit$logMarginalLikelihood, independent_laplace,
                 tolerance = 2e-5,
                 info = paste(dimensions, collapse = " x "))
  }
})


test_that("uncollapse preserves fractional inverse-Wishart df for P = 1 and P > 1", {
  iterations <- 12000L
  for (P in c(1L, 2L)) {
    z <- iw_problem(P = P, N = 3)
    moments <- posterior_moments_iw(
      z$eta, z$X, z$Theta, z$Gamma, z$Xi, z$v
    )
    eta_one <- array(z$eta, dim = c(z$P, z$N, 1))
    means <- uncollapsePibble(
      eta_one, z$X, z$Theta, z$Gamma, z$Xi, z$v,
      seed = 911, ret_mean = TRUE, ncores = 1
    )

    expect_equal(matrix(means$Lambda[, , 1], z$P, nrow(z$X)),
                 moments$Lambda, tolerance = 1e-12, info = paste("P =", P))
    expect_equal(matrix(means$Sigma[, , 1], z$P, z$P),
                 moments$Sigma, tolerance = 1e-12, info = paste("P =", P))

    eta_many <- array(rep(z$eta, iterations),
                      dim = c(z$P, z$N, iterations))
    draws <- uncollapsePibble(
      eta_many, z$X, z$Theta, z$Gamma, z$Xi, z$v,
      seed = 912, ret_mean = FALSE, ncores = 1
    )
    lambda_draws <- matrix(draws$Lambda, nrow = z$P * nrow(z$X))
    sigma_mean <- apply(draws$Sigma, c(1, 2), mean)
    precision_mean <- Reduce(`+`, lapply(seq_len(iterations), function(i) {
      solve(matrix(draws$Sigma[, , i], z$P, z$P))
    })) / iterations
    lambda_covariance <- stats::cov(t(lambda_draws))
    expected_lambda_covariance <- kronecker(moments$Gamma, moments$Sigma)

    expect_equal(matrix(rowMeans(lambda_draws), z$P, nrow(z$X)),
                 moments$Lambda, tolerance = 0.04, info = paste("P =", P))
    expect_equal(matrix(sigma_mean, z$P, z$P), moments$Sigma,
                 tolerance = 0.05, info = paste("P =", P))
    expect_equal(precision_mean, moments$precision,
                 tolerance = 0.03, info = paste("P =", P))
    expect_equal(lambda_covariance, expected_lambda_covariance,
                 tolerance = 0.08, info = paste("P =", P))
  }
})


test_that("serial conjugate sampler preserves fractional inverse-Wishart df", {
  z <- iw_problem(P = 2, N = 3)
  moments <- posterior_moments_iw(
    z$eta, z$X, z$Theta, z$Gamma, z$Xi, z$v
  )
  iterations <- 12000L
  set.seed(4401)
  draws <- conjugateLinearModel(
    z$eta, z$X, z$Theta, z$Gamma, z$Xi, z$v,
    n_samples = iterations
  )
  precision_mean <- Reduce(`+`, lapply(seq_len(iterations), function(i) {
    solve(draws$Sigma[, , i])
  })) / iterations

  expect_equal(precision_mean, moments$precision, tolerance = 0.03)
})


test_that("pibble converts public IW df and returns conditional posterior means", {
  z <- iw_problem(P = 2, N = 3)
  nu <- z$v - z$P + 1
  common <- list(
    n_samples = 20, calcGradHess = TRUE, ret_mean = TRUE,
    max_iter = 2000, eps_f = 1e-12, eps_g = 1e-8, ncores = 1,
    seed = 1401
  )
  fit <- do.call(pibble, c(list(
    Y = z$Y, X = z$X, upsilon = z$v, Theta = z$Theta,
    Gamma = z$Gamma, Xi = z$Xi, init = z$eta, use_names = FALSE
  ), common))
  low <- do.call(optimPibbleCollapsed, c(list(
    Y = z$Y, upsilon = nu, ThetaX = z$Theta %*% z$X,
    KInv = solve(z$Xi),
    AInv = solve(diag(z$N) + crossprod(z$X, z$Gamma %*% z$X)),
    init = z$eta
  ), common[c("n_samples", "calcGradHess", "max_iter", "eps_f", "eps_g",
              "ncores", "seed")]))
  moments <- posterior_mean_arrays_iw(
    low$Samples, z$X, z$Theta, z$Gamma, z$Xi, z$v
  )

  expect_equal(fit$Eta, low$Samples, tolerance = 1e-10)
  expect_equal(fit$logLik, low$LogLik, tolerance = 1e-10)
  expect_equal(fit$logMarginalLikelihood, low$logMarginalLikelihood,
               tolerance = 1e-10)
  expect_equal(fit$Lambda, moments$Lambda, tolerance = 1e-10)
  expect_equal(fit$Sigma, moments$Sigma, tolerance = 1e-10)
  expect_equal(fit$upsilon, z$v)
})


test_that("single-GP basset uses public IW df and conditional means", {
  z <- iw_problem(P = 2, N = 3)
  nu <- z$v - z$P + 1
  x_gp <- matrix(seq(-0.8, 0.9, length.out = z$N), nrow = 1)
  theta_gp <- function(X) {
    rbind(0.1 + 0.04 * X[1, ], -0.08 + 0.03 * X[1, ])
  }
  gamma_gp <- function(X) {
    distance <- outer(X[1, ], X[1, ], `-`)
    0.55 * exp(-distance^2 / 0.7) + diag(0.12, ncol(X))
  }
  theta_train <- theta_gp(x_gp)
  gamma_train <- gamma_gp(x_gp)
  fit <- basset(
    z$Y, x_gp, upsilon = z$v, Theta = theta_gp, Gamma = gamma_gp,
    Xi = z$Xi, linear = 1, init = z$eta, n_samples = 20,
    calcGradHess = TRUE, ret_mean = TRUE, max_iter = 2000,
    eps_f = 1e-12, eps_g = 1e-8, ncores = 1, seed = 1402,
    use_names = FALSE
  )
  low <- optimPibbleCollapsed(
    z$Y, nu, theta_train, solve(z$Xi),
    solve(diag(z$N) + gamma_train), z$eta,
    n_samples = 20, calcGradHess = TRUE, max_iter = 2000,
    eps_f = 1e-12, eps_g = 1e-8, ncores = 1, seed = 1402
  )
  moments <- posterior_mean_arrays_iw(
    low$Samples, diag(z$N), theta_train, gamma_train, z$Xi, z$v
  )

  expect_equal(fit$Eta, low$Samples, tolerance = 1e-10)
  expect_equal(fit$logMarginalLikelihood, low$logMarginalLikelihood,
               tolerance = 1e-10)
  expect_equal(fit$Lambda, moments$Lambda, tolerance = 1e-10)
  expect_equal(fit$Sigma, moments$Sigma, tolerance = 1e-10)
  expect_equal(fit$upsilon, z$v)
})


test_that("additive basset collapses with public IW df", {
  z <- iw_problem(P = 2, N = 3)
  x_gp <- matrix(seq(-0.8, 0.9, length.out = z$N), nrow = 1)
  theta_one <- function(X) {
    rbind(0.07 + 0.02 * X[1, ], -0.03 + 0.01 * X[1, ])
  }
  theta_two <- function(X) {
    rbind(-0.02 + 0.01 * X[1, ]^2, 0.04 - 0.015 * X[1, ]^2)
  }
  gamma_one <- function(X) {
    distance <- outer(X[1, ], X[1, ], `-`)
    0.35 * exp(-distance^2 / 0.8) + diag(0.1, ncol(X))
  }
  gamma_two <- function(X) {
    distance <- outer(X[1, ], X[1, ], `-`)
    0.2 * exp(-abs(distance) / 0.9) + diag(0.08, ncol(X))
  }
  theta_sum <- theta_one(x_gp) + theta_two(x_gp)
  gamma_sum <- gamma_one(x_gp) + gamma_two(x_gp)
  common <- list(
    upsilon = z$v, Xi = z$Xi, init = z$eta,
    n_samples = 10,
    calcGradHess = TRUE, ret_mean = TRUE, max_iter = 2000,
    eps_f = 1e-12, eps_g = 1e-8, ncores = 1, seed = 1403,
    use_names = FALSE
  )
  fit <- do.call(basset, c(list(
    Y = z$Y, X = x_gp,
    Theta = list(theta_one, theta_two),
    Gamma = list(gamma_one, gamma_two), linear = 1
  ), common))
  expected <- do.call(pibble, c(list(
    Y = z$Y, X = diag(z$N), Theta = theta_sum, Gamma = gamma_sum
  ), common))

  expect_equal(fit$Eta, expected$Eta, tolerance = 1e-10)
  expect_equal(fit$Sigma, expected$Sigma, tolerance = 1e-10)
  expect_equal(fit$logMarginalLikelihood, expected$logMarginalLikelihood,
               tolerance = 1e-10)
  expect_equal(fit$upsilon, z$v)
})
