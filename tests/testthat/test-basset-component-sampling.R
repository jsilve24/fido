context("test-basset-component-sampling")


repeat_matrix <- function(x, iterations) {
  array(rep(x, iterations), dim = c(dim(x), iterations))
}


conditional_component_moments <- function(total, component_mean,
                                          component_covariance, design,
                                          other_mean, other_covariance,
                                          sigma) {
  centered_total <- total - other_mean
  total_covariance <- t(design) %*% component_covariance %*% design +
    other_covariance
  posterior_mean <- component_mean +
    (centered_total - component_mean %*% design) %*%
      solve(total_covariance, t(design) %*% component_covariance)
  posterior_covariance <- component_covariance -
    component_covariance %*% design %*%
      solve(total_covariance, t(design) %*% component_covariance)

  list(
    mean = posterior_mean,
    covariance = kronecker(posterior_covariance, sigma)
  )
}


test_that("additive basset samples GP-first and linear-first components correctly", {
  P <- 2L
  N <- 3L
  Q <- 2L
  iterations <- 5000L

  Y <- matrix(c(8, 5, 3, 4, 7, 6, 2, 3, 9), nrow = P + 1L)
  X <- rbind(1, c(-0.8, 0.2, 1.1))
  theta_linear <- matrix(c(0.12, -0.08, 0.18, 0.09), nrow = P)
  gamma_linear <- matrix(c(0.7, 0.12, 0.12, 0.5), nrow = Q)
  theta_gp <- function(x) {
    rbind(0.05 + 0.08 * x[2, ], -0.04 + 0.03 * x[2, ])
  }
  gamma_gp <- function(x) {
    distance <- outer(x[2, ], x[2, ], `-`)
    0.45 * exp(-distance^2 / 0.9) + diag(0.15, ncol(x))
  }
  theta_gp_extra <- function(x) {
    rbind(-0.03 + 0.02 * x[2, ], 0.06 - 0.04 * x[2, ])
  }
  gamma_gp_extra <- function(x) {
    distance <- outer(x[2, ], x[2, ], `-`)
    0.25 * exp(-abs(distance) / 0.7) + diag(0.1, ncol(x))
  }
  combined_lambda <- matrix(
    c(0.55, -0.30, -0.20, 0.25, 0.35, 0.10), nrow = P
  )
  sigma <- matrix(c(0.8, 0.22, 0.22, 0.55), nrow = P)
  combined_draws <- repeat_matrix(combined_lambda, iterations)
  sigma_draws <- repeat_matrix(sigma, iterations)

  collapsed_fit <- list(
    Eta = array(0, dim = c(P, N, iterations)),
    Lambda = combined_draws,
    Sigma = sigma_draws,
    N = N,
    Q = N,
    D = P + 1L,
    init = matrix(0, P, N),
    logMarginalLikelihood = 0
  )
  basset_with_fixed_joint_draws <- basset
  environment(basset_with_fixed_joint_draws) <- list2env(
    list(pibble = function(...) collapsed_fit),
    parent = environment(basset)
  )

  gp_mean <- theta_gp(X)
  gp_covariance <- gamma_gp(X)
  extra_gp_mean <- theta_gp_extra(X)
  extra_gp_covariance <- gamma_gp_extra(X)
  linear_mean <- theta_linear %*% X
  linear_covariance <- t(X) %*% gamma_linear %*% X
  branches <- list(
    `GP first` = list(
      Theta = list(theta_gp, theta_linear, theta_gp_extra),
      Gamma = list(gamma_gp, gamma_linear, gamma_gp_extra),
      component_mean = gp_mean,
      component_covariance = gp_covariance,
      design = diag(N),
      other_mean = linear_mean + extra_gp_mean,
      other_covariance = linear_covariance + extra_gp_covariance,
      to_sample_scale = list(
        function(lambda) lambda,
        function(lambda) {
          out <- array(NA_real_, dim = c(P, N, dim(lambda)[3]))
          for (i in seq_len(dim(lambda)[3])) {
            out[, , i] <- lambda[, , i] %*% X
          }
          out
        },
        function(lambda) lambda
      )
    ),
    `linear first` = list(
      Theta = list(theta_linear, theta_gp, theta_gp_extra),
      Gamma = list(gamma_linear, gamma_gp, gamma_gp_extra),
      component_mean = theta_linear,
      component_covariance = gamma_linear,
      design = X,
      other_mean = gp_mean + extra_gp_mean,
      other_covariance = gp_covariance + extra_gp_covariance,
      to_sample_scale = list(
        function(lambda) {
          out <- array(NA_real_, dim = c(P, N, dim(lambda)[3]))
          for (i in seq_len(dim(lambda)[3])) {
            out[, , i] <- lambda[, , i] %*% X
          }
          out
        },
        function(lambda) lambda,
        function(lambda) lambda
      )
    )
  )

  fit_component <- function(branch, ret_mean, ncores, seed) {
    basset_with_fixed_joint_draws(
      Y = Y, X = X, upsilon = 7, Theta = branch$Theta,
      Gamma = branch$Gamma, Xi = diag(P), linear = seq_len(Q),
      n_samples = iterations, ret_mean = ret_mean,
      ncores = ncores, seed = seed
    )
  }

  for (branch_name in names(branches)) {
    branch <- branches[[branch_name]]
    expected <- conditional_component_moments(
      combined_lambda, branch$component_mean,
      branch$component_covariance, branch$design,
      branch$other_mean, branch$other_covariance, sigma
    )

    mean_fit <- fit_component(branch, ret_mean = TRUE, ncores = 0, seed = 9101)
    expect_equal(
      mean_fit$Lambda[[1]][, , 1], expected$mean,
      tolerance = 1e-10, info = branch_name
    )
    expect_equal(
      mean_fit$Lambda[[1]],
      repeat_matrix(expected$mean, iterations),
      tolerance = 1e-10, info = branch_name
    )

    sampled_one_core <- fit_component(
      branch, ret_mean = FALSE, ncores = 1, seed = 9201
    )
    sampled_two_cores <- fit_component(
      branch, ret_mean = FALSE, ncores = 2, seed = 9201
    )
    sampled_other_seed <- fit_component(
      branch, ret_mean = FALSE, ncores = 1, seed = 9202
    )

    component_draws <- matrix(
      sampled_one_core$Lambda[[1]],
      nrow = length(expected$mean), ncol = iterations
    )
    mean_standard_error <- sqrt(diag(expected$covariance) / iterations)
    covariance_estimate <- cov(t(component_draws))
    covariance_standard_error <- sqrt(
      (expected$covariance^2 +
         outer(diag(expected$covariance), diag(expected$covariance))) /
        (iterations - 1)
    )
    expect_true(
      max(abs(rowMeans(component_draws) - c(expected$mean)) /
            mean_standard_error) < 5,
      info = branch_name
    )
    expect_true(
      max(abs(covariance_estimate - expected$covariance) /
            covariance_standard_error) < 5,
      info = branch_name
    )
    expect_identical(
      sampled_one_core$Lambda, sampled_two_cores$Lambda,
      info = branch_name
    )
    expect_false(
      identical(sampled_one_core$Lambda, sampled_other_seed$Lambda),
      info = branch_name
    )

    components_on_sample_scale <- Map(
      function(project, lambda) project(lambda),
      branch$to_sample_scale, sampled_one_core$Lambda
    )
    expect_equal(
      Reduce(`+`, components_on_sample_scale),
      combined_draws, tolerance = 1e-12, info = branch_name
    )
  }
})


test_that("additive GP components have the joint conditional cross-covariance", {
  P <- 2L
  N <- 3L
  iterations <- 5000L
  Y <- matrix(c(6, 4, 2, 3, 7, 5, 4, 2, 8), nrow = P + 1L)
  X <- rbind(1, c(-0.7, 0.1, 0.9))
  sigma <- matrix(c(0.75, 0.18, 0.18, 0.5), nrow = P)
  total <- matrix(c(0.4, -0.2, -0.1, 0.3, 0.25, 0.05), nrow = P)
  kernels <- list(
    matrix(c(0.7, 0.16, 0.08, 0.16, 0.6, 0.12,
             0.08, 0.12, 0.55), nrow = N),
    matrix(c(0.45, 0.09, 0.04, 0.09, 0.5, 0.11,
             0.04, 0.11, 0.4), nrow = N),
    matrix(c(0.35, 0.06, 0.03, 0.06, 0.3, 0.07,
             0.03, 0.07, 0.38), nrow = N)
  )
  zero_mean <- function(x) matrix(0, P, ncol(x))
  gamma_functions <- lapply(kernels, function(kernel) {
    force(kernel)
    function(x) kernel
  })
  combined_draws <- repeat_matrix(total, iterations)
  collapsed_fit <- list(
    Eta = array(0, dim = c(P, N, iterations)),
    Lambda = combined_draws,
    Sigma = repeat_matrix(sigma, iterations),
    N = N,
    Q = N,
    D = P + 1L,
    init = matrix(0, P, N),
    logMarginalLikelihood = 0
  )
  basset_with_fixed_joint_draws <- basset
  environment(basset_with_fixed_joint_draws) <- list2env(
    list(pibble = function(...) collapsed_fit),
    parent = environment(basset)
  )

  fit <- basset_with_fixed_joint_draws(
    Y = Y, X = X, upsilon = 7,
    Theta = rep(list(zero_mean), 3), Gamma = gamma_functions,
    Xi = diag(P), linear = seq_len(nrow(X)),
    n_samples = iterations, ret_mean = FALSE, ncores = 1, seed = 9301
  )

  kernel_sum <- Reduce(`+`, kernels)
  column_covariance_11 <- kernels[[1]] -
    kernels[[1]] %*% solve(kernel_sum, kernels[[1]])
  column_covariance_22 <- kernels[[2]] -
    kernels[[2]] %*% solve(kernel_sum, kernels[[2]])
  column_cross_covariance <-
    -kernels[[1]] %*% solve(kernel_sum, kernels[[2]])
  covariance_11 <- kronecker(column_covariance_11, sigma)
  covariance_22 <- kronecker(column_covariance_22, sigma)
  expected_cross_covariance <- kronecker(column_cross_covariance, sigma)
  draws_one <- matrix(fit$Lambda[[1]], nrow = P * N, ncol = iterations)
  draws_two <- matrix(fit$Lambda[[2]], nrow = P * N, ncol = iterations)
  estimated_cross_covariance <- cov(t(draws_one), t(draws_two))
  cross_covariance_standard_error <- sqrt(
    (expected_cross_covariance^2 +
       outer(diag(covariance_11), diag(covariance_22))) /
      (iterations - 1)
  )

  expect_true(
    max(abs(estimated_cross_covariance - expected_cross_covariance) /
          cross_covariance_standard_error) < 5
  )
  expect_equal(Reduce(`+`, fit$Lambda), combined_draws, tolerance = 1e-12)
})
