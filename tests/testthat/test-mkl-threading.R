test_that("sampling and Hessian regions preserve MKL limits and results", {
  design <- rbind(1, c(-1, 0, 1, 2))
  prior_mean <- matrix(0, 2, 2)
  prior_covariance <- diag(2)
  draws <- array(seq_len(32) / 100, c(2, 4, 4))
  known_covariance <- array(rep(diag(2), 4), c(2, 2, 4))
  counts <- matrix(seq_len(12), 3, 4)
  column_precision <- solve(diag(4) + crossprod(design))
  original_limit <- fido_build_info()$mkl_max_threads
  results <- list()

  for (cores in c(1L, 2L)) {
    sampled <- uncollapsePibble(draws, design, prior_mean, prior_covariance,
                                prior_covariance, 7, seed = 4201,
                                ret_mean = FALSE, ncores = cores)
    expect_identical(fido_build_info()$mkl_max_threads, original_limit)
    known <- uncollapsePibble_sigmaKnown(
      draws, design, prior_mean, prior_covariance, diag(ncol(design)),
      prior_covariance, known_covariance, 7, seed = 4201,
      ret_mean = FALSE, ncores = cores
    )
    expect_identical(fido_build_info()$mkl_max_threads, original_limit)
    pibble_hessian <- hessPibbleCollapsed(
      counts, 6, prior_mean %*% design, prior_covariance,
      column_precision, draws[, , 1]
    )
    expect_identical(fido_build_info()$mkl_max_threads, original_limit)
    maltipoo_hessian <- hessMaltipooCollapsed(
      counts, 6, prior_mean, design, prior_covariance,
      prior_covariance, draws[, , 1], 0
    )
    expect_identical(fido_build_info()$mkl_max_threads, original_limit)
    results[[cores]] <- list(sampled$Lambda, sampled$Sigma, known$Lambda,
                             pibble_hessian, maltipoo_hessian)
    expect_true(all(is.finite(unlist(results[[cores]]))))
  }
  expect_equal(results[[1L]], results[[2L]], tolerance = 1e-12)
})
