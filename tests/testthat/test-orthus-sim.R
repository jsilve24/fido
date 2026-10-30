test_that("orthus true priors have the simulated covariance as their mean", {
  for (arguments in list(list(), list(D = 4, P = 3), list(D = 4, P = 12))) {
    set.seed(42)
    sim <- do.call(orthus_sim, c(arguments, list(true_priors = TRUE)))
    dimension <- nrow(sim$Sigma)

    expect_gt(sim$upsilon, dimension + 1)
    expect_true(all(is.finite(chol(sim$Xi))))
    expect_equal(sim$Xi / (sim$upsilon - dimension - 1), sim$Sigma)
    expect_equal(sim$KInv %*% sim$Xi, diag(dimension), tolerance = 1e-12)
  }
})

test_that("orthus simulations preserve singleton matrix dimensions", {
  cases <- list(
    c(D = 2L, P = 2L, N = 4L, Q = 2L),
    c(D = 3L, P = 1L, N = 4L, Q = 2L),
    c(D = 3L, P = 2L, N = 1L, Q = 2L),
    c(D = 3L, P = 2L, N = 4L, Q = 1L),
    c(D = 2L, P = 1L, N = 1L, Q = 1L)
  )
  for (dimensions in cases) {
    for (true_priors in c(FALSE, TRUE)) {
      set.seed(42)
      sim <- do.call(orthus_sim, c(as.list(dimensions), list(true_priors = true_priors)))
      D <- dimensions[["D"]]
      P <- dimensions[["P"]]
      N <- dimensions[["N"]]
      Q <- dimensions[["Q"]]

      expect_identical(dim(sim$Y), c(D, N))
      expect_identical(dim(sim$Z), c(P, N))
      expect_identical(dim(sim$Eta), c(D - 1L, N))
      expect_identical(dim(sim$X), c(Q, N))
      expect_identical(dim(sim$Gamma), c(Q, Q))
      expect_identical(dim(sim$Sigma), rep(D - 1L + P, 2))
      expect_true(all(is.finite(chol(sim$Sigma))))
      expect_true(all(is.finite(sim$Y)) && all(sim$Y >= 0))
      expect_true(all(is.finite(sim$Z)))
      if (true_priors) {
        expect_equal(sim$Xi / (sim$upsilon - D - P), sim$Sigma)
      }
    }
  }
})

test_that("orthus simulation rejects invalid dimensions before drawing data", {
  for (name in c("D", "P", "N", "Q")) {
    for (value in list(0, -1, 1.5, NA_real_, Inf, numeric(), c(2, 3))) {
      expect_error(do.call(orthus_sim, setNames(list(value), name)),
                   paste0(name, " must be a single integer"))
    }
  }
  expect_error(orthus_sim(D = 1), "D must be a single integer >= 2")
})

test_that("orthus simulated covariances remain positive definite at unit variances", {
  set.seed(34)
  sim <- orthus_sim(D = 2, P = 2)
  expect_equal(diag(sim$Sigma)[2:3], c(1, 1))
  expect_true(all(is.finite(chol(sim$Sigma))))
})
