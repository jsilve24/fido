context("test-laplaceapproximation")

test_that("eigen_lap gets correct answer", {
  n_samples <- 100000
  m <- 1:3
  S <- diag(4:6)
  S[1,2] <- S[2,1] <- -1
  #S <- -S
  z <- eigen_lap_test(n_samples, m, S, 0)
  
  expect_equal(var(t(z)), solve(S), tolerance=0.005)
  expect_equal(rowMeans(z), m, tolerance=.01)
})


test_that("cholesky_lap gets correct answer", {
  n_samples <- 100000
  m <- 1:3
  S <- diag(4:6)
  S[1,2] <- S[2,1] <- -1
  #S <- -S
  z <- cholesky_lap_test(n_samples, m, S, 0)
  
  expect_equal(var(t(z)), solve(S), tolerance=0.005)
  expect_equal(rowMeans(z), m, tolerance=.01)
})

test_that("LaplaceApproximation rejects a non-positive-definite Hessian", {
  samples <- NULL
  expect_warning(
    expect_error(
      samples <- LaplaceApproximation_test(10, c(0, 0), diag(c(1, -1)),
                                           "cholesky", 0),
      "decomposition failed"
    ),
    "Cholesky of Hessian failed"
  )
  expect_null(samples)
})

test_that("LaplaceApproximation gets correct result for full hessian", {
  n_samples <- 1000000
  m <- 1:3
  S <- diag(1:3)
  S[1,2] <- S[2,1] <- -1
  #S <- -S
  
  z <- LaplaceApproximation_test(n_samples, m, S, "eigen", 0)
  expect_equal(var(t(z)), solve(S), tolerance=0.005)
  expect_equal(rowMeans(z), m, tolerance=.01)
  
  z <- LaplaceApproximation_test(n_samples, m, S, "cholesky", 0)
  expect_equal(var(t(z)), solve(S), tolerance=0.005)
  expect_equal(rowMeans(z), m, tolerance=.01)
})

test_that("LaplaceApproximation handles row-stacked Hessian blocks", {
  first_block <- matrix(c(4, 1, 1, 3), 2)
  second_block <- matrix(c(6, -2, -2, 5), 2)
  block_hessian <- rbind(first_block, second_block)
  full_hessian <- matrix(0, 4, 4)
  full_hessian[1:2, 1:2] <- first_block
  full_hessian[3:4, 3:4] <- second_block
  location <- c(1, -2, 3, 4)

  for (method in c("cholesky", "eigen")) {
    if (fido_build_info()$fido_use_mkl) {
      expect_error(
        LaplaceApproximation_test(10, location, block_hessian, method, 0,
                                  seed = 2401),
        "Partial Hessian Not Implemented for MKL"
      )
    } else {
      samples <- LaplaceApproximation_test(100000, location, block_hessian,
                                           method, 0, seed = 2401)
      expect_identical(dim(samples), c(4L, 100000L))
      expect_true(all(is.finite(samples)))
      expect_equal(rowMeans(samples), location, tolerance = 0.01)
      expect_equal(var(t(samples)), solve(full_hessian), tolerance = 0.005)
      if (method == "cholesky") {
        full_samples <- LaplaceApproximation_test(100000, location, full_hessian,
                                                  method, 0, seed = 2401)
        expect_equal(samples, full_samples, tolerance = 1e-12)
      }
    }
  }
})

test_that("LaplaceApproximation rejects incomplete Hessian blocks", {
  expect_error(
    LaplaceApproximation_test(10, rep(0, 3), matrix(1, 3, 2), "cholesky", 0),
    "Rectangular Hessian of wrong dimension passed"
  )
})

test_that("LaplaceApproximation rejects a non-positive-definite later block", {
  block_hessian <- rbind(diag(2), diag(c(1, -1)))
  samples <- NULL
  if (fido_build_info()$fido_use_mkl) {
    expect_error(
      samples <- LaplaceApproximation_test(10, rep(0, 4), block_hessian,
                                           "cholesky", 0, seed = 2401),
      "Partial Hessian Not Implemented for MKL"
    )
  } else {
    expect_warning(
      expect_error(
        samples <- LaplaceApproximation_test(10, rep(0, 4), block_hessian,
                                             "cholesky", 0, seed = 2401),
        "decomposition failed"
      ),
      "Cholesky of Hessian failed"
    )
  }
  expect_null(samples)
})

