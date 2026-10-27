check_build_expectation <- function(name, observed) {
  expected <- Sys.getenv(name, unset = "")
  if (nzchar(expected)) {
    expect_true(
      expected %in% c("true", "false"),
      info = paste0(name, " must be 'true' or 'false'")
    )
    if (expected %in% c("true", "false")) {
      expect_identical(observed, expected == "true")
    }
  }
}

test_that("fido_build_info has a stable schema", {
  info <- fido_build_info()
  fields <- c(
    "openmp_compiled",
    "fido_use_parallel",
    "fido_use_mkl",
    "openmp_max_threads",
    "eigen_max_threads",
    "mkl_max_threads",
    "r_blas"
  )

  expect_type(info, "list")
  expect_identical(names(info), fields)
  expect_type(info$openmp_compiled, "logical")
  expect_type(info$fido_use_parallel, "logical")
  expect_type(info$fido_use_mkl, "logical")
  expect_length(info$openmp_compiled, 1)
  expect_length(info$fido_use_parallel, 1)
  expect_length(info$fido_use_mkl, 1)
  expect_type(info$openmp_max_threads, "integer")
  expect_type(info$eigen_max_threads, "integer")
  expect_type(info$mkl_max_threads, "integer")
  expect_type(info$r_blas, "character")
  expect_length(info$r_blas, 1)
  expect_identical(info$r_blas, unname(extSoftVersion()["BLAS"]))

  if (any(nzchar(Sys.getenv(c("FIDO_EXPECT_OPENMP", "FIDO_EXPECT_MKL"), unset = "")))) {
    print(info)
  }
  check_build_expectation("FIDO_EXPECT_OPENMP", info$openmp_compiled)
  check_build_expectation("FIDO_EXPECT_OPENMP", info$fido_use_parallel)
  check_build_expectation("FIDO_EXPECT_MKL", info$fido_use_mkl)
})

test_that("fido_build_info reports native thread capabilities", {
  info <- fido_build_info()

  if (isTRUE(info$openmp_compiled)) {
    expect_gt(info$openmp_max_threads, 0L)
  } else {
    expect_true(is.na(info$openmp_max_threads))
  }
  expect_gt(info$eigen_max_threads, 0L)

  if (isTRUE(info$fido_use_mkl)) {
    expect_gt(info$mkl_max_threads, 0L)
  } else {
    expect_true(is.na(info$mkl_max_threads))
  }
})

test_that("fido_build_info does not change R's random stream", {
  set.seed(42)
  seed <- .Random.seed
  fido_build_info()
  expect_identical(.Random.seed, seed)
})
