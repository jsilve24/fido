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


test_that("fresh processes limit check threads without changing user defaults", {
  script <- tempfile(fileext = ".R")
  result <- tempfile(fileext = ".rds")
  on.exit(unlink(c(script, result)), add = TRUE)
  writeLines(c(
    paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
    "saveRDS(fido::fido_build_info(), commandArgs(trailingOnly = TRUE)[1])"
  ), script)

  run_info <- function(check_limit, threads) {
    variables <- c(
      `_R_CHECK_LIMIT_CORES_` = check_limit,
      OMP_NUM_THREADS = as.character(threads),
      MKL_NUM_THREADS = as.character(threads)
    )
    previous <- Sys.getenv(names(variables), unset = NA_character_)
    on.exit({
      Sys.unsetenv(names(previous)[is.na(previous)])
      do.call(Sys.setenv, as.list(previous[!is.na(previous)]))
    }, add = TRUE)
    do.call(Sys.setenv, as.list(variables))
    unlink(result)
    output <- system2(
      file.path(R.home("bin"), "Rscript"),
      c("--vanilla", shQuote(script), shQuote(result)),
      stdout = TRUE, stderr = TRUE
    )
    expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
    expect_true(file.exists(result))
    readRDS(result)
  }

  regular <- run_info("FALSE", 4)
  limited <- run_info("TRUE", 4)
  lower <- run_info("true", 1)

  for (field in c("openmp_max_threads", "eigen_max_threads", "mkl_max_threads")) {
    if (!is.na(regular[[field]])) {
      expect_equal(limited[[field]], min(regular[[field]], 2L), info = field)
      expect_equal(lower[[field]], 1L, info = field)
    }
  }
  if (regular$openmp_compiled) {
    expect_equal(regular$openmp_max_threads, 4L)
    expect_equal(regular$eigen_max_threads, 4L)
  }
})
