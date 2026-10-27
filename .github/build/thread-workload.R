args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop("Usage: thread-workload.R FIDO_LIBRARY READY START CALLS_STARTED FINISHED RESULT CORES")
}

fido_library <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
ready_file <- args[[2L]]
start_file <- args[[3L]]
calls_started_file <- args[[4L]]
finished_file <- args[[5L]]
result_file <- args[[6L]]
cores <- as.integer(args[[7L]])
if (!is.finite(cores) || cores < 1L) stop("cores must be a positive integer")

package_directory <- fido_library
if (basename(package_directory) != "fido" ||
    !file.exists(file.path(package_directory, "DESCRIPTION"))) {
  if (basename(package_directory) == "fido.so" && file.exists(package_directory)) {
    package_directory <- dirname(dirname(package_directory))
  } else if (file.exists(file.path(package_directory, "fido", "DESCRIPTION"))) {
    package_directory <- file.path(package_directory, "fido")
  } else if (file.exists(file.path(package_directory, "fido.so"))) {
    package_directory <- dirname(package_directory)
  } else if (file.exists(file.path(package_directory, "libs", "fido.so"))) {
    package_directory <- fido_library
  } else {
    stop("FIDO_LIBRARY must be an R library root, fido package directory, or fido.so path")
  }
}
library_root <- dirname(package_directory)
if (!file.exists(file.path(package_directory, "DESCRIPTION"))) {
  stop("Could not resolve an installed fido package")
}
library(fido, lib.loc = library_root)

touch <- function(path, value) {
  writeLines(value, path, useBytes = TRUE)
}

wait_for <- function(path) {
  while (!file.exists(path)) Sys.sleep(0.01)
}

scalar_text <- function(value) {
  if (length(value) == 0L || is.na(value[[1L]])) return("NA")
  if (is.logical(value[[1L]])) return(if (isTRUE(value[[1L]])) "TRUE" else "FALSE")
  if (is.numeric(value[[1L]])) return(format(value[[1L]], digits = 17, scientific = TRUE, trim = TRUE))
  as.character(value[[1L]])
}

build_info_lines <- function(prefix, info) {
  fields <- c(
    "openmp_compiled", "fido_use_parallel", "fido_use_mkl",
    "openmp_max_threads", "eigen_max_threads", "mkl_max_threads", "r_blas"
  )
  paste0(prefix, ".", fields, "=", vapply(fields, function(field) scalar_text(info[[field]]), character(1L)))
}

make_fixture <- function(dminus1, sample_count, covariate_count, draws) {
  sample_index <- seq_len(sample_count)
  design_matrix <- matrix(0, nrow = covariate_count, ncol = sample_count)
  design_matrix[1L, ] <- 1
  for (row in 2L:covariate_count) {
    design_matrix[row, ] <- cos((row - 1L) * sample_index * pi / (sample_count + 1))
  }
  if (qr(design_matrix)$rank != covariate_count) stop("fixture design matrix is not full rank")
  theta <- outer(seq_len(dminus1), seq_len(covariate_count), function(row, column) {
    0.01 * sin(row * column / 17)
  })
  gamma <- diag(seq(1, 2, length.out = covariate_count))
  xi <- diag(seq(1, 2, length.out = dminus1))
  eta_values <- seq_len(dminus1 * sample_count * draws)
  eta <- array(sin(eta_values * 0.0017) + 0.1 * cos(eta_values * 0.00031),
               dim = c(dminus1, sample_count, draws))
  list(
    eta = eta,
    design_matrix = design_matrix,
    theta = theta,
    gamma = gamma,
    xi = xi,
    upsilon = dminus1 + 12
  )
}

dminus1 <- 128L
sample_count <- 256L
covariate_count <- 8L
draws <- max(32L, cores * 4L)
fixture <- make_fixture(dminus1, sample_count, covariate_count, draws)
before <- fido_build_info()
touch(ready_file, "ready\n")
wait_for(start_file)

touch(calls_started_file, "calls-started\n")
parallel_output <- NULL
repeat_count <- 0L
measured_start <- proc.time()[["elapsed"]]
repeat {
  parallel_output <- fido:::uncollapsePibble(
    fixture$eta, fixture$design_matrix, fixture$theta, fixture$gamma, fixture$xi,
    fixture$upsilon, seed = 2849L, ret_mean = FALSE, ncores = cores
  )
  repeat_count <- repeat_count + 1L
  if (proc.time()[["elapsed"]] - measured_start >= 2) break
}
measured_elapsed <- proc.time()[["elapsed"]] - measured_start
touch(finished_file, "finished\n")
after_parallel <- fido_build_info()

serial_output <- fido:::uncollapsePibble(
  fixture$eta, fixture$design_matrix, fixture$theta, fixture$gamma, fixture$xi,
  fixture$upsilon, seed = 2849L, ret_mean = FALSE, ncores = 1L
)
after <- fido_build_info()

parallel_values <- c(as.numeric(parallel_output$Lambda), as.numeric(parallel_output$Sigma))
serial_values <- c(as.numeric(serial_output$Lambda), as.numeric(serial_output$Sigma))
finite_parallel_count <- sum(is.finite(parallel_values))
finite_serial_count <- sum(is.finite(serial_values))
finite_parallel <- finite_parallel_count == length(parallel_values)
finite_serial <- finite_serial_count == length(serial_values)
differences <- abs(parallel_values - serial_values)
relative_differences <- differences / pmax(abs(serial_values), 1)
max_absolute_difference <- max(differences)
max_relative_difference <- max(relative_differences)
absolute_tolerance <- 1e-10
relative_tolerance <- 1e-10
comparison_passed <- length(parallel_values) == length(serial_values) &&
  all(is.finite(differences)) &&
  all(differences <= absolute_tolerance + relative_tolerance * abs(serial_values))

same_runtime_value <- function(first, second) {
  if (all(is.na(c(first, second)))) return(TRUE)
  isTRUE(all.equal(as.numeric(first), as.numeric(second), tolerance = 0))
}

mkl_restored <- same_runtime_value(before$mkl_max_threads, after$mkl_max_threads) &&
  same_runtime_value(before$mkl_max_threads, after_parallel$mkl_max_threads)
actual_ncores_passed <- identical(as.integer(parallel_output$NoCores), cores) &&
  identical(as.integer(serial_output$NoCores), 1L)
build_info_passed <- all(vapply(
  list(before, after_parallel, after),
  function(info) isTRUE(info$openmp_compiled) &&
    isTRUE(info$fido_use_parallel) && isTRUE(info$fido_use_mkl),
  logical(1L)
))
measured_seconds_per_call <- measured_elapsed / repeat_count

result_lines <- c(
  build_info_lines("before", before),
  build_info_lines("after_parallel", after_parallel),
  build_info_lines("after", after),
  paste0("requested_ncores=", cores),
  paste0("parallel_actual_ncores=", scalar_text(parallel_output$NoCores)),
  paste0("serial_actual_ncores=", scalar_text(serial_output$NoCores)),
  paste0("draws=", draws),
  paste0("repeat_count=", repeat_count),
  paste0("measured_elapsed_seconds=", scalar_text(measured_elapsed)),
  paste0("measured_seconds_per_call=", scalar_text(measured_seconds_per_call)),
  paste0("finite_parallel=", scalar_text(finite_parallel)),
  paste0("finite_parallel_count=", finite_parallel_count),
  paste0("finite_serial=", scalar_text(finite_serial)),
  paste0("finite_serial_count=", finite_serial_count),
  paste0("max_absolute_difference=", scalar_text(max_absolute_difference)),
  paste0("max_relative_difference=", scalar_text(max_relative_difference)),
  paste0("absolute_tolerance=", scalar_text(absolute_tolerance)),
  paste0("relative_tolerance=", scalar_text(relative_tolerance)),
  paste0("comparison_passed=", scalar_text(comparison_passed)),
  paste0("mkl_restored=", scalar_text(mkl_restored)),
  paste0("actual_ncores_passed=", scalar_text(actual_ncores_passed)),
  paste0("build_info_passed=", scalar_text(build_info_passed))
)
writeLines(result_lines, result_file, useBytes = TRUE)
if (!finite_parallel || !finite_serial || !comparison_passed || !mkl_restored ||
    !actual_ncores_passed || !build_info_passed) {
  stop("workload correctness or runtime restoration check failed")
}
