args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: install-intel-comparison.R SOURCE MKL_PREFIX LIBRARY")
}
source <- normalizePath(args[[1L]], mustWork = TRUE)
prefix <- normalizePath(args[[2L]], mustWork = TRUE)
library_path <- args[[3L]]
dir.create(library_path, recursive = TRUE, showWarnings = FALSE)
library_path <- normalizePath(library_path, mustWork = TRUE)
libraries <- file.path(prefix, "lib", c(
  "libmkl_intel_lp64.so.2", "libmkl_intel_thread.so.2",
  "libmkl_core.so.2", "libiomp5.so"
))
stopifnot(all(file.exists(libraries)))
temporary <- tempfile("fido-intel-comparison-")
dir.create(temporary)
stopifnot(file.copy(source, temporary, recursive = TRUE, copy.date = TRUE))
comparison_source <- file.path(temporary, basename(source))
unlink(file.path(comparison_source, "src", paste0("fido", .Platform$dynlib.ext)))
makevars <- file.path(temporary, "Makevars")
writeLines(paste(
  "FIDO_MKL_LIBS = -Wl,--no-as-needed",
  paste(shQuote(libraries), collapse = " "),
  paste0("-Wl,-rpath,", shQuote(file.path(prefix, "lib"))),
  "-lpthread -lm -ldl"
), makevars)
Sys.setenv(R_MAKEVARS_USER = makevars, MKLROOT = prefix)
cat("TEST ONLY: Intel MKL threading with GNU OpenMP; not a supported install mode\n")
status <- system2(file.path(R.home("bin"), "R"), c(
  "CMD", "INSTALL", "--install-tests",
  paste0("--library=", shQuote(library_path)), shQuote(comparison_source)
))
if (status != 0L) stop("Intel comparison installation failed")
unlink(temporary, recursive = TRUE)
