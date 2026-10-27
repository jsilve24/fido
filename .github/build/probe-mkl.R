source_file <- normalizePath(".github/build/mkl-probe.cpp", mustWork = TRUE)
stopifnot(nzchar(Sys.getenv("R_MAKEVARS_USER")))
build_dir <- tempfile("fido-mkl-probe-")
dir.create(build_dir)
file.copy(source_file, file.path(build_dir, "mkl-probe.cpp"))
writeLines(c(
  paste0("PKG_CPPFLAGS = -I", shQuote(system.file("include", package = "RcppEigen"))),
  "PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS)",
  "PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS)"
), file.path(build_dir, "Makevars"))
setwd(build_dir)
status <- system2(file.path(R.home("bin"), "R"), c("CMD", "SHLIB", "mkl-probe.cpp"))
if (status != 0L) stop("MKL/Eigen probe failed to compile or link")
library_path <- paste0("mkl-probe", .Platform$dynlib.ext)
probe <- dyn.load(library_path)
result <- .C(getNativeSymbolInfo("fido_mkl_probe", probe), status = integer(1), version = integer(3))
print(result)
stopifnot(identical(result$status, 0L), identical(result$version, c(2024L, 0L, 2L)))
dyn.unload(library_path)
cat("MKL probe passed: 2024.2.2, LP64, GNU threading, R dyn.load, dpotrf and Eigen eigensolver\n")
