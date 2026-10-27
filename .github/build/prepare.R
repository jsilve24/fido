args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || !args[[1L]] %in% c("default", "mkl", "mac-stock", "mac-libomp")) {
  stop("Usage: prepare.R PROFILE MAKEVARS [MKL_OR_LIBOMP_PREFIX]")
}
profile <- args[[1L]]
makevars <- args[[2L]]
settings <- character()
if (profile %in% c("mkl", "mac-libomp")) {
  if (length(args) != 3L) stop("This profile requires a library prefix")
  prefix <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
  if (grepl("[[:space:]#$]", prefix)) stop("Unsupported characters in CI prefix")
}
if (profile == "mkl") {
  if (Sys.info()[["sysname"]] != "Linux" || R.version$arch != "x86_64") {
    stop("The MKL CI profile requires Linux x86_64")
  }
  libraries <- file.path(prefix, "lib", c(
    "libmkl_intel_lp64.so.2", "libmkl_gnu_thread.so.2", "libmkl_core.so.2"
  ))
  stopifnot(all(file.exists(c(file.path(prefix, "include", "mkl.h"), libraries))))
  settings <- c(
    paste("PKG_CPPFLAGS += -DFIDO_USE_MKL -DEIGEN_USE_MKL_ALL", paste0("-I", prefix, "/include")),
    paste("PKG_LIBS += -Wl,--no-as-needed", paste(libraries, collapse = " "),
          paste0("-Wl,-rpath,", prefix, "/lib"), "-lgomp -lpthread -lm -ldl")
  )
}
if (profile %in% c("mac-stock", "mac-libomp")) {
  if (Sys.info()[["sysname"]] != "Darwin") stop("This profile requires macOS")
  settings <- paste(c("CXX", "CXX11", "CXX14", "CXX17", "CXX20", "CXX23"), "= /usr/bin/clang++")
  openmp <- ""
  if (profile == "mac-libomp") {
    stopifnot(all(file.exists(file.path(prefix, c("include/omp.h", "lib/libomp.dylib")))))
    openmp <- paste("-Xclang -fopenmp", paste0("-I", prefix, "/include"))
    settings <- c(settings, paste("PKG_LIBS +=", paste0("-L", prefix, "/lib"),
                                 paste0("-Wl,-rpath,", prefix, "/lib"), "-lomp"))
  }
  settings <- c(settings, paste("SHLIB_OPENMP_CXXFLAGS =", openmp))
}
writeLines(settings, makevars)
cat("CI build profile:", profile, "\n")
cat("R_MAKEVARS_USER:", normalizePath(makevars, winslash = "/"), "\n")
cat(settings, sep = "\n")
