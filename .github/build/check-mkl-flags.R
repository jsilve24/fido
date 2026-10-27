args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: check-mkl-flags.R CHECK_DIR")

check_dir <- normalizePath(args[[1L]], mustWork = TRUE)
install_logs <- list.files(check_dir, "00install\\.out$", recursive = TRUE, full.names = TRUE)
if (!length(install_logs)) stop("No 00install.out files found under CHECK_DIR")

cpp_names <- list.files("src", "\\.cpp$", full.names = FALSE)
if (!length(cpp_names)) stop("No repository src/*.cpp files found")

compiler_lines <- unlist(lapply(install_logs, function(path) {
  lines <- readLines(path, warn = FALSE)
  has_cpp_name <- Reduce("|", lapply(cpp_names, function(name) {
    grepl(name, lines, fixed = TRUE)
  }))
  lines[grepl("g\\+\\+", lines) & has_cpp_name]
}), use.names = FALSE)

missing_sources <- cpp_names[!vapply(cpp_names, function(name) {
  any(grepl("g\\+\\+", compiler_lines) & grepl(name, compiler_lines, fixed = TRUE))
}, logical(1L))]
if (length(missing_sources)) {
  stop("Missing package g++ compiler lines for: ", paste(missing_sources, collapse = ", "))
}

required_flags <- c(
  "-DFIDO_USE_MKL",
  "-DEIGEN_USE_MKL_ALL",
  "-I../inst/include/",
  "-fopenmp"
)
missing_flags <- vapply(required_flags, function(flag) {
  !all(grepl(flag, compiler_lines, fixed = TRUE))
}, logical(1L))
if (any(missing_flags)) {
  stop("Missing required MKL flags on compiler lines: ",
       paste(required_flags[missing_flags], collapse = ", "))
}

cat("MKL compiler flag check passed:", length(compiler_lines),
    "package g++ compiler lines across", length(install_logs), "install log(s)\n")
