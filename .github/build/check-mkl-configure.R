args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript .github/build/check-mkl-configure.R MKL_PREFIX OUTPUT_DIR",
       call. = FALSE)
}

if (!identical(unname(Sys.info()[["sysname"]]), "Linux")) {
  stop("check-mkl-configure.R is Linux-only", call. = FALSE)
}

repo_root <- normalizePath(getwd(), mustWork = TRUE)
mkl_prefix <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
if (!dir.create(output_dir, recursive = TRUE, showWarnings = FALSE) &&
    !dir.exists(output_dir)) {
  stop("Could not create OUTPUT_DIR: ", output_dir, call. = FALSE)
}

source_files <- c(
  "configure",
  "DESCRIPTION",
  file.path("src", "Makevars.in"),
  file.path("inst", "include", "fidoGenerated.h.in")
)
missing_source_files <- source_files[!file.exists(file.path(repo_root, source_files))]
if (length(missing_source_files)) {
  stop("Missing tracked configure input(s): ",
       paste(missing_source_files, collapse = ", "), call. = FALSE)
}

find_existing_directory <- function(paths, required_file) {
  paths <- unique(paths[dir.exists(paths)])
  matches <- paths[vapply(paths, function(path) {
    file.exists(file.path(path, required_file))
  }, logical(1L))]
  if (length(matches)) matches[[1L]] else character()
}

include_candidates <- c(
  file.path(mkl_prefix, "include"),
  file.path(mkl_prefix, "include", "mkl")
)
include_dir <- find_existing_directory(include_candidates, "mkl.h")
if (is.null(include_dir) || !length(include_dir)) {
  stop("MKL_PREFIX does not contain include/mkl.h: ", mkl_prefix,
       call. = FALSE)
}

mkl_headers <- list.files(mkl_prefix, pattern = "^mkl\\.h$", recursive = TRUE,
                          full.names = TRUE)
if (!length(mkl_headers)) {
  stop("MKL_PREFIX does not contain mkl.h: ", mkl_prefix, call. = FALSE)
}

library_candidates <- c(
  file.path(mkl_prefix, "lib"),
  file.path(mkl_prefix, "lib", "intel64"),
  file.path(mkl_prefix, "lib", "intel64_lin")
)
library_candidates <- unique(library_candidates[dir.exists(library_candidates)])
library_candidates <- library_candidates[vapply(library_candidates, function(path) {
  length(list.files(path, pattern = "^libmkl_core\\.so", full.names = TRUE)) > 0L
}, logical(1L))]
if (!length(library_candidates)) {
  stop("MKL_PREFIX does not contain a usable MKL library directory: ",
       mkl_prefix, call. = FALSE)
}
mkl_lib_dir <- library_candidates[[1L]]

pkg_config_files <- list.files(mkl_prefix, pattern = "\\.pc$", recursive = TRUE,
                               full.names = TRUE)
pkg_config_files <- pkg_config_files[vapply(pkg_config_files, function(path) {
  lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
  grepl("mkl", basename(path), ignore.case = TRUE) ||
    any(grepl("mkl", lines, ignore.case = TRUE))
}, logical(1L))]
pkg_config_dir <- if (length(pkg_config_files)) dirname(pkg_config_files[[1L]]) else ""
pkg_config_candidates <- c(
  file.path(mkl_prefix, "bin", "pkg-config"),
  file.path(mkl_prefix, "bin", "pkgconfig"),
  unname(Sys.which("pkg-config"))
)
pkg_config_candidates <- pkg_config_candidates[nzchar(pkg_config_candidates)]
pkg_config_candidates <- pkg_config_candidates[file.access(pkg_config_candidates, 1L) == 0L]
if (!length(pkg_config_candidates)) {
  stop("Could not find a usable pkg-config executable", call. = FALSE)
}
pkg_config_command <- normalizePath(pkg_config_candidates[[1L]], mustWork = TRUE)

work_root <- tempfile("fido-mkl-configure-")
dir.create(work_root, recursive = TRUE)
on.exit(unlink(work_root, recursive = TRUE, force = TRUE), add = TRUE)

empty_pkg_config_dir <- file.path(work_root, "empty-pkg-config")
dir.create(empty_pkg_config_dir, recursive = TRUE)
makevars_user <- file.path(work_root, "Makevars.user")
makevars_site <- file.path(work_root, "Makevars.site")
file.create(makevars_user)
file.create(makevars_site)

make_broken_prefix <- function() {
  broken_prefix <- file.path(work_root, "broken-mkl")
  dir.create(broken_prefix, recursive = TRUE)
  include_relative <- substring(include_dir, nchar(mkl_prefix) + 2L)
  include_destination <- file.path(broken_prefix, include_relative)
  dir.create(dirname(include_destination), recursive = TRUE,
             showWarnings = FALSE)
  if (!file.symlink(include_dir, include_destination)) {
    stop("Could not symlink the real MKL include directory for probe testing",
         call. = FALSE)
  }
  library_relative <- substring(mkl_lib_dir, nchar(mkl_prefix) + 2L)
  broken_lib_dir <- file.path(broken_prefix, library_relative)
  dir.create(broken_lib_dir, recursive = TRUE)
  library_files <- list.files(mkl_lib_dir, all.files = TRUE, full.names = TRUE,
                              no.. = TRUE)
  library_files <- library_files[!vapply(library_files, dir.exists, logical(1L))]
  if (!length(library_files)) {
    stop("Could not find MKL libraries for probe testing", call. = FALSE)
  }
  broken_names <- grepl("^libmkl_(core|gnu_thread)(\\.so|\\.a)",
                        basename(library_files))
  for (index in seq_along(library_files)) {
    destination <- file.path(broken_lib_dir, basename(library_files[[index]]))
    if (broken_names[[index]]) {
      writeLines("not an ELF shared library", destination)
    } else if (!file.symlink(library_files[[index]], destination)) {
      stop("Could not symlink a real MKL library for probe testing",
           call. = FALSE)
    }
  }
  if (!any(broken_names)) {
    writeLines("not an ELF shared library",
               file.path(broken_lib_dir, "libmkl_core.so"))
  }
  broken_prefix
}

broken_prefix <- make_broken_prefix()
make_missing_thread_api_prefix <- function() {
  prefix <- file.path(work_root, "missing-thread-api-mkl")
  include_destination <- file.path(prefix, "include")
  dir.create(include_destination, recursive = TRUE)
  header_path <- file.path(include_destination, "mkl.h")
  writeLines(c(
    paste0('#include "', normalizePath(file.path(include_dir, "mkl.h"),
                                       winslash = "/", mustWork = TRUE), '"'),
    "#undef mkl_set_num_threads_local",
    "#define mkl_set_num_threads_local fido_missing_mkl_thread_api"
  ), header_path)
  include_entries <- list.files(include_dir, all.files = TRUE, full.names = TRUE,
                                no.. = TRUE)
  include_entries <- include_entries[basename(include_entries) != "mkl.h"]
  for (entry in include_entries) {
    if (!file.symlink(entry, file.path(include_destination, basename(entry)))) {
      stop("Could not symlink an MKL header for thread API testing: ", entry,
           call. = FALSE)
    }
  }
  if (!file.symlink(mkl_lib_dir, file.path(prefix, "lib"))) {
    stop("Could not symlink the real MKL library directory for thread API testing",
         call. = FALSE)
  }
  prefix
}

missing_thread_api_prefix <- make_missing_thread_api_prefix()
space_prefix <- file.path(work_root, "mkl prefix")
dir.create(space_prefix)
if (!file.symlink(include_dir, file.path(space_prefix, "include")) ||
    !file.symlink(mkl_lib_dir, file.path(space_prefix, "lib"))) {
  stop("Could not create the MKL path-with-space test fixture", call. = FALSE)
}

copy_source_tree <- function(destination) {
  for (relative_path in source_files) {
    source_path <- file.path(repo_root, relative_path)
    destination_path <- file.path(destination, relative_path)
    dir.create(dirname(destination_path), recursive = TRUE,
               showWarnings = FALSE)
    if (!file.copy(source_path, destination_path, overwrite = TRUE,
                   copy.date = TRUE)) {
      stop("Could not copy configure input: ", relative_path, call. = FALSE)
    }
  }
}

summary_path <- file.path(output_dir, "mkl-configure-summary.tsv")
summary_rows <- character()
write_summary <- function() {
  writeLines(c("case\texit\texpected", summary_rows), summary_path)
}

configure_line <- function(text, prefix) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  lines[startsWith(lines, prefix)]
}

read_makevar <- function(path, name) {
  lines <- readLines(path, warn = FALSE)
  pattern <- paste0("^", name, "[[:space:]]*:?=[[:space:]]*")
  matches <- grep(pattern, lines)
  if (length(matches) != 1L) {
    stop("Expected exactly one ", name, " assignment in ", path,
         call. = FALSE)
  }
  value <- sub(pattern, "", lines[[matches[[1L]]]])
  next_line <- matches[[1L]] + 1L
  while (grepl("\\\\[[:space:]]*$", value) && next_line <= length(lines)) {
    value <- paste(value, lines[[next_line]])
    next_line <- next_line + 1L
  }
  value
}

check_makevars <- function(source_dir, enabled) {
  makevars_path <- file.path(source_dir, "src", "Makevars")
  if (!file.exists(makevars_path)) {
    stop("configure did not generate src/Makevars", call. = FALSE)
  }
  cppflags <- read_makevar(makevars_path, "FIDO_MKL_CPPFLAGS")
  libs <- read_makevar(makevars_path, "FIDO_MKL_LIBS")
  if (!enabled) {
    if (nzchar(trimws(cppflags)) || nzchar(trimws(libs))) {
      stop("Disabled configure case generated non-empty MKL Makevars", call. = FALSE)
    }
    return(invisible(NULL))
  }
  combined <- tolower(paste(cppflags, libs))
  required_libs <- c("mkl_intel_lp64", "mkl_gnu_thread", "mkl_core")
  if (!grepl("-I", cppflags, fixed = TRUE) || !grepl("/include", cppflags, fixed = TRUE) ||
      !all(vapply(required_libs, grepl, logical(1L), x = combined,
                  fixed = TRUE)) ||
      !grepl("rpath", combined, fixed = TRUE)) {
    stop("Enabled configure case generated incomplete MKL Makevars", call. = FALSE)
  }
  if (grepl("mkl_rt", combined, fixed = TRUE) ||
      grepl("intel_thread", combined, fixed = TRUE)) {
    stop("Enabled configure case selected a forbidden MKL threading library",
         call. = FALSE)
  }
  invisible(NULL)
}

run_case <- function(case_name, option, mklroot_mode, mklroot_value,
                     pkg_path, expected, enabled, required_status,
                     require_probe_evidence = FALSE,
                     require_missing_thread_api = FALSE) {
  source_dir <- file.path(work_root, paste0("source-", case_name))
  dir.create(source_dir, recursive = TRUE)
  copy_source_tree(source_dir)
  case_tmp <- file.path(source_dir, "tmp")
  dir.create(case_tmp, recursive = TRUE)
  stdout_path <- file.path(output_dir, paste0(case_name, ".stdout.log"))
  stderr_path <- file.path(output_dir, paste0(case_name, ".stderr.log"))
  env_args <- c(
    paste0("R_MAKEVARS_USER=", makevars_user),
    paste0("R_MAKEVARS_SITE=", makevars_site),
    paste0("PKG_CONFIG_PATH=", pkg_path),
    paste0("PKG_CONFIG_LIBDIR=", pkg_path),
    paste0("PKG_CONFIG=", pkg_config_command),
    paste0("TMPDIR=", case_tmp)
  )
  if (identical(mklroot_mode, "set")) {
    env_args <- c(env_args, paste0("MKLROOT=", mklroot_value))
  } else {
    env_args <- c("-u", "MKLROOT", env_args)
  }
  configure_args <- c("sh", "configure")
  if (nzchar(option)) configure_args <- c(configure_args, option)
  old_dir <- setwd(source_dir)
  status <- tryCatch(
    system2("env", args = shQuote(c(env_args, configure_args)),
            stdout = stdout_path, stderr = stderr_path),
    error = function(error) {
      writeLines(conditionMessage(error), stderr_path)
      127L
    }
  )
  setwd(old_dir)
  if (file.exists(file.path(source_dir, "config.log"))) {
    file.copy(file.path(source_dir, "config.log"),
              file.path(output_dir, paste0(case_name, ".config.log")), overwrite = TRUE)
  }
  if (!length(status) || is.null(status)) status <- 0L
  status <- as.integer(status[[1L]])
  output <- paste(c(readLines(stdout_path, warn = FALSE),
                    readLines(stderr_path, warn = FALSE)), collapse = "\n")
  enabled_lines <- configure_line(output, "configure: MKL enabled:")
  disabled_lines <- configure_line(output, "configure: MKL disabled:")
  error_lines <- configure_line(output, "configure: error:")
  if (required_status == 0L && status != 0L) {
    stop("Case ", case_name, " exited ", status, " instead of 0; logs: ",
         stdout_path, ", ", stderr_path, call. = FALSE)
  }
  if (required_status != 0L && status == 0L) {
    stop("Case ", case_name, " unexpectedly succeeded; logs: ",
         stdout_path, ", ", stderr_path, call. = FALSE)
  }
  if (enabled) {
    if (length(enabled_lines) != 1L || !identical(enabled_lines[[1L]], expected) ||
        length(disabled_lines)) {
      stop("Case ", case_name, " did not print the required enabled line; logs: ",
           stdout_path, ", ", stderr_path, call. = FALSE)
    }
    check_makevars(source_dir, enabled = TRUE)
  } else if (required_status == 0L) {
    if (length(disabled_lines) != 1L ||
        !identical(disabled_lines[[1L]], expected) || length(enabled_lines)) {
      stop("Case ", case_name, " did not print the required disabled line; logs: ",
           stdout_path, ", ", stderr_path, call. = FALSE)
    }
    check_makevars(source_dir, enabled = FALSE)
  } else if (length(error_lines) != 1L || !identical(error_lines[[1L]], expected)) {
    stop("Case ", case_name, " did not print the required error line; logs: ",
         stdout_path, ", ", stderr_path, call. = FALSE)
  }
  if (isTRUE(require_probe_evidence)) {
    config_log <- file.path(output_dir, paste0(case_name, ".config.log"))
    if (!file.exists(config_log)) {
      stop("Case ", case_name, " did not produce config.log probe evidence",
           call. = FALSE)
    }
    config_lines <- readLines(config_log, warn = FALSE)
    probe_compiler_lines <- config_lines[
      grepl("fido-mkl-probe\\.cpp", config_lines) &
        grepl("g\\+\\+|clang\\+\\+|c\\+\\+", config_lines)
    ]
    if (!length(probe_compiler_lines)) {
      stop("Case ", case_name,
           " config.log lacks fido-mkl-probe.cpp compiler evidence",
           call. = FALSE)
    }
  }
  if (isTRUE(require_missing_thread_api)) {
    config_log <- file.path(output_dir, paste0(case_name, ".config.log"))
    config_lines <- readLines(config_log, warn = FALSE)
    api_error_lines <- config_lines[
      grepl("fido_missing_mkl_thread_api", config_lines, fixed = TRUE) &
        grepl("not declared|undeclared", config_lines, ignore.case = TRUE)
    ]
    if (!length(api_error_lines)) {
      stop("Case ", case_name,
           " config.log lacks the expected undeclared MKL thread API diagnostic",
           call. = FALSE)
    }
  }
  list(status = status, stdout = stdout_path, stderr = stderr_path)
}

record_case <- function(case_name, result, expected) {
  summary_rows <<- c(summary_rows, paste(case_name, result$status, expected,
                                         sep = "\t"))
  write_summary()
}

cases <- list(
  list(name = "valid-mklroot-auto", option = "", mklroot_mode = "set",
       mklroot_value = mkl_prefix, pkg_path = empty_pkg_config_dir,
       expected = "configure: MKL enabled: LP64, GNU threading (MKLROOT)",
       enabled = TRUE, required_status = 0L),
  list(name = "bogus-mklroot-auto", option = "", mklroot_mode = "set",
       mklroot_value = file.path(work_root, "does-not-exist"),
       pkg_path = empty_pkg_config_dir,
       expected = "configure: MKL disabled: MKLROOT directory does not exist",
       enabled = FALSE, required_status = 0L),
  list(name = "no-mkl-auto", option = "", mklroot_mode = "unset",
       mklroot_value = "", pkg_path = empty_pkg_config_dir,
       expected = "configure: MKL disabled: not found via MKLROOT or pkg-config",
       enabled = FALSE,
       required_status = 0L),
  list(name = "valid-mklroot-no", option = "--with-mkl=no",
       mklroot_mode = "set", mklroot_value = mkl_prefix,
       pkg_path = empty_pkg_config_dir,
       expected = "configure: MKL disabled: disabled by --with-mkl=no",
       enabled = FALSE, required_status = 0L),
  list(name = "explicit-valid-path", option = paste0("--with-mkl=", mkl_prefix),
       mklroot_mode = "set", mklroot_value = file.path(work_root, "bogus-mklroot"),
       pkg_path = empty_pkg_config_dir,
       expected = "configure: MKL enabled: LP64, GNU threading (--with-mkl)",
       enabled = TRUE, required_status = 0L),
  list(name = "explicit-space-path", option = paste0("--with-mkl=", space_prefix),
       mklroot_mode = "unset", mklroot_value = "", pkg_path = empty_pkg_config_dir,
       expected = "configure: MKL enabled: LP64, GNU threading (--with-mkl)",
       enabled = TRUE, required_status = 0L),
  list(name = "explicit-bad-path", option = paste0("--with-mkl=", work_root,
                                                    "/bad-path"),
       mklroot_mode = "unset", mklroot_value = "", pkg_path = empty_pkg_config_dir,
       expected = "configure: error: MKL requested but unavailable: --with-mkl directory does not exist",
       enabled = FALSE,
       required_status = 1L),
  list(name = "pkg-config-only", option = "", mklroot_mode = "unset",
       mklroot_value = "", pkg_path = pkg_config_dir,
       expected = "configure: MKL enabled: LP64, GNU threading (pkg-config)",
       enabled = TRUE, required_status = 0L),
  list(name = "broken-probe-auto", option = "", mklroot_mode = "set",
       mklroot_value = broken_prefix, pkg_path = empty_pkg_config_dir,
       expected = "configure: MKL disabled: compile, link, or R load/run probe failed; see config.log",
       enabled = FALSE, required_status = 0L, require_probe_evidence = TRUE),
  list(name = "broken-probe-explicit", option = paste0("--with-mkl=", broken_prefix),
       mklroot_mode = "unset", mklroot_value = "", pkg_path = empty_pkg_config_dir,
       expected = "configure: error: MKL requested but unavailable: compile, link, or R load/run probe failed; see config.log",
       enabled = FALSE, required_status = 1L, require_probe_evidence = TRUE),
  list(name = "missing-thread-api-auto", option = "", mklroot_mode = "set",
       mklroot_value = missing_thread_api_prefix, pkg_path = empty_pkg_config_dir,
       expected = "configure: MKL disabled: compile, link, or R load/run probe failed; see config.log",
       enabled = FALSE, required_status = 0L, require_probe_evidence = TRUE,
       require_missing_thread_api = TRUE),
  list(name = "missing-thread-api-explicit",
       option = paste0("--with-mkl=", missing_thread_api_prefix),
       mklroot_mode = "unset", mklroot_value = "", pkg_path = empty_pkg_config_dir,
       expected = "configure: error: MKL requested but unavailable: compile, link, or R load/run probe failed; see config.log",
       enabled = FALSE, required_status = 1L, require_probe_evidence = TRUE,
       require_missing_thread_api = TRUE)
)

if (!nzchar(pkg_config_dir)) {
  stop("No MKL pkg-config metadata found below MKL_PREFIX", call. = FALSE)
}

for (case in cases) {
  result <- run_case(case$name, case$option, case$mklroot_mode,
                     case$mklroot_value, case$pkg_path, case$expected,
                     case$enabled, case$required_status,
                     case$require_probe_evidence,
                     case$require_missing_thread_api)
  record_case(case$name, result, case$expected)
}

cat("MKL configure check passed for", length(cases), "case(s)\n")
cat("Summary:", summary_path, "\n")
