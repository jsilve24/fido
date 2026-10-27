#' Report fido native build and runtime information
#'
#' Thread counts are the current runtime limits, not a count of active OS
#' threads. R's BLAS is separate from the backend used by Eigen in fido.
#'
#' @return A named list with stable fields:
#' \describe{
#'   \item{openmp_compiled}{Logical scalar indicating whether the native code
#'   was compiled with the `_OPENMP` macro.}
#'   \item{fido_use_parallel}{Logical scalar indicating whether `fido.h`
#'   defines `FIDO_USE_PARALLEL`.}
#'   \item{fido_use_mkl}{Logical scalar indicating whether `fido.h` defines
#'   `FIDO_USE_MKL`.}
#'   \item{openmp_max_threads}{Maximum OpenMP threads at runtime, or `NA`
#'   when `_OPENMP` is unavailable.}
#'   \item{eigen_max_threads}{Maximum Eigen threads reported at runtime.}
#'   \item{mkl_max_threads}{Maximum MKL threads at runtime, or `NA` when MKL
#'   is not enabled.}
#'   \item{r_blas}{The BLAS reported by R's `extSoftVersion()`, or `NA` if R
#'   does not report one.}
#' }
#' @export
fido_build_info <- function() {
  info <- fido_build_info_cpp()
  info$r_blas <- unname(extSoftVersion()["BLAS"])
  info
}
