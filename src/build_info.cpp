#include <fido.h>
#include <RcppEigen.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

#ifdef FIDO_USE_MKL
#include <mkl.h>
#endif

// [[Rcpp::export(rng = false)]]
Rcpp::List fido_build_info_cpp() {
#if defined(_OPENMP)
  const int openmp_max_threads = omp_get_max_threads();
#else
  const int openmp_max_threads = NA_INTEGER;
#endif

#ifdef FIDO_USE_MKL
  const int mkl_max_threads = mkl_get_max_threads();
#else
  const int mkl_max_threads = NA_INTEGER;
#endif

  return Rcpp::List::create(
    Rcpp::_["openmp_compiled"] = static_cast<bool>(
#if defined(_OPENMP)
      true
#else
      false
#endif
    ),
    Rcpp::_["fido_use_parallel"] = static_cast<bool>(
#ifdef FIDO_USE_PARALLEL
      true
#else
      false
#endif
    ),
    Rcpp::_["fido_use_mkl"] = static_cast<bool>(
#ifdef FIDO_USE_MKL
      true
#else
      false
#endif
    ),
    Rcpp::_["openmp_max_threads"] = openmp_max_threads,
    Rcpp::_["eigen_max_threads"] = Eigen::nbThreads(),
    Rcpp::_["mkl_max_threads"] = mkl_max_threads
  );
}
