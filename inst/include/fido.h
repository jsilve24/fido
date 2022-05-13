// File used to import others in order

// On Windows do not read autoconf-updated header
#if defined(WIN32) || defined(_WIN32)
  // R can be built with its own Rlapack library, or use an external
  // one. Only the latter has zgesdd, a complex-valued SVD using divide-and-conquer
  // on Windows we do not assume ZGESDD
  /* #define ARMA_CRIPPLED_LAPACK 1 */
  // on Windows we can now assume OpenMP with Rtools / gcc 4.9.3
  // note that performance is said to still be poor
  // cf https://cran.r-project.org/doc/manuals/r-devel/R-admin.html#The-MinGW_002dw64-toolchain
  #define ARMA_USE_OPENMP
#else
  // on the other OSs we test via LAPACK_LIBS (in configure) which
  // updates this include file
  #include <fidoGenerated.h>
#endif


#ifdef FIDO_USE_MKL // requres openmp support
  #define FIDO_USE_PARALLEL
  #define EIGEN_USE_MKL_ALL
  // #define EIGEN_DONT_PARALLELIZE
#else
  #ifdef ARMA_USE_OPENMP
    #define FIDO_USE_PARALLEL
    #include <omp.h>
  #endif
#endif 
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(openmp)]]

#include "MatrixAlgebra.h"
#include "MatDist_thread.h"
#include "MatDist.h"
#include "MultDirichletBoot.h"
#include "SpecialFunctions.h"
#include "LaplaceApproximation.h"
#include "PibbleCollapsed.h"
#include "MaltipooCollapsed.h"
#include "AdamOptim.h"
