// File used to import others in order

#ifdef FIDO_USE_MKL // requres openmp support
  #define FIDO_USE_PARALLEL
  #define EIGEN_USE_MKL_ALL
  // #define EIGEN_DONT_PARALLELIZE
#else
  #ifdef _OPENMP
    #define FIDO_USE_PARALLEL
    #include <omp.h>
  #endif
#endif 
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]


#include "MatrixAlgebra.h"
#include "MatDist_thread.h"
#include "MatDist.h"
#include "MultDirichletBoot.h"
#include "SpecialFunctions.h"
#include "LaplaceApproximation.h"
#include "PibbleCollapsed.h"
#include "MaltipooCollapsed.h"
#include "AdamOptim.h"