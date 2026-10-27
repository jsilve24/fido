#include <mkl.h>
#include <Eigen/Dense>
#include <cmath>

static_assert(sizeof(MKL_INT) == 4, "The probe requires LP64");

extern "C" void fido_mkl_probe(int* status, int* version) {
  MKLVersion runtime_version;
  mkl_get_version(&runtime_version);
  version[0] = runtime_version.MajorVersion;
  version[1] = runtime_version.MinorVersion;
  version[2] = runtime_version.UpdateVersion;
  double matrix[] = {4.0, 0.0, 0.0, 9.0};
  *status = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'U', 2, matrix, 2);
  if (*status != 0) return;
  if (matrix[0] != 2.0 || matrix[3] != 3.0) {
    *status = -1;
    return;
  }
  Eigen::MatrixXd covariance = Eigen::MatrixXd::Identity(3, 3);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(covariance);
  if (solver.info() != Eigen::Success ||
      (solver.eigenvalues().array() - 1.0).abs().maxCoeff() > 1e-12) {
    *status = -2;
  }
}
