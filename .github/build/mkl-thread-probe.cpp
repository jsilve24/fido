#include <MklThreads.h>
#include <omp.h>

extern "C" void fido_mkl_thread_probe(int* status, int* workers) {
  int failures = 0;
  int team_size = 0;
  #pragma omp parallel num_threads(2) reduction(+:failures) reduction(max:team_size)
  {
    team_size = omp_get_num_threads();
#ifdef FIDO_USE_MKL
    const int original_limit = mkl_set_num_threads_local(3);
    {
      fido_threads::MklThreadLimit outer_limit;
      if (mkl_get_max_threads() != 1) ++failures;
      {
        fido_threads::MklThreadLimit inner_limit;
        if (mkl_get_max_threads() != 1) ++failures;
      }
      if (mkl_get_max_threads() != 1) ++failures;
    }
    if (mkl_set_num_threads_local(0) != 3) ++failures;
    {
      fido_threads::MklThreadLimit default_limit;
      if (mkl_get_max_threads() != 1) ++failures;
    }
    if (mkl_set_num_threads_local(original_limit) != 0) ++failures;
#else
    fido_threads::MklThreadLimit default_limit;
#endif
  }
  *status = failures;
  *workers = team_size;
}
