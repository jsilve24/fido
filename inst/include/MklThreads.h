#ifndef FIDO_MKL_THREADS_H
#define FIDO_MKL_THREADS_H

#ifdef FIDO_USE_MKL
#include <mkl.h>
#endif

namespace fido_threads {
class MklThreadLimit {
public:
  MklThreadLimit()
#ifdef FIDO_USE_MKL
    : previous_limit(mkl_set_num_threads_local(1))
#endif
  {}

  ~MklThreadLimit() {
#ifdef FIDO_USE_MKL
    mkl_set_num_threads_local(previous_limit);
#endif
  }

  MklThreadLimit(const MklThreadLimit&) = delete;
  MklThreadLimit& operator=(const MklThreadLimit&) = delete;

private:
#ifdef FIDO_USE_MKL
  int previous_limit;
#endif
};
}

#endif
