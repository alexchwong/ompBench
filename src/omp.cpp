#include "omp.h"

unsigned int use_threads_omp(int n_threads) {
  #ifdef _OPENMP
  if(n_threads <= 1) return(1);
  if(n_threads > omp_get_max_threads()) {
    return((unsigned int)omp_get_max_threads());
  }
  return((unsigned int)n_threads);
  #else
  return(1);
  #endif
}
