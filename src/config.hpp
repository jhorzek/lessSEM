//RcppArmadillo uses all threads in Ubuntu, no matter what I change in the makefile. I tried setting the configs
//here manually, but that didn't work either. Ignoring the issue for the moment
#include <RcppArmadillo.h>
#undef ARMA_DONT_USE_OPENMP
#define ARMA_DONT_USE_OPENMP 1
#undef ARMA_USE_OPENMP
#undef ARMA_OPENMP_THREADS
#define ARMA_OPENMP_THREADS 1
