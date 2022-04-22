//RcppArmadillo uses all threads in Ubuntu, no matter what I change in the makefile. I tried setting the configs
//here manually, but that didn't work either. Ignoring the issue for the moment
/*
#ifdef ARMA_DONT_USE_OPENMP
#pragma message("Not using OpenMP.")
#endif

#ifdef ARMA_USE_OPENMP
#undef ARMA_USE_OPENMP
#define ARMA_DONT_USE_OPENMP
#pragma message("Undefined ARMA_USE_OPENMP.")
#endif
*/