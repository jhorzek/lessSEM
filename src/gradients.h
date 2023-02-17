#ifndef GRADIENTS_H
#define GRADIENTS_H

#include <RcppArmadillo.h>
#include "SEM.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

void initializeGradients(SEMCpp& SEM, const bool raw);

arma::rowvec gradientsByGroup(const SEMCpp& SEM, bool raw);

#endif