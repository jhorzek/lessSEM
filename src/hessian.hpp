#ifndef HESSIAN_H
#define HESSIAN_H

#include <RcppArmadillo.h>
#include "config.hpp"
#include "SEM.hpp"

// [[Rcpp :: depends ( RcppArmadillo )]]

arma::mat approximateHessian(SEMCpp& SEM,
                             Rcpp::StringVector label_,
                             arma::vec value_,
                             bool raw,
                             double eps);

#endif