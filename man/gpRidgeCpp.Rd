% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gpOptimizationInterfaceCpp.R
\name{gpRidgeCpp}
\alias{gpRidgeCpp}
\title{gpRidgeCpp}
\usage{
gpRidgeCpp(
  par,
  regularized,
  fn,
  gr,
  lambdas,
  additionalArguments,
  method = "glmnet",
  control = lessSEM::controlGlmnet()
)
}
\arguments{
\item{par}{labeled vector with starting values}

\item{regularized}{vector with names of parameters which are to be regularized.}

\item{fn}{R function which takes the parameters as input and returns the
fit value (a single value)}

\item{gr}{R function which takes the parameters as input and returns the
gradients of the objective function. If set to NULL, numDeriv will be used
to approximate the gradients}

\item{lambdas}{numeric vector: values for the tuning parameter lambda}

\item{additionalArguments}{list with additional arguments passed to fn and gr}

\item{method}{which optimizer should be used? Currently implemented are ista
and glmnet.}

\item{control}{used to control the optimizer. This element is generated with
the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
for more details.}
}
\value{
Object of class gpRegularized
}
\description{
Implements ridge regularization for general purpose optimization problems with C++ functions.
The penalty function is given by:
\deqn{p( x_j) = \lambda x_j^2}
Note that ridge regularization will not set any of the parameters to zero
but result in a shrinkage towards zero.
}
\details{
The interface is inspired by optim, but a bit more restrictive. Users have to supply a vector
with starting values (important: This vector \emph{must} have labels), a fitting
function, and a gradient function. These fitting functions \emph{must} take an const Rcpp::NumericVector& with parameter
values as first argument and an Rcpp::List& as second argument

Ridge regularization:
\itemize{
\item Hoerl, A. E., & Kennard, R. W. (1970). Ridge Regression: Biased Estimation
for Nonorthogonal Problems. Technometrics, 12(1), 55–67.
https://doi.org/10.1080/00401706.1970.10488634
}

For more details on GLMNET, see:
\itemize{
\item Friedman, J., Hastie, T., & Tibshirani, R. (2010).
Regularization Paths for Generalized Linear Models via Coordinate Descent.
Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
\item Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
A Comparison of Optimization Methods and Software for Large-scale
L1-regularized Linear Classification. Journal of Machine Learning Research, 11, 3183–3234.
\item Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012).
An improved GLMNET for l1-regularized logistic regression.
The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
}

For more details on ISTA, see:
\itemize{
\item Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding
Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
183–202. https://doi.org/10.1137/080716542
\item Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013).
A General Iterative Shrinkage and Thresholding Algorithm for Non-convex
Regularized Optimization Problems. Proceedings of the 30th International
Conference on Machine Learning, 28(2)(2), 37–45.
\item Parikh, N., & Boyd, S. (2013). Proximal Algorithms. Foundations and
Trends in Optimization, 1(3), 123–231.
}
}
\examples{
\donttest{
# This example shows how to use the optimizers
# for C++ objective functions. We will use
# a linear regression as an example. Note that
# this is not a useful application of the optimizers
# as there are specialized packages for linear regression
# (e.g., glmnet)

library(Rcpp)
library(lessSEM)

linreg <- '
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
double fitfunction(const Rcpp::NumericVector& parameters, Rcpp::List& data){
  // extract all required elements:
  arma::colvec b = Rcpp::as<arma::colvec>(parameters);
  arma::colvec y = Rcpp::as<arma::colvec>(data["y"]); // the dependent variable
  arma::mat X = Rcpp::as<arma::mat>(data["X"]); // the design matrix
  
  // compute the sum of squared errors:
    arma::mat sse = arma::trans(y-X*b)*(y-X*b);
    
    // other packages, such as glmnet, scale the sse with 
    // 1/(2*N), where N is the sample size. We will do that here as well
    
    sse *= 1.0/(2.0 * y.n_elem);
    
    // note: We must return a double, but the sse is a matrix
    // To get a double, just return the single value that is in 
    // this matrix:
      return(sse(0,0));
}

// [[Rcpp::export]]
arma::rowvec gradientfunction(const Rcpp::NumericVector& parameters, Rcpp::List& data){
  // extract all required elements:
  arma::colvec b = Rcpp::as<arma::colvec>(parameters);
  arma::colvec y = Rcpp::as<arma::colvec>(data["y"]); // the dependent variable
  arma::mat X = Rcpp::as<arma::mat>(data["X"]); // the design matrix
  
  // note: we want to return our gradients as row-vector; therefore,
  // we have to transpose the resulting column-vector:
    arma::rowvec gradients = arma::trans(-2.0*X.t() * y + 2.0*X.t()*X*b);
    
    // other packages, such as glmnet, scale the sse with 
    // 1/(2*N), where N is the sample size. We will do that here as well
    
    gradients *= (.5/y.n_rows);
    
    return(gradients);
}

// https://gallery.rcpp.org/articles/passing-cpp-function-pointers/
typedef double (*fitFunPtr)(const Rcpp::NumericVector&, //parameters
                Rcpp::List& //additional elements
);
typedef Rcpp::XPtr<fitFunPtr> fitFunPtr_t;

typedef arma::rowvec (*gradientFunPtr)(const Rcpp::NumericVector&, //parameters
                      Rcpp::List& //additional elements
);
typedef Rcpp::XPtr<gradientFunPtr> gradientFunPtr_t;

// [[Rcpp::export]]
fitFunPtr_t fitfunPtr() {
        return(fitFunPtr_t(new fitFunPtr(&fitfunction)));
}

// [[Rcpp::export]]
gradientFunPtr_t gradfunPtr() {
        return(gradientFunPtr_t(new gradientFunPtr(&gradientfunction)));
}
'

Rcpp::sourceCpp(code = linreg)

ffp <- fitfunPtr()
gfp <- gradfunPtr()

N <- 100 # number of persons
p <- 10 # number of predictors
X <- matrix(rnorm(N*p),	nrow = N, ncol = p) # design matrix
b <- c(rep(1,4), 
       rep(0,6)) # true regression weights
y <- X\%*\%matrix(b,ncol = 1) + rnorm(N,0,.2)

data <- list("y" = y,
             "X" = cbind(1,X))
parameters <- rep(0, ncol(data$X))
names(parameters) <- paste0("b", 0:(length(parameters)-1))

r <- gpRidgeCpp(par = parameters, 
                 regularized = paste0("b", 1:(length(b)-1)),
                 fn = ffp, 
                 gr = gfp, 
                 lambdas = seq(0,1,.1), 
                 additionalArguments = data)

r@parameters
}
}
