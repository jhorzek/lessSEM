test_that("testing C++ general purpose optimization", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("Rcpp")
  testthat::skip_if_not_installed("glmnet")
  testthat::skip_if_not_installed("ncvreg")
  library(Rcpp)
  
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

  // Dirk Eddelbuettel at
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

library(lessSEM)
library(glmnet)
library(ncvreg)
set.seed(123)

# first, we simulate data for our
# linear regression.
N <- 100 # number of persons
p <- 10 # number of predictors
X <- matrix(stats::rnorm(N*p),	nrow = N, ncol = p) # design matrix
b <- c(rep(1,4),
       rep(0,6)) # true regression weights
y <- X%*%matrix(b,ncol = 1) + stats::rnorm(N,0,.2)

data <- list("y" = y,
             "X" = cbind(1,X))
parameters <- rep(0, ncol(data$X))
names(parameters) <- paste0("b", 0:(length(parameters)-1))

# optimize
regularized <- paste0("b", 1:(length(b)))
lambdas <- seq(0,1,.1)
ridgePen <- gpRidgeCpp(
  par = parameters,
  regularized = regularized,
  fn = ffp, 
  gr = gfp,
  lambdas = lambdas,
  additionalArguments = data)

# for comparison:
fittingFunction <- function(par, y, X, N, lambda){
  pred <- X %*% matrix(par, ncol = 1)
  sse <- sum((y - pred)^2)
  return((.5/N)*sse + lambda * sum(par[2:length(par)]^2))
}

opt <- optim(par = parameters,
             fn = fittingFunction,
             y = y,
             X = cbind(1,X),
             N = N,
             lambda =  ridgePen@fits$lambda[3],
             method = "BFGS")
fitfunction(opt$par, data) +
  ridgePen@fits$lambda[3] * sum(opt$par[2:length(opt$par)]^2) -
opt$value
ridgePen@fits$regM2LL[3] - opt$value

testthat::expect_equal(
  all(round(opt$par - 
              ridgePen@parameters[3,ridgePen@parameterLabels], 4) == 0), TRUE)

## Testing lasso
Xsc <- scale(X)
out <- glmnet(x = Xsc, 
              y = y, 
              lambda = seq(0,.8,.1),
              standardize = FALSE)
lambdas <- out$lambda

# with lessSEM
data <- list("y" = y,
             "X" = cbind(1,Xsc))
lassoPen <- gpLassoCpp(
  par = parameters,
  regularized = regularized,
  fn = ffp, 
  gr = gfp,
  lambdas = lambdas,
  additionalArguments = data)

testthat::expect_equal(all(abs(t(lassoPen@parameters[,lassoPen@parameterLabels]) -
  coef(out)) <.0001), TRUE)

Xstd <- ncvreg::std(X)
attr(Xstd, "center") <- NULL
attr(Xstd, "scale") <- NULL
attr(Xstd, "nonsingular") <- NULL

mcpFit <- ncvreg(X = Xstd, y = y, penalty = "MCP", nlambda = 30)
lambdas <- mcpFit$lambda
thetas <- mcpFit$gamma
coefs <- coef(mcpFit)

data <- list("y" = y,
             "X" = cbind(1,Xstd))
mcpFitGp <- gpMcpCpp(par = parameters,
                     regularized = regularized,
                     fn = ffp, 
                     gr = gfp,
                     lambdas = lambdas,
                     thetas = thetas,
                     additionalArguments = data)


testthat::expect_equal(all(abs(t(mcpFitGp@parameters[,mcpFitGp@parameterLabels]) -
                                 coef(mcpFit)) <.001), TRUE)

# plot(mcpFit)
# plot(mcpFitGp)

scadFit <- ncvreg(X = Xstd, y = y, penalty = "SCAD", 
                  nlambda = 30, gamma = 3.4)
lambdas <- scadFit$lambda
thetas <- scadFit$gamma
coefs <- coef(scadFit)

scadFitGp <-  gpScadCpp(par = parameters,
                       regularized = regularized,
                       fn = ffp, 
                       gr = gfp,
                       lambdas = lambdas,
                       thetas = thetas,
                       additionalArguments = data)

testthat::expect_equal(all(abs(t(scadFitGp@parameters[,scadFitGp@parameterLabels]) -
                                 coef(scadFit)) <.001), TRUE)
# plot(scadFit)
# plot(scadFitGp)
})

