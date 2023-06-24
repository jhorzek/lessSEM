#' gpLassoCpp 
#' 
#' Implements lasso regularization for general purpose optimization problems with C++ functions.
#' The penalty function is given by:
#' \deqn{p( x_j) = \lambda |x_j|}
#' Lasso regularization will set parameters to zero if \eqn{\lambda} is large enough
#' 
#' The interface is inspired by optim, but a bit more restrictive. Users have to supply a vector 
#' with starting values (important: This vector _must_ have labels), a fitting
#' function, and a gradient function. These fitting functions _must_ take an const Rcpp::NumericVector& with parameter
#' values as first argument and an Rcpp::List& as second argument
#' 
#' Lasso regularization:
#' 
#' * Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. Journal of the Royal Statistical 
#' Society. Series B (Methodological), 58(1), 267–288.
#'  
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classification. Journal of Machine Learning Research, 11, 3183–3234.
#' * Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
#' An improved GLMNET for l1-regularized logistic regression. 
#' The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
#' 
#' For more details on ISTA, see:
#' 
#' * Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
#' 183–202. https://doi.org/10.1137/080716542
#' * Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
#' Regularized Optimization Problems. Proceedings of the 30th International 
#' Conference on Machine Learning, 28(2)(2), 37–45.
#' * Parikh, N., & Boyd, S. (2013). Proximal Algorithms. Foundations and 
#' Trends in Optimization, 1(3), 123–231.
#'  
#' @param par labeled vector with starting values
#' @param regularized vector with names of parameters which are to be regularized.
#' @param fn pointer to Rcpp function which takes the parameters as input and returns the 
#' fit value (a single value)
#' @param gr pointer to Rcpp function which takes the parameters as input and returns the 
#' gradients of the objective function.
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param nLambdas alternative to lambda: If alpha = 1, lessSEM can automatically
#' compute the first lambda value which sets all regularized parameters to zero.
#' It will then generate nLambda values between 0 and the computed lambda.
#' @param curve Allows for unequally spaced lambda steps (e.g., .01,.02,.05,1,5,20). 
#' If curve is close to 1 all lambda values will be equally spaced, if curve is large 
#' lambda values will be more concentrated close to 0. See ?lessSEM::curveLambda for more information.
#' @param additionalArguments list with additional arguments passed to fn and gr
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. 
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
#' @returns Object of class gpRegularized

#' @examples 
#' \donttest{
#' # This example shows how to use the optimizers
#' # for C++ objective functions. We will use
#' # a linear regression as an example. Note that
#' # this is not a useful application of the optimizers
#' # as there are specialized packages for linear regression
#' # (e.g., glmnet)
#' 
#' library(Rcpp)
#' library(lessSEM)
#' 
#' linreg <- '
#' // [[Rcpp::depends(RcppArmadillo)]]
#' #include <RcppArmadillo.h>
#' 
#' // [[Rcpp::export]]
#' double fitfunction(const Rcpp::NumericVector& parameters, Rcpp::List& data){
#'   // extract all required elements:
#'   arma::colvec b = Rcpp::as<arma::colvec>(parameters);
#'   arma::colvec y = Rcpp::as<arma::colvec>(data["y"]); // the dependent variable
#'   arma::mat X = Rcpp::as<arma::mat>(data["X"]); // the design matrix
#'   
#'   // compute the sum of squared errors:
#'     arma::mat sse = arma::trans(y-X*b)*(y-X*b);
#'     
#'     // other packages, such as glmnet, scale the sse with 
#'     // 1/(2*N), where N is the sample size. We will do that here as well
#'     
#'     sse *= 1.0/(2.0 * y.n_elem);
#'     
#'     // note: We must return a double, but the sse is a matrix
#'     // To get a double, just return the single value that is in 
#'     // this matrix:
#'       return(sse(0,0));
#' }
#' 
#' // [[Rcpp::export]]
#' arma::rowvec gradientfunction(const Rcpp::NumericVector& parameters, Rcpp::List& data){
#'   // extract all required elements:
#'   arma::colvec b = Rcpp::as<arma::colvec>(parameters);
#'   arma::colvec y = Rcpp::as<arma::colvec>(data["y"]); // the dependent variable
#'   arma::mat X = Rcpp::as<arma::mat>(data["X"]); // the design matrix
#'   
#'   // note: we want to return our gradients as row-vector; therefore,
#'   // we have to transpose the resulting column-vector:
#'     arma::rowvec gradients = arma::trans(-2.0*X.t() * y + 2.0*X.t()*X*b);
#'     
#'     // other packages, such as glmnet, scale the sse with 
#'     // 1/(2*N), where N is the sample size. We will do that here as well
#'     
#'     gradients *= (.5/y.n_rows);
#'     
#'     return(gradients);
#' }
#' 
#' // Dirk Eddelbuettel at
#' // https://gallery.rcpp.org/articles/passing-cpp-function-pointers/
#' typedef double (*fitFunPtr)(const Rcpp::NumericVector&, //parameters
#'                 Rcpp::List& //additional elements
#' );
#' typedef Rcpp::XPtr<fitFunPtr> fitFunPtr_t;
#' 
#' typedef arma::rowvec (*gradientFunPtr)(const Rcpp::NumericVector&, //parameters
#'                       Rcpp::List& //additional elements
#' );
#' typedef Rcpp::XPtr<gradientFunPtr> gradientFunPtr_t;
#' 
#' // [[Rcpp::export]]
#' fitFunPtr_t fitfunPtr() {
#'         return(fitFunPtr_t(new fitFunPtr(&fitfunction)));
#' }
#' 
#' // [[Rcpp::export]]
#' gradientFunPtr_t gradfunPtr() {
#'         return(gradientFunPtr_t(new gradientFunPtr(&gradientfunction)));
#' }
#' '
#' 
#' Rcpp::sourceCpp(code = linreg)
#' 
#' ffp <- fitfunPtr()
#' gfp <- gradfunPtr()
#' 
#' N <- 100 # number of persons
#' p <- 10 # number of predictors
#' X <- matrix(rnorm(N*p),	nrow = N, ncol = p) # design matrix
#' b <- c(rep(1,4), 
#'        rep(0,6)) # true regression weights
#' y <- X%*%matrix(b,ncol = 1) + rnorm(N,0,.2)
#' 
#' data <- list("y" = y,
#'              "X" = cbind(1,X))
#' parameters <- rep(0, ncol(data$X))
#' names(parameters) <- paste0("b", 0:(length(parameters)-1))
#' 
#' l1 <- gpLassoCpp(par = parameters, 
#'                  regularized = paste0("b", 1:(length(b)-1)),
#'                  fn = ffp, 
#'                  gr = gfp, 
#'                  lambdas = seq(0,1,.1), 
#'                  additionalArguments = data)
#' 
#' l1@parameters
#' 
#' }
#' @export
gpLassoCpp <- function(par,
                    regularized,
                    fn,
                    gr,
                    lambdas = NULL,
                    nLambdas = NULL,
                    curve = 1,
                    additionalArguments,
                    method = "glmnet", 
                    control = lessSEM::controlGlmnet()
){
  
  if(!is(fn, "externalptr") | !is(gr, "externalptr")){
    stop("fn and gr must be pointers to C++ functions.")
  }
  if(!is(additionalArguments, "list")){
    stop("additionalArguments must be of class list.")
  }
  
  
  if(is.null(lambdas) && is.null(nLambdas)){
    stop("Specify either lambdas or nLambdas")
  }
  
  if(!is.null(nLambdas)){
    tuningParameters <- data.frame(nLambdas = nLambdas,
                                   curve = curve)
  }else{
    tuningParameters <- data.frame(lambda = lambdas,
                                   alpha = 1,
                                   theta = 0)
  }
  
  result <- .gpOptimizationInternal(par = par,
                                            fn = fn,
                                            gr = gr,
                                            additionalArguments = additionalArguments,
                                            isCpp = TRUE,
                                            penalty = "lasso", 
                                            weights = regularized,
                                            tuningParameters = tuningParameters, 
                                            method = method,
                                            control = control
  )
  
  return(result)
  
}

#' gpAdaptiveLassoCpp
#' 
#' Implements adaptive lasso regularization for general purpose optimization problems with C++ functions.
#' The penalty function is given by:
#' \deqn{p( x_j) = p( x_j) = \frac{1}{w_j}\lambda| x_j|}
#' Adaptive lasso regularization will set parameters to zero if \eqn{\lambda} is large enough.
#' 
#' The interface is inspired by optim, but a bit more restrictive. Users have to supply a vector 
#' with starting values (important: This vector _must_ have labels), a fitting
#' function, and a gradient function. These fitting functions _must_ take an const Rcpp::NumericVector& with parameter
#' values as first argument and an Rcpp::List& as second argument
#' 
#' Adaptive lasso regularization:
#' 
#' * Zou, H. (2006). The adaptive lasso and its oracle properties. Journal of the American Statistical Association, 
#' 101(476), 1418–1429. https://doi.org/10.1198/016214506000000735
#'  
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classification. Journal of Machine Learning Research, 11, 3183–3234.
#' * Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
#' An improved GLMNET for l1-regularized logistic regression. 
#' The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
#' 
#' For more details on ISTA, see:
#' 
#' * Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
#' 183–202. https://doi.org/10.1137/080716542
#' * Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
#' Regularized Optimization Problems. Proceedings of the 30th International 
#' Conference on Machine Learning, 28(2)(2), 37–45.
#' * Parikh, N., & Boyd, S. (2013). Proximal Algorithms. Foundations and 
#' Trends in Optimization, 1(3), 123–231.
#'  
#' @param par labeled vector with starting values
#' @param regularized vector with names of parameters which are to be regularized.
#' @param fn R function which takes the parameters as input and returns the 
#' fit value (a single value)
#' @param gr R function which takes the parameters as input and returns the 
#' gradients of the objective function. If set to NULL, numDeriv will be used
#' to approximate the gradients 
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param nLambdas alternative to lambda: If alpha = 1, lessSEM can automatically
#' compute the first lambda value which sets all regularized parameters to zero.
#' It will then generate nLambda values between 0 and the computed lambda.
#' @param curve Allows for unequally spaced lambda steps (e.g., .01,.02,.05,1,5,20). 
#' If curve is close to 1 all lambda values will be equally spaced, if curve is large 
#' lambda values will be more concentrated close to 0. See ?lessSEM::curveLambda for more information.
#' @param weights labeled vector with adaptive lasso weights. NULL will use 1/abs(par)
#' @param additionalArguments list with additional arguments passed to fn and gr
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. 
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
#' @returns Object of class gpRegularized

#' @examples 
#' \donttest{
#' # This example shows how to use the optimizers
#' # for C++ objective functions. We will use
#' # a linear regression as an example. Note that
#' # this is not a useful application of the optimizers
#' # as there are specialized packages for linear regression
#' # (e.g., glmnet)
#' 
#' library(Rcpp)
#' library(lessSEM)
#' 
#' linreg <- '
#' // [[Rcpp::depends(RcppArmadillo)]]
#' #include <RcppArmadillo.h>
#' 
#' // [[Rcpp::export]]
#' double fitfunction(const Rcpp::NumericVector& parameters, Rcpp::List& data){
#'   // extract all required elements:
#'   arma::colvec b = Rcpp::as<arma::colvec>(parameters);
#'   arma::colvec y = Rcpp::as<arma::colvec>(data["y"]); // the dependent variable
#'   arma::mat X = Rcpp::as<arma::mat>(data["X"]); // the design matrix
#'   
#'   // compute the sum of squared errors:
#'     arma::mat sse = arma::trans(y-X*b)*(y-X*b);
#'     
#'     // other packages, such as glmnet, scale the sse with 
#'     // 1/(2*N), where N is the sample size. We will do that here as well
#'     
#'     sse *= 1.0/(2.0 * y.n_elem);
#'     
#'     // note: We must return a double, but the sse is a matrix
#'     // To get a double, just return the single value that is in 
#'     // this matrix:
#'       return(sse(0,0));
#' }
#' 
#' // [[Rcpp::export]]
#' arma::rowvec gradientfunction(const Rcpp::NumericVector& parameters, Rcpp::List& data){
#'   // extract all required elements:
#'   arma::colvec b = Rcpp::as<arma::colvec>(parameters);
#'   arma::colvec y = Rcpp::as<arma::colvec>(data["y"]); // the dependent variable
#'   arma::mat X = Rcpp::as<arma::mat>(data["X"]); // the design matrix
#'   
#'   // note: we want to return our gradients as row-vector; therefore,
#'   // we have to transpose the resulting column-vector:
#'     arma::rowvec gradients = arma::trans(-2.0*X.t() * y + 2.0*X.t()*X*b);
#'     
#'     // other packages, such as glmnet, scale the sse with 
#'     // 1/(2*N), where N is the sample size. We will do that here as well
#'     
#'     gradients *= (.5/y.n_rows);
#'     
#'     return(gradients);
#' }
#' 
#' // Dirk Eddelbuettel at
#' // https://gallery.rcpp.org/articles/passing-cpp-function-pointers/
#' typedef double (*fitFunPtr)(const Rcpp::NumericVector&, //parameters
#'                 Rcpp::List& //additional elements
#' );
#' typedef Rcpp::XPtr<fitFunPtr> fitFunPtr_t;
#' 
#' typedef arma::rowvec (*gradientFunPtr)(const Rcpp::NumericVector&, //parameters
#'                       Rcpp::List& //additional elements
#' );
#' typedef Rcpp::XPtr<gradientFunPtr> gradientFunPtr_t;
#' 
#' // [[Rcpp::export]]
#' fitFunPtr_t fitfunPtr() {
#'         return(fitFunPtr_t(new fitFunPtr(&fitfunction)));
#' }
#' 
#' // [[Rcpp::export]]
#' gradientFunPtr_t gradfunPtr() {
#'         return(gradientFunPtr_t(new gradientFunPtr(&gradientfunction)));
#' }
#' '
#' 
#' Rcpp::sourceCpp(code = linreg)
#' 
#' ffp <- fitfunPtr()
#' gfp <- gradfunPtr()
#' 
#' N <- 100 # number of persons
#' p <- 10 # number of predictors
#' X <- matrix(rnorm(N*p),	nrow = N, ncol = p) # design matrix
#' b <- c(rep(1,4), 
#'        rep(0,6)) # true regression weights
#' y <- X%*%matrix(b,ncol = 1) + rnorm(N,0,.2)
#' 
#' data <- list("y" = y,
#'              "X" = cbind(1,X))
#' parameters <- rep(0, ncol(data$X))
#' names(parameters) <- paste0("b", 0:(length(parameters)-1))
#' 
#' al1 <- gpAdaptiveLassoCpp(par = parameters, 
#'                  regularized = paste0("b", 1:(length(b)-1)),
#'                  fn = ffp, 
#'                  gr = gfp, 
#'                  lambdas = seq(0,1,.1), 
#'                  additionalArguments = data)
#' 
#' al1@parameters
#' }
#' @export
gpAdaptiveLassoCpp <- function(par,
                            regularized,
                            weights = NULL,
                            fn,
                            gr,
                            lambdas = NULL,
                            nLambdas = NULL,
                            curve = 1,
                            additionalArguments,
                            method = "glmnet", 
                            control = lessSEM::controlGlmnet()){
  
  if(!is(fn, "externalptr") | !is(gr, "externalptr")){
    stop("fn and gr must be pointers to C++ functions.")
  }
  if(!is(additionalArguments, "list")){
    stop("additionalArguments must be of class list.")
  }
  
  if(is.null(weights)){
    weights <- 1/abs(par)
    weights[!names(weights) %in% regularized] <- 0
    cat("\n")
    rlang::inform(c("Note","Building weights based on par as weights = 1/abs(par)."))
  }
  
  if(! all(regularized %in% names(weights))) stop(paste0(
    "You specified that the following parameters should be regularized:\n",
    paste0(regularized, collapse = ", "), 
    ". Not all of these parameters could be found in the model.\n",
    "The model has the following parameters:\n",
    names(weights)
  ))
  
  if(is.null(lambdas) && is.null(nLambdas)){
    stop("Specify either lambdas or nLambdas")
  }
  
  if(!is.null(nLambdas)){
    tuningParameters <- data.frame(nLambdas = nLambdas,
                                   curve = curve)
  }else{
    tuningParameters <- data.frame(lambda = lambdas,
                                   alpha = 1,
                                   theta = 0)
  }
  
  result <- .gpOptimizationInternal(par = par,
                                            fn = fn,
                                            gr = gr,
                                            additionalArguments = additionalArguments,
                                            isCpp = TRUE,
                                            penalty = "adaptiveLasso", 
                                            weights = weights,
                                            tuningParameters = tuningParameters, 
                                            method = method,
                                            control = control
  )
  
  return(result)
  
}

#' gpRidgeCpp
#' 
#' Implements ridge regularization for general purpose optimization problems with C++ functions.
#' The penalty function is given by:
#' \deqn{p( x_j) = \lambda x_j^2}
#' Note that ridge regularization will not set any of the parameters to zero
#' but result in a shrinkage towards zero. 
#' 
#' The interface is inspired by optim, but a bit more restrictive. Users have to supply a vector 
#' with starting values (important: This vector _must_ have labels), a fitting
#' function, and a gradient function. These fitting functions _must_ take an const Rcpp::NumericVector& with parameter
#' values as first argument and an Rcpp::List& as second argument
#' 
#' Ridge regularization:
#' 
#' * Hoerl, A. E., & Kennard, R. W. (1970). Ridge Regression: Biased Estimation 
#' for Nonorthogonal Problems. Technometrics, 12(1), 55–67. 
#' https://doi.org/10.1080/00401706.1970.10488634
#' 
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classification. Journal of Machine Learning Research, 11, 3183–3234.
#' * Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
#' An improved GLMNET for l1-regularized logistic regression. 
#' The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
#' 
#' For more details on ISTA, see:
#' 
#' * Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
#' 183–202. https://doi.org/10.1137/080716542
#' * Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
#' Regularized Optimization Problems. Proceedings of the 30th International 
#' Conference on Machine Learning, 28(2)(2), 37–45.
#' * Parikh, N., & Boyd, S. (2013). Proximal Algorithms. Foundations and 
#' Trends in Optimization, 1(3), 123–231.
#'  
#' @param par labeled vector with starting values
#' @param regularized vector with names of parameters which are to be regularized.
#' @param fn R function which takes the parameters as input and returns the 
#' fit value (a single value)
#' @param gr R function which takes the parameters as input and returns the 
#' gradients of the objective function. If set to NULL, numDeriv will be used
#' to approximate the gradients 
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param additionalArguments list with additional arguments passed to fn and gr
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. 
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
#' @returns Object of class gpRegularized

#' @examples 
#' \donttest{
#' # This example shows how to use the optimizers
#' # for C++ objective functions. We will use
#' # a linear regression as an example. Note that
#' # this is not a useful application of the optimizers
#' # as there are specialized packages for linear regression
#' # (e.g., glmnet)
#' 
#' library(Rcpp)
#' library(lessSEM)
#' 
#' linreg <- '
#' // [[Rcpp::depends(RcppArmadillo)]]
#' #include <RcppArmadillo.h>
#' 
#' // [[Rcpp::export]]
#' double fitfunction(const Rcpp::NumericVector& parameters, Rcpp::List& data){
#'   // extract all required elements:
#'   arma::colvec b = Rcpp::as<arma::colvec>(parameters);
#'   arma::colvec y = Rcpp::as<arma::colvec>(data["y"]); // the dependent variable
#'   arma::mat X = Rcpp::as<arma::mat>(data["X"]); // the design matrix
#'   
#'   // compute the sum of squared errors:
#'     arma::mat sse = arma::trans(y-X*b)*(y-X*b);
#'     
#'     // other packages, such as glmnet, scale the sse with 
#'     // 1/(2*N), where N is the sample size. We will do that here as well
#'     
#'     sse *= 1.0/(2.0 * y.n_elem);
#'     
#'     // note: We must return a double, but the sse is a matrix
#'     // To get a double, just return the single value that is in 
#'     // this matrix:
#'       return(sse(0,0));
#' }
#' 
#' // [[Rcpp::export]]
#' arma::rowvec gradientfunction(const Rcpp::NumericVector& parameters, Rcpp::List& data){
#'   // extract all required elements:
#'   arma::colvec b = Rcpp::as<arma::colvec>(parameters);
#'   arma::colvec y = Rcpp::as<arma::colvec>(data["y"]); // the dependent variable
#'   arma::mat X = Rcpp::as<arma::mat>(data["X"]); // the design matrix
#'   
#'   // note: we want to return our gradients as row-vector; therefore,
#'   // we have to transpose the resulting column-vector:
#'     arma::rowvec gradients = arma::trans(-2.0*X.t() * y + 2.0*X.t()*X*b);
#'     
#'     // other packages, such as glmnet, scale the sse with 
#'     // 1/(2*N), where N is the sample size. We will do that here as well
#'     
#'     gradients *= (.5/y.n_rows);
#'     
#'     return(gradients);
#' }
#' 
#' // https://gallery.rcpp.org/articles/passing-cpp-function-pointers/
#' typedef double (*fitFunPtr)(const Rcpp::NumericVector&, //parameters
#'                 Rcpp::List& //additional elements
#' );
#' typedef Rcpp::XPtr<fitFunPtr> fitFunPtr_t;
#' 
#' typedef arma::rowvec (*gradientFunPtr)(const Rcpp::NumericVector&, //parameters
#'                       Rcpp::List& //additional elements
#' );
#' typedef Rcpp::XPtr<gradientFunPtr> gradientFunPtr_t;
#' 
#' // [[Rcpp::export]]
#' fitFunPtr_t fitfunPtr() {
#'         return(fitFunPtr_t(new fitFunPtr(&fitfunction)));
#' }
#' 
#' // [[Rcpp::export]]
#' gradientFunPtr_t gradfunPtr() {
#'         return(gradientFunPtr_t(new gradientFunPtr(&gradientfunction)));
#' }
#' '
#' 
#' Rcpp::sourceCpp(code = linreg)
#' 
#' ffp <- fitfunPtr()
#' gfp <- gradfunPtr()
#' 
#' N <- 100 # number of persons
#' p <- 10 # number of predictors
#' X <- matrix(rnorm(N*p),	nrow = N, ncol = p) # design matrix
#' b <- c(rep(1,4), 
#'        rep(0,6)) # true regression weights
#' y <- X%*%matrix(b,ncol = 1) + rnorm(N,0,.2)
#' 
#' data <- list("y" = y,
#'              "X" = cbind(1,X))
#' parameters <- rep(0, ncol(data$X))
#' names(parameters) <- paste0("b", 0:(length(parameters)-1))
#' 
#' r <- gpRidgeCpp(par = parameters, 
#'                  regularized = paste0("b", 1:(length(b)-1)),
#'                  fn = ffp, 
#'                  gr = gfp, 
#'                  lambdas = seq(0,1,.1), 
#'                  additionalArguments = data)
#' 
#' r@parameters
#' }
#' @export
gpRidgeCpp <- function(par,
                    regularized,
                    fn,
                    gr,
                    lambdas,
                    additionalArguments,
                    method = "glmnet", 
                    control = lessSEM::controlGlmnet()){
  
  if(!is(fn, "externalptr") | !is(gr, "externalptr")){
    stop("fn and gr must be pointers to C++ functions.")
  }
  if(!is(additionalArguments, "list")){
    stop("additionalArguments must be of class list.")
  }
  
  
  tuningParameters <- data.frame(lambda = lambdas,
                                 alpha = 0,
                                 theta = 0)
  
  result <- .gpOptimizationInternal(par = par,
                                            fn = fn,
                                            gr = gr,
                                            additionalArguments = additionalArguments,
                                            isCpp = TRUE,
                                            penalty = "ridge", 
                                            weights = regularized,
                                            tuningParameters = tuningParameters, 
                                            method = method,
                                            control = control
  )
  
  return(result)
  
}

#' gpElasticNetCpp
#' 
#' Implements elastic net regularization for general purpose optimization problems with C++ functions.
#' The penalty function is given by:
#' \deqn{p( x_j) = p( x_j) = \frac{1}{w_j}\lambda| x_j|}
#' Note that the elastic net combines ridge and lasso regularization. If \eqn{\alpha = 0}, 
#' the elastic net reduces to ridge regularization. If \eqn{\alpha = 1} it reduces
#' to lasso regularization. In between, elastic net is a compromise between the shrinkage of
#' the lasso and the ridge penalty. 
#' 
#' The interface is inspired by optim, but a bit more restrictive. Users have to supply a vector 
#' with starting values (important: This vector _must_ have labels), a fitting
#' function, and a gradient function. These fitting functions _must_ take an const Rcpp::NumericVector& with parameter
#' values as first argument and an Rcpp::List& as second argument
#' 
#' Elastic net regularization:
#' 
#' * Zou, H., & Hastie, T. (2005). Regularization and variable selection via the elastic net. 
#' Journal of the Royal Statistical Society: Series B, 67(2), 301–320. https://doi.org/10.1111/j.1467-9868.2005.00503.x
#'  
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classification. Journal of Machine Learning Research, 11, 3183–3234.
#' * Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
#' An improved GLMNET for l1-regularized logistic regression. 
#' The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
#' 
#' For more details on ISTA, see:
#' 
#' * Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
#' 183–202. https://doi.org/10.1137/080716542
#' * Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
#' Regularized Optimization Problems. Proceedings of the 30th International 
#' Conference on Machine Learning, 28(2)(2), 37–45.
#' * Parikh, N., & Boyd, S. (2013). Proximal Algorithms. Foundations and 
#' Trends in Optimization, 1(3), 123–231.
#'  
#' @param par labeled vector with starting values
#' @param regularized vector with names of parameters which are to be regularized.
#' @param fn R function which takes the parameters AND their labels 
#' as input and returns the fit value (a single value)
#' @param gr R function which takes the parameters AND their labels
#' as input and returns the gradients of the objective function. 
#' If set to NULL, numDeriv will be used to approximate the gradients 
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param alphas numeric vector with values of the tuning parameter alpha. Must be
#' between 0 and 1. 0 = ridge, 1 = lasso.
#' @param additionalArguments list with additional arguments passed to fn and gr
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. 
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
#' @returns Object of class gpRegularized

#' @examples
#' \donttest{
#' # This example shows how to use the optimizers
#' # for C++ objective functions. We will use
#' # a linear regression as an example. Note that
#' # this is not a useful application of the optimizers
#' # as there are specialized packages for linear regression
#' # (e.g., glmnet)
#' 
#' library(Rcpp)
#' library(lessSEM)
#' 
#' linreg <- '
#' // [[Rcpp::depends(RcppArmadillo)]]
#' #include <RcppArmadillo.h>
#' 
#' // [[Rcpp::export]]
#' double fitfunction(const Rcpp::NumericVector& parameters, Rcpp::List& data){
#'   // extract all required elements:
#'   arma::colvec b = Rcpp::as<arma::colvec>(parameters);
#'   arma::colvec y = Rcpp::as<arma::colvec>(data["y"]); // the dependent variable
#'   arma::mat X = Rcpp::as<arma::mat>(data["X"]); // the design matrix
#'   
#'   // compute the sum of squared errors:
#'     arma::mat sse = arma::trans(y-X*b)*(y-X*b);
#'     
#'     // other packages, such as glmnet, scale the sse with 
#'     // 1/(2*N), where N is the sample size. We will do that here as well
#'     
#'     sse *= 1.0/(2.0 * y.n_elem);
#'     
#'     // note: We must return a double, but the sse is a matrix
#'     // To get a double, just return the single value that is in 
#'     // this matrix:
#'       return(sse(0,0));
#' }
#' 
#' // [[Rcpp::export]]
#' arma::rowvec gradientfunction(const Rcpp::NumericVector& parameters, Rcpp::List& data){
#'   // extract all required elements:
#'   arma::colvec b = Rcpp::as<arma::colvec>(parameters);
#'   arma::colvec y = Rcpp::as<arma::colvec>(data["y"]); // the dependent variable
#'   arma::mat X = Rcpp::as<arma::mat>(data["X"]); // the design matrix
#'   
#'   // note: we want to return our gradients as row-vector; therefore,
#'   // we have to transpose the resulting column-vector:
#'     arma::rowvec gradients = arma::trans(-2.0*X.t() * y + 2.0*X.t()*X*b);
#'     
#'     // other packages, such as glmnet, scale the sse with 
#'     // 1/(2*N), where N is the sample size. We will do that here as well
#'     
#'     gradients *= (.5/y.n_rows);
#'     
#'     return(gradients);
#' }
#' 
#' // Dirk Eddelbuettel at
#' // https://gallery.rcpp.org/articles/passing-cpp-function-pointers/
#' typedef double (*fitFunPtr)(const Rcpp::NumericVector&, //parameters
#'                 Rcpp::List& //additional elements
#' );
#' typedef Rcpp::XPtr<fitFunPtr> fitFunPtr_t;
#' 
#' typedef arma::rowvec (*gradientFunPtr)(const Rcpp::NumericVector&, //parameters
#'                       Rcpp::List& //additional elements
#' );
#' typedef Rcpp::XPtr<gradientFunPtr> gradientFunPtr_t;
#' 
#' // [[Rcpp::export]]
#' fitFunPtr_t fitfunPtr() {
#'         return(fitFunPtr_t(new fitFunPtr(&fitfunction)));
#' }
#' 
#' // [[Rcpp::export]]
#' gradientFunPtr_t gradfunPtr() {
#'         return(gradientFunPtr_t(new gradientFunPtr(&gradientfunction)));
#' }
#' '
#' 
#' Rcpp::sourceCpp(code = linreg)
#' 
#' ffp <- fitfunPtr()
#' gfp <- gradfunPtr()
#' 
#' N <- 100 # number of persons
#' p <- 10 # number of predictors
#' X <- matrix(rnorm(N*p),	nrow = N, ncol = p) # design matrix
#' b <- c(rep(1,4), 
#'        rep(0,6)) # true regression weights
#' y <- X%*%matrix(b,ncol = 1) + rnorm(N,0,.2)
#' 
#' data <- list("y" = y,
#'              "X" = cbind(1,X))
#' parameters <- rep(0, ncol(data$X))
#' names(parameters) <- paste0("b", 0:(length(parameters)-1))
#' 
#' en <- gpElasticNetCpp(par = parameters, 
#'                  regularized = paste0("b", 1:(length(b)-1)),
#'                  fn = ffp, 
#'                  gr = gfp, 
#'                  lambdas = seq(0,1,.1), 
#'                  alphas = c(0,.5,1),
#'                  additionalArguments = data)
#' 
#' en@parameters
#' }
#' @export
gpElasticNetCpp <- function(par,
                         regularized,
                         fn,
                         gr,
                         lambdas,
                         alphas,
                         additionalArguments,
                         method = "glmnet", 
                         control = lessSEM::controlGlmnet()){
  
  if(!is(fn, "externalptr") | !is(gr, "externalptr")){
    stop("fn and gr must be pointers to C++ functions.")
  }
  if(!is(additionalArguments, "list")){
    stop("additionalArguments must be of class list.")
  }
  
  tuningParameters <- expand.grid(lambda = lambdas,
                                  alpha = alphas,
                                  theta = 0)
  
  
  result <- .gpOptimizationInternal(par = par,
                                            fn = fn,
                                            gr = gr,
                                            additionalArguments = additionalArguments,
                                            isCpp = TRUE,
                                            penalty = "elasticNet", 
                                            weights = regularized,
                                            tuningParameters = tuningParameters, 
                                            method = method,
                                            control = control
  )
  
  return(result)
}

#' gpCappedL1Cpp
#' 
#' Implements cappedL1 regularization for general purpose optimization problems with C++ functions.
#' The penalty function is given by:
#' \deqn{p( x_j) = \lambda \min(| x_j|, \theta)}
#' where \eqn{\theta > 0}. The cappedL1 penalty is identical to the lasso for 
#' parameters which are below \eqn{\theta} and identical to a constant for parameters
#' above \eqn{\theta}. As adding a constant to the fitting function will not change its
#' minimum, larger parameters can stay unregularized while smaller ones are set to zero.
#' 
#' The interface is inspired by optim, but a bit more restrictive. Users have to supply a vector 
#' with starting values (important: This vector _must_ have labels), a fitting
#' function, and a gradient function. These fitting functions _must_ take an const Rcpp::NumericVector& with parameter
#' values as first argument and an Rcpp::List& as second argument
#' 
#' CappedL1 regularization:
#' 
#' * Zhang, T. (2010). Analysis of Multi-stage Convex Relaxation for Sparse Regularization. 
#' Journal of Machine Learning Research, 11, 1081–1107.
#'  
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classification. Journal of Machine Learning Research, 11, 3183–3234.
#' * Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
#' An improved GLMNET for l1-regularized logistic regression. 
#' The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
#' 
#' For more details on ISTA, see:
#' 
#' * Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
#' 183–202. https://doi.org/10.1137/080716542
#' * Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
#' Regularized Optimization Problems. Proceedings of the 30th International 
#' Conference on Machine Learning, 28(2)(2), 37–45.
#' * Parikh, N., & Boyd, S. (2013). Proximal Algorithms. Foundations and 
#' Trends in Optimization, 1(3), 123–231.
#'  
#' @param par labeled vector with starting values
#' @param fn R function which takes the parameters AND their labels 
#' as input and returns the fit value (a single value)
#' @param gr R function which takes the parameters AND their labels
#' as input and returns the gradients of the objective function. 
#' If set to NULL, numDeriv will be used to approximate the gradients 
#' @param additionalArguments list with additional arguments passed to fn and gr
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param thetas parameters whose absolute value is above this threshold will be penalized with
#' a constant (theta)
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. 
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
#' @returns Object of class gpRegularized

#' @examples 
#' \donttest{
#' # This example shows how to use the optimizers
#' # for C++ objective functions. We will use
#' # a linear regression as an example. Note that
#' # this is not a useful application of the optimizers
#' # as there are specialized packages for linear regression
#' # (e.g., glmnet)
#' 
#' library(Rcpp)
#' library(lessSEM)
#' 
#' linreg <- '
#' // [[Rcpp::depends(RcppArmadillo)]]
#' #include <RcppArmadillo.h>
#' 
#' // [[Rcpp::export]]
#' double fitfunction(const Rcpp::NumericVector& parameters, Rcpp::List& data){
#'   // extract all required elements:
#'   arma::colvec b = Rcpp::as<arma::colvec>(parameters);
#'   arma::colvec y = Rcpp::as<arma::colvec>(data["y"]); // the dependent variable
#'   arma::mat X = Rcpp::as<arma::mat>(data["X"]); // the design matrix
#'   
#'   // compute the sum of squared errors:
#'     arma::mat sse = arma::trans(y-X*b)*(y-X*b);
#'     
#'     // other packages, such as glmnet, scale the sse with 
#'     // 1/(2*N), where N is the sample size. We will do that here as well
#'     
#'     sse *= 1.0/(2.0 * y.n_elem);
#'     
#'     // note: We must return a double, but the sse is a matrix
#'     // To get a double, just return the single value that is in 
#'     // this matrix:
#'       return(sse(0,0));
#' }
#' 
#' // [[Rcpp::export]]
#' arma::rowvec gradientfunction(const Rcpp::NumericVector& parameters, Rcpp::List& data){
#'   // extract all required elements:
#'   arma::colvec b = Rcpp::as<arma::colvec>(parameters);
#'   arma::colvec y = Rcpp::as<arma::colvec>(data["y"]); // the dependent variable
#'   arma::mat X = Rcpp::as<arma::mat>(data["X"]); // the design matrix
#'   
#'   // note: we want to return our gradients as row-vector; therefore,
#'   // we have to transpose the resulting column-vector:
#'     arma::rowvec gradients = arma::trans(-2.0*X.t() * y + 2.0*X.t()*X*b);
#'     
#'     // other packages, such as glmnet, scale the sse with 
#'     // 1/(2*N), where N is the sample size. We will do that here as well
#'     
#'     gradients *= (.5/y.n_rows);
#'     
#'     return(gradients);
#' }
#' 
#' // https://gallery.rcpp.org/articles/passing-cpp-function-pointers/
#' typedef double (*fitFunPtr)(const Rcpp::NumericVector&, //parameters
#'                 Rcpp::List& //additional elements
#' );
#' typedef Rcpp::XPtr<fitFunPtr> fitFunPtr_t;
#' 
#' typedef arma::rowvec (*gradientFunPtr)(const Rcpp::NumericVector&, //parameters
#'                       Rcpp::List& //additional elements
#' );
#' typedef Rcpp::XPtr<gradientFunPtr> gradientFunPtr_t;
#' 
#' // [[Rcpp::export]]
#' fitFunPtr_t fitfunPtr() {
#'         return(fitFunPtr_t(new fitFunPtr(&fitfunction)));
#' }
#' 
#' // [[Rcpp::export]]
#' gradientFunPtr_t gradfunPtr() {
#'         return(gradientFunPtr_t(new gradientFunPtr(&gradientfunction)));
#' }
#' '
#' 
#' Rcpp::sourceCpp(code = linreg)
#' 
#' ffp <- fitfunPtr()
#' gfp <- gradfunPtr()
#' 
#' N <- 100 # number of persons
#' p <- 10 # number of predictors
#' X <- matrix(rnorm(N*p),	nrow = N, ncol = p) # design matrix
#' b <- c(rep(1,4), 
#'        rep(0,6)) # true regression weights
#' y <- X%*%matrix(b,ncol = 1) + rnorm(N,0,.2)
#' 
#' data <- list("y" = y,
#'              "X" = cbind(1,X))
#' parameters <- rep(0, ncol(data$X))
#' names(parameters) <- paste0("b", 0:(length(parameters)-1))
#' 
#' cL1 <- gpCappedL1Cpp(par = parameters, 
#'                  regularized = paste0("b", 1:(length(b)-1)),
#'                  fn = ffp, 
#'                  gr = gfp, 
#'                  lambdas = seq(0,1,.1), 
#'                  thetas = seq(0.1,1,.1),
#'                  additionalArguments = data)
#' 
#' cL1@parameters
#' }
#' @export
gpCappedL1Cpp <- function(par,
                       fn,
                       gr,
                       additionalArguments,
                       regularized,
                       lambdas,
                       thetas,
                       method = "glmnet",
                       control = lessSEM::controlGlmnet()){
  
  if(!is(fn, "externalptr") | !is(gr, "externalptr")){
    stop("fn and gr must be pointers to C++ functions.")
  }
  if(!is(additionalArguments, "list")){
    stop("additionalArguments must be of class list.")
  }
  
  if(any(thetas <= 0)) stop("Theta must be > 0")
  
  result <- .gpOptimizationInternal(par = par,
                                            fn = fn,
                                            gr = gr,
                                            additionalArguments = additionalArguments,
                                            isCpp = TRUE,
                                            penalty = "cappedL1", 
                                            weights = regularized,
                                            tuningParameters = expand.grid(lambda = lambdas, 
                                                                           theta = thetas,
                                                                           alpha = 1), 
                                            method = method,
                                            control = control
  )
  
  return(result)
  
}

#' gpLspCpp
#' 
#' Implements lsp regularization for general purpose optimization problems with C++ functions.
#' The penalty function is given by:
#' \deqn{p( x_j) = \lambda \log(1 + |x_j|/\theta)}
#' where \eqn{\theta > 0}. 
#' 
#' The interface is inspired by optim, but a bit more restrictive. Users have to supply a vector 
#' with starting values (important: This vector must have labels), a fitting
#' function, and a gradient function. These fitting functions must take an const Rcpp::NumericVector& with parameter
#' values as first argument and an Rcpp::List& as second argument 
#' 
#' lsp regularization:
#' 
#' * Candès, E. J., Wakin, M. B., & Boyd, S. P. (2008). Enhancing Sparsity by 
#' Reweighted l1 Minimization. Journal of Fourier Analysis and Applications, 14(5–6), 
#' 877–905. https://doi.org/10.1007/s00041-008-9045-x
#'  
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classification. Journal of Machine Learning Research, 11, 3183–3234.
#' * Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
#' An improved GLMNET for l1-regularized logistic regression. 
#' The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
#' 
#' For more details on ISTA, see:
#' 
#' * Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
#' 183–202. https://doi.org/10.1137/080716542
#' * Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
#' Regularized Optimization Problems. Proceedings of the 30th International 
#' Conference on Machine Learning, 28(2)(2), 37–45.
#' * Parikh, N., & Boyd, S. (2013). Proximal Algorithms. Foundations and 
#' Trends in Optimization, 1(3), 123–231.
#' 
#' @param par labeled vector with starting values
#' @param fn R function which takes the parameters AND their labels 
#' as input and returns the fit value (a single value)
#' @param gr R function which takes the parameters AND their labels
#' as input and returns the gradients of the objective function. 
#' If set to NULL, numDeriv will be used to approximate the gradients 
#' @param additionalArguments list with additional arguments passed to fn and gr
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param thetas numeric vector: values for the tuning parameter theta
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. 
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
#' @returns Object of class gpRegularized
#' @examples
#' \donttest{
#' # This example shows how to use the optimizers
#' # for C++ objective functions. We will use
#' # a linear regression as an example. Note that
#' # this is not a useful application of the optimizers
#' # as there are specialized packages for linear regression
#' # (e.g., glmnet)
#' 
#' library(Rcpp)
#' library(lessSEM)
#' 
#' linreg <- '
#' // [[Rcpp::depends(RcppArmadillo)]]
#' #include <RcppArmadillo.h>
#' 
#' // [[Rcpp::export]]
#' double fitfunction(const Rcpp::NumericVector& parameters, Rcpp::List& data){
#'   // extract all required elements:
#'   arma::colvec b = Rcpp::as<arma::colvec>(parameters);
#'   arma::colvec y = Rcpp::as<arma::colvec>(data["y"]); // the dependent variable
#'   arma::mat X = Rcpp::as<arma::mat>(data["X"]); // the design matrix
#'   
#'   // compute the sum of squared errors:
#'     arma::mat sse = arma::trans(y-X*b)*(y-X*b);
#'     
#'     // other packages, such as glmnet, scale the sse with 
#'     // 1/(2*N), where N is the sample size. We will do that here as well
#'     
#'     sse *= 1.0/(2.0 * y.n_elem);
#'     
#'     // note: We must return a double, but the sse is a matrix
#'     // To get a double, just return the single value that is in 
#'     // this matrix:
#'       return(sse(0,0));
#' }
#' 
#' // [[Rcpp::export]]
#' arma::rowvec gradientfunction(const Rcpp::NumericVector& parameters, Rcpp::List& data){
#'   // extract all required elements:
#'   arma::colvec b = Rcpp::as<arma::colvec>(parameters);
#'   arma::colvec y = Rcpp::as<arma::colvec>(data["y"]); // the dependent variable
#'   arma::mat X = Rcpp::as<arma::mat>(data["X"]); // the design matrix
#'   
#'   // note: we want to return our gradients as row-vector; therefore,
#'   // we have to transpose the resulting column-vector:
#'     arma::rowvec gradients = arma::trans(-2.0*X.t() * y + 2.0*X.t()*X*b);
#'     
#'     // other packages, such as glmnet, scale the sse with 
#'     // 1/(2*N), where N is the sample size. We will do that here as well
#'     
#'     gradients *= (.5/y.n_rows);
#'     
#'     return(gradients);
#' }
#' 
#' // Dirk Eddelbuettel at
#' // https://gallery.rcpp.org/articles/passing-cpp-function-pointers/
#' typedef double (*fitFunPtr)(const Rcpp::NumericVector&, //parameters
#'                 Rcpp::List& //additional elements
#' );
#' typedef Rcpp::XPtr<fitFunPtr> fitFunPtr_t;
#' 
#' typedef arma::rowvec (*gradientFunPtr)(const Rcpp::NumericVector&, //parameters
#'                       Rcpp::List& //additional elements
#' );
#' typedef Rcpp::XPtr<gradientFunPtr> gradientFunPtr_t;
#' 
#' // [[Rcpp::export]]
#' fitFunPtr_t fitfunPtr() {
#'         return(fitFunPtr_t(new fitFunPtr(&fitfunction)));
#' }
#' 
#' // [[Rcpp::export]]
#' gradientFunPtr_t gradfunPtr() {
#'         return(gradientFunPtr_t(new gradientFunPtr(&gradientfunction)));
#' }
#' '
#' 
#' Rcpp::sourceCpp(code = linreg)
#' 
#' ffp <- fitfunPtr()
#' gfp <- gradfunPtr()
#' 
#' N <- 100 # number of persons
#' p <- 10 # number of predictors
#' X <- matrix(rnorm(N*p),	nrow = N, ncol = p) # design matrix
#' b <- c(rep(1,4), 
#'        rep(0,6)) # true regression weights
#' y <- X%*%matrix(b,ncol = 1) + rnorm(N,0,.2)
#' 
#' data <- list("y" = y,
#'              "X" = cbind(1,X))
#' parameters <- rep(0, ncol(data$X))
#' names(parameters) <- paste0("b", 0:(length(parameters)-1))
#' 
#' l <- gpLspCpp(par = parameters, 
#'                  regularized = paste0("b", 1:(length(b)-1)),
#'                  fn = ffp, 
#'                  gr = gfp, 
#'                  lambdas = seq(0,1,.1), 
#'                  thetas = seq(0.1,1,.1),
#'                  additionalArguments = data)
#' 
#' l@parameters
#' }
#' @export
gpLspCpp <- function(par,
                  fn,
                  gr,
                  additionalArguments,
                  regularized,
                  lambdas,
                  thetas,
                  method = "glmnet",
                  control = lessSEM::controlGlmnet()){
  if(!is(fn, "externalptr") | !is(gr, "externalptr")){
    stop("fn and gr must be pointers to C++ functions.")
  }
  if(!is(additionalArguments, "list")){
    stop("additionalArguments must be of class list.")
  }
  
  if(any(thetas <= 0)) stop("Theta must be > 0")
  
  result <- .gpOptimizationInternal(par = par,
                                            fn = fn,
                                            gr = gr,
                                            additionalArguments = additionalArguments,
                                            isCpp = TRUE,
                                            penalty = "lsp", 
                                            weights = regularized,
                                            tuningParameters = expand.grid(lambda = lambdas, 
                                                                           theta = thetas,
                                                                           alpha = 1), 
                                            method = method,
                                            control = control
  )
  
  return(result)
  
}

#' gpMcpCpp 
#' 
#' Implements mcp regularization for general purpose optimization problems with C++ functions.
#' The penalty function is given by:
#' \ifelse{html}{\deqn{p( x_j) = \begin{cases}
#' \lambda |x_j| - x_j^2/(2\theta) & \text{if } |x_j| \leq \theta\lambda\\
#' \theta\lambda^2/2 & \text{if } |x_j| > \lambda\theta
#' \end{cases}} where \eqn{\theta > 0}.}{
#' Equation Omitted in Pdf Documentation.}
#' 
#' The interface is inspired by optim, but a bit more restrictive. Users have to supply a vector 
#' with starting values (important: This vector _must_ have labels), a fitting
#' function, and a gradient function. These fitting functions _must_ take an const Rcpp::NumericVector& with parameter
#' values as first argument and an Rcpp::List& as second argument
#' 
#' mcp regularization:
#' 
#' * Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax concave penalty. 
#' The Annals of Statistics, 38(2), 894–942. https://doi.org/10.1214/09-AOS729
#'  
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classification. Journal of Machine Learning Research, 11, 3183–3234.
#' * Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
#' An improved GLMNET for l1-regularized logistic regression. 
#' The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
#' 
#' For more details on ISTA, see:
#' 
#' * Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
#' 183–202. https://doi.org/10.1137/080716542
#' * Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
#' Regularized Optimization Problems. Proceedings of the 30th International 
#' Conference on Machine Learning, 28(2)(2), 37–45.
#' * Parikh, N., & Boyd, S. (2013). Proximal Algorithms. Foundations and 
#' Trends in Optimization, 1(3), 123–231.
#' 
#' @param par labeled vector with starting values
#' @param fn R function which takes the parameters AND their labels 
#' as input and returns the fit value (a single value)
#' @param gr R function which takes the parameters AND their labels
#' as input and returns the gradients of the objective function. 
#' If set to NULL, numDeriv will be used to approximate the gradients 
#' @param additionalArguments list with additional arguments passed to fn and gr
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param thetas numeric vector: values for the tuning parameter theta
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. 
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
#' @returns Object of class gpRegularized

#' @examples 
#' \donttest{
#' # This example shows how to use the optimizers
#' # for C++ objective functions. We will use
#' # a linear regression as an example. Note that
#' # this is not a useful application of the optimizers
#' # as there are specialized packages for linear regression
#' # (e.g., glmnet)
#' 
#' library(Rcpp)
#' library(lessSEM)
#' 
#' linreg <- '
#' // [[Rcpp::depends(RcppArmadillo)]]
#' #include <RcppArmadillo.h>
#' 
#' // [[Rcpp::export]]
#' double fitfunction(const Rcpp::NumericVector& parameters, Rcpp::List& data){
#'   // extract all required elements:
#'   arma::colvec b = Rcpp::as<arma::colvec>(parameters);
#'   arma::colvec y = Rcpp::as<arma::colvec>(data["y"]); // the dependent variable
#'   arma::mat X = Rcpp::as<arma::mat>(data["X"]); // the design matrix
#'   
#'   // compute the sum of squared errors:
#'     arma::mat sse = arma::trans(y-X*b)*(y-X*b);
#'     
#'     // other packages, such as glmnet, scale the sse with 
#'     // 1/(2*N), where N is the sample size. We will do that here as well
#'     
#'     sse *= 1.0/(2.0 * y.n_elem);
#'     
#'     // note: We must return a double, but the sse is a matrix
#'     // To get a double, just return the single value that is in 
#'     // this matrix:
#'       return(sse(0,0));
#' }
#' 
#' // [[Rcpp::export]]
#' arma::rowvec gradientfunction(const Rcpp::NumericVector& parameters, Rcpp::List& data){
#'   // extract all required elements:
#'   arma::colvec b = Rcpp::as<arma::colvec>(parameters);
#'   arma::colvec y = Rcpp::as<arma::colvec>(data["y"]); // the dependent variable
#'   arma::mat X = Rcpp::as<arma::mat>(data["X"]); // the design matrix
#'   
#'   // note: we want to return our gradients as row-vector; therefore,
#'   // we have to transpose the resulting column-vector:
#'     arma::rowvec gradients = arma::trans(-2.0*X.t() * y + 2.0*X.t()*X*b);
#'     
#'     // other packages, such as glmnet, scale the sse with 
#'     // 1/(2*N), where N is the sample size. We will do that here as well
#'     
#'     gradients *= (.5/y.n_rows);
#'     
#'     return(gradients);
#' }
#' 
#' // Dirk Eddelbuettel at
#' // https://gallery.rcpp.org/articles/passing-cpp-function-pointers/
#' typedef double (*fitFunPtr)(const Rcpp::NumericVector&, //parameters
#'                 Rcpp::List& //additional elements
#' );
#' typedef Rcpp::XPtr<fitFunPtr> fitFunPtr_t;
#' 
#' typedef arma::rowvec (*gradientFunPtr)(const Rcpp::NumericVector&, //parameters
#'                       Rcpp::List& //additional elements
#' );
#' typedef Rcpp::XPtr<gradientFunPtr> gradientFunPtr_t;
#' 
#' // [[Rcpp::export]]
#' fitFunPtr_t fitfunPtr() {
#'         return(fitFunPtr_t(new fitFunPtr(&fitfunction)));
#' }
#' 
#' // [[Rcpp::export]]
#' gradientFunPtr_t gradfunPtr() {
#'         return(gradientFunPtr_t(new gradientFunPtr(&gradientfunction)));
#' }
#' '
#' 
#' Rcpp::sourceCpp(code = linreg)
#' 
#' ffp <- fitfunPtr()
#' gfp <- gradfunPtr()
#' 
#' N <- 100 # number of persons
#' p <- 10 # number of predictors
#' X <- matrix(rnorm(N*p),	nrow = N, ncol = p) # design matrix
#' b <- c(rep(1,4), 
#'        rep(0,6)) # true regression weights
#' y <- X%*%matrix(b,ncol = 1) + rnorm(N,0,.2)
#' 
#' data <- list("y" = y,
#'              "X" = cbind(1,X))
#' parameters <- rep(0, ncol(data$X))
#' names(parameters) <- paste0("b", 0:(length(parameters)-1))
#' 
#' m <- gpMcpCpp(par = parameters, 
#'                  regularized = paste0("b", 1:(length(b)-1)),
#'                  fn = ffp, 
#'                  gr = gfp, 
#'                  lambdas = seq(0,1,.1), 
#'                  thetas = seq(.1,1,.1),
#'                  additionalArguments = data)
#' 
#' m@parameters
#' }
#' @export
gpMcpCpp <- function(par,
                  fn,
                  gr,
                  additionalArguments,
                  regularized,
                  lambdas,
                  thetas,
                  method = "glmnet",
                  control = lessSEM::controlGlmnet()){
  
  if(!is(fn, "externalptr") | !is(gr, "externalptr")){
    stop("fn and gr must be pointers to C++ functions.")
  }
  if(!is(additionalArguments, "list")){
    stop("additionalArguments must be of class list.")
  }
  
  if(any(thetas <= 0)) stop("Theta must be > 0")
  
  result <- .gpOptimizationInternal(par = par,
                                            fn = fn,
                                            gr = gr,
                                            additionalArguments = additionalArguments,
                                            isCpp = TRUE,
                                            penalty = "mcp", 
                                            weights = regularized,
                                            tuningParameters = expand.grid(lambda = lambdas, 
                                                                           theta = thetas,
                                                                           alpha = 1), 
                                            method = method,
                                            control = control
  )
  
  return(result)
  
}

#' gpScadCpp 
#' 
#' Implements scad regularization for general purpose optimization problems with C++ functions.
#' The penalty function is given by:
#' \ifelse{html}{
#' \deqn{p( x_j) = \begin{cases}
#' \lambda |x_j| & \text{if } |x_j| \leq \theta\\
#' \frac{-x_j^2 + 2\theta\lambda |x_j| - \lambda^2}{2(\theta -1)} & 
#' \text{if } \lambda < |x_j| \leq \lambda\theta \\
#' (\theta + 1) \lambda^2/2 & \text{if } |x_j| \geq \theta\lambda\\
#' \end{cases}}
#' where \eqn{\theta > 2}.}{
#' Equation Omitted in Pdf Documentation.
#' } 
#' 
#' The interface is inspired by optim, but a bit more restrictive. Users have to supply a vector 
#' with starting values (important: This vector _must_ have labels), a fitting
#' function, and a gradient function. These fitting functions _must_ take an const Rcpp::NumericVector& with parameter
#' values as first argument and an Rcpp::List& as second argument
#' 
#' scad regularization:
#' 
#' * Fan, J., & Li, R. (2001). Variable selection via nonconcave penalized 
#' likelihood and its oracle properties. Journal of the American Statistical Association, 
#' 96(456), 1348–1360. https://doi.org/10.1198/016214501753382273
#'  
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classification. Journal of Machine Learning Research, 11, 3183–3234.
#' * Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
#' An improved GLMNET for l1-regularized logistic regression. 
#' The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
#' 
#' For more details on ISTA, see:
#' 
#' * Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
#' 183–202. https://doi.org/10.1137/080716542
#' * Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
#' Regularized Optimization Problems. Proceedings of the 30th International 
#' Conference on Machine Learning, 28(2)(2), 37–45.
#' * Parikh, N., & Boyd, S. (2013). Proximal Algorithms. Foundations and 
#' Trends in Optimization, 1(3), 123–231.
#'  
#' @param par labeled vector with starting values
#' @param fn R function which takes the parameters AND their labels 
#' as input and returns the fit value (a single value)
#' @param gr R function which takes the parameters AND their labels
#' as input and returns the gradients of the objective function. 
#' If set to NULL, numDeriv will be used to approximate the gradients 
#' @param additionalArguments list with additional arguments passed to fn and gr
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param thetas numeric vector: values for the tuning parameter theta
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. 
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
#' @returns Object of class gpRegularized

#' @examples 
#' \donttest{
#' # This example shows how to use the optimizers
#' # for C++ objective functions. We will use
#' # a linear regression as an example. Note that
#' # this is not a useful application of the optimizers
#' # as there are specialized packages for linear regression
#' # (e.g., glmnet)
#' 
#' library(Rcpp)
#' library(lessSEM)
#' 
#' linreg <- '
#' // [[Rcpp::depends(RcppArmadillo)]]
#' #include <RcppArmadillo.h>
#' 
#' // [[Rcpp::export]]
#' double fitfunction(const Rcpp::NumericVector& parameters, Rcpp::List& data){
#'   // extract all required elements:
#'   arma::colvec b = Rcpp::as<arma::colvec>(parameters);
#'   arma::colvec y = Rcpp::as<arma::colvec>(data["y"]); // the dependent variable
#'   arma::mat X = Rcpp::as<arma::mat>(data["X"]); // the design matrix
#'   
#'   // compute the sum of squared errors:
#'     arma::mat sse = arma::trans(y-X*b)*(y-X*b);
#'     
#'     // other packages, such as glmnet, scale the sse with 
#'     // 1/(2*N), where N is the sample size. We will do that here as well
#'     
#'     sse *= 1.0/(2.0 * y.n_elem);
#'     
#'     // note: We must return a double, but the sse is a matrix
#'     // To get a double, just return the single value that is in 
#'     // this matrix:
#'       return(sse(0,0));
#' }
#' 
#' // [[Rcpp::export]]
#' arma::rowvec gradientfunction(const Rcpp::NumericVector& parameters, Rcpp::List& data){
#'   // extract all required elements:
#'   arma::colvec b = Rcpp::as<arma::colvec>(parameters);
#'   arma::colvec y = Rcpp::as<arma::colvec>(data["y"]); // the dependent variable
#'   arma::mat X = Rcpp::as<arma::mat>(data["X"]); // the design matrix
#'   
#'   // note: we want to return our gradients as row-vector; therefore,
#'   // we have to transpose the resulting column-vector:
#'     arma::rowvec gradients = arma::trans(-2.0*X.t() * y + 2.0*X.t()*X*b);
#'     
#'     // other packages, such as glmnet, scale the sse with 
#'     // 1/(2*N), where N is the sample size. We will do that here as well
#'     
#'     gradients *= (.5/y.n_rows);
#'     
#'     return(gradients);
#' }
#' 
#' // Dirk Eddelbuettel at
#' // https://gallery.rcpp.org/articles/passing-cpp-function-pointers/
#' typedef double (*fitFunPtr)(const Rcpp::NumericVector&, //parameters
#'                 Rcpp::List& //additional elements
#' );
#' typedef Rcpp::XPtr<fitFunPtr> fitFunPtr_t;
#' 
#' typedef arma::rowvec (*gradientFunPtr)(const Rcpp::NumericVector&, //parameters
#'                       Rcpp::List& //additional elements
#' );
#' typedef Rcpp::XPtr<gradientFunPtr> gradientFunPtr_t;
#' 
#' // [[Rcpp::export]]
#' fitFunPtr_t fitfunPtr() {
#'         return(fitFunPtr_t(new fitFunPtr(&fitfunction)));
#' }
#' 
#' // [[Rcpp::export]]
#' gradientFunPtr_t gradfunPtr() {
#'         return(gradientFunPtr_t(new gradientFunPtr(&gradientfunction)));
#' }
#' '
#' 
#' Rcpp::sourceCpp(code = linreg)
#' 
#' ffp <- fitfunPtr()
#' gfp <- gradfunPtr()
#' 
#' N <- 100 # number of persons
#' p <- 10 # number of predictors
#' X <- matrix(rnorm(N*p),	nrow = N, ncol = p) # design matrix
#' b <- c(rep(1,4), 
#'        rep(0,6)) # true regression weights
#' y <- X%*%matrix(b,ncol = 1) + rnorm(N,0,.2)
#' 
#' data <- list("y" = y,
#'              "X" = cbind(1,X))
#' parameters <- rep(0, ncol(data$X))
#' names(parameters) <- paste0("b", 0:(length(parameters)-1))
#' 
#' s <- gpScadCpp(par = parameters, 
#'                  regularized = paste0("b", 1:(length(b)-1)),
#'                  fn = ffp, 
#'                  gr = gfp, 
#'                  lambdas = seq(0,1,.1), 
#'                  thetas = seq(2.1,3,.1),
#'                  additionalArguments = data)
#' 
#' s@parameters
#' }
#' @export
gpScadCpp <- function(par,
                   fn,
                   gr,
                   additionalArguments,
                   regularized,
                   lambdas,
                   thetas,
                   method = "glmnet",                      
                   control = lessSEM::controlGlmnet()){
  
  if(!is(fn, "externalptr") | !is(gr, "externalptr")){
    stop("fn and gr must be pointers to C++ functions.")
  }
  if(!is(additionalArguments, "list")){
    stop("additionalArguments must be of class list.")
  }
  
  if(any(thetas <= 2)) stop("Theta must be > 2")
  
  result <- .gpOptimizationInternal(par = par,
                                            fn = fn,
                                            gr = gr,
                                            additionalArguments = additionalArguments,
                                            isCpp = TRUE,
                                            penalty = "scad", 
                                            weights = regularized,
                                            tuningParameters = expand.grid(lambda = lambdas, 
                                                                           theta = thetas,
                                                                           alpha = 1), 
                                            method = method,
                                            control = control
  )
  
  return(result)
}