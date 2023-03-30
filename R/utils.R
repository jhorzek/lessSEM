#' simulateExampleData
#' 
#' simulate data for a simple CFA model
#' @param N number of persons in the data set
#' @param loadings loadings of the latent variable on the manifest observations
#' @param percentMissing percentage of missing data
#' @returns data set for a single-factor CFA.
#' @examples 
#' y <- lessSEM::simulateExampleData()
#' @export
simulateExampleData <- function(N = 100, # sample size
                                loadings = c(rep(1,5), rep(.4,5), rep(0,5)),
                                percentMissing = 0
){
  
  f <- matrix(stats::rnorm(N, 0, 1), ncol = 1) # latent factor
  L <- matrix(loadings, 
              nrow = 1) # loadings
  # covariances
  covs <- diag(max(L^2)+.2, length(loadings))
  
  y <- matrix(NA, nrow = N, ncol = ncol(L))
  
  for(i in 1:N){
    y[i,] <- L*f[i,] +  mvtnorm::rmvnorm(1, sigma = covs)
  }
  
  yNames <- paste0("y", 1:ncol(y))
  colnames(y) <- yNames
  
  if(percentMissing != 0){
    
    allPoints <- expand.grid(row = 1:nrow(y),
                             col = 1:ncol(y))
    missing <- sample(size = round(nrow(allPoints)*percentMissing/100),
                      x = 1:nrow(allPoints))
    for(i in missing){
      y[allPoints$row[i], allPoints$col[i]] <- NA
    }
  }
  
  return(y)
  
}

#' .noDotDotDot
#' 
#' remplaces the dot dot dot part of the fitting and gradient fuction
#' @param fn fit or gradient function. IMPORTANT: THE FIRST ARGUMENT TO
#' THE FUNCTION MUST BE THE PARAMETER VECTOR
#' @param fnName name of the function fn
#' @param ... additional arguments
#' @return list with (1) new function which wraps fn and (2) list with arguments passed to fn
#' @keywords internal
.noDotDotDot <- function(fn, fnName, ...){
  
  # see Barranka, https://stackoverflow.com/questions/26164078/r-define-a-function-from-character-string
  
  dotdotdot <- list(...)
  fnUser <- fn
  
  argsAre <- formalArgs(fnUser)
  
  if(any(!names(dotdotdot) %in% argsAre)){
    warning(paste0(
      "You passed the following argument(s) using ...: ", 
      paste0(names(dotdotdot)[!names(dotdotdot) %in% argsAre], collapse = ", "),
      ". The fitting function or gradient function seems to not use all of these arguments."
    )
    )
  }
  
  namesOfPar <- argsAre[1] # the first argument is the parameter vector
  
  if(length(argsAre) == 1) {
    additionalArguments <- list(fnUser)
    names(additionalArguments) <- paste0(fnName, "User")
    eval(
      parse(
        text = paste('fn <- function(par, parameterLabels, additionalArguments) { 
                          names(par) <- parameterLabels
                          return(additionalArguments$', fnName, 'User(', namesOfPar, ' = par))}', 
                     sep='')
      )
    )
    
    ret <- list(fn,
                additionalArguments
    )
    names(ret) <- c(paste0(fnName, "User"), "additionalArguments")
    
    return(ret)
    
    
  }
  
  
  namesOfAdditional <- argsAre[2:length(argsAre)]
  
  fnArgs <- c(
    paste0(namesOfPar[1]," = par"), 
    paste0(namesOfAdditional, " = ", "additionalArguments$", namesOfAdditional)
  )
  
  
  body <- paste0(
    "additionalArguments$", fnName, "User(", 
    paste0(fnArgs, collapse = ", "), 
    ")")
  
  additionalArguments <- c(dotdotdot, fnUser)
  names(additionalArguments) <- c(names(dotdotdot), paste0(fnName, "User"))
  
  rm(fn)
  eval(parse(text = paste('fn <- function(par, parameterLabels, additionalArguments) { 
                          names(par) <- parameterLabels
                          return(' , body , ')}', sep='')))
  
  ret <- list(fn,
              additionalArguments
  )
  names(ret) <- c(paste0(fnName, "User"), "additionalArguments")
  
  return(ret)
}

#' makePtrs
#' 
#' This function helps you create the pointers necessary to use the Cpp interface
#' 
#' @param fitFunName name of your C++ fit function (IMPORTANT: This must be the name
#' used in C++)
#' @param gradFunName name of your C++ gradient function (IMPORTANT: This must be the name
#' used in C++)
#' @returns a string which can be copied in the C++ function to create the pointers.
#' @examples 
#' # see vignette("General-Purpose-Optimization", package = "lessSEM") for an example
#' @export
makePtrs <- function(fitFunName, gradFunName){
  
  makePtrsSyntax <- paste0(
    '
// INSTRUCTIONS: ADD THE FOLLOWING LINES TO YOUR C++ FUNCTIONS

// IF RCPPARMADILLO IS NOT IMPORTED YET, UNCOMMENT THE FOLLOWING TWO LINES
// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>

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
fitFunPtr_t ', fitFunName,'Ptr() {
        return(fitFunPtr_t(new fitFunPtr(&',fitFunName,')));
}

// [[Rcpp::export]]
gradientFunPtr_t ', gradFunName,'Ptr() {
        return(gradientFunPtr_t(new gradientFunPtr(&',gradFunName,')));
}
'
  )
  return(makePtrsSyntax)
}