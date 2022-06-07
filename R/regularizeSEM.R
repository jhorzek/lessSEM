#' regularizeSEM
#' 
#' This function provides optimization for regularized structural equation models with ridge,
#' lasso, adaptive lasso, or elastic net penalty. In case of lasso and adaptive lasso, the lambda values
#' can be selected automatically. The returned object is an S4 class; its elements can be accessed
#' with the "@" operator (see examples).
#' 
#' # References
#' 
#' Jacobucci, R., Grimm, K. J., & McArdle, J. J. (2016). Regularized Structural Equation Modeling. 
#' Structural Equation Modeling: A Multidisciplinary Journal, 23(4), 555–566. https://doi.org/10.1080/10705511.2016.1154793
#' 
#' Huang, P.-H., Chen, H., & Weng, L.-J. (2017). A Penalized Likelihood Method for Structural Equation Modeling. 
#' Psychometrika, 82(2), 329–354. https://doi.org/10.1007/s11336-017-9566-9
#' 
#' @param lavaanModel model of class lavaan 
#' @param regularizedParameterLabels labels of regularized parameters
#' @param penalty which penalty should be used? Available are "ridge", "lasso", "adaptiveLasso", and "elasticNet"
#' @param lambdas vector with lambda values. Higher values = higher penalty
#' @param nLambdas if penalty == "lasso" or penalty == "adaptiveLasso", one can specify the number of lambda values to use. In this case, set lambdas = NULL.
#' @param alphas 0<alpha<1. only required for elastic net. Controls the weight of ridge and lasso terms. alpha = 1 is lasso, alpha = 0 ridge
#' @param adaptiveLassoWeights vector with weights for adaptive LASSO. Set to NULL if not using adaptive LASSO. Default is inverse of absolute unregularized parameter estimates
#' @param raw controls if the internal transformations of aCV4SEM is used. Recommended to set to TRUE!
#' @param control option to set parameters of the GLMNET optimizer. See ?controlGLMNET
#' @examples  
#' library(aCV4SEM)
#' 
#' # Identical to regsem, aCV4SEM builds on the lavaan
#' # package for model specification. The first step
#' # therefore is to implement the model in lavaan.
#' 
#' dataset <- simulateExampleData()
#' 
#' lavaanSyntax <- "
#' f =~ l1*y1 + l2*y2 + l3*y3 + l4*y4 + l5*y5 + 
#'      l6*y6 + l7*y7 + l8*y8 + l9*y9 + l10*y10 + 
#'      l11*y11 + l12*y12 + l13*y13 + l14*y14 + l15*y15
#' f ~~ 1*f
#' "
#' 
#' lavaanModel <- lavaan::sem(lavaanSyntax,
#'                            data = dataset,
#'                            meanstructure = TRUE,
#'                            std.lv = TRUE)
#' 
#' # Optional: Plot the model
#' # semPlot::semPaths(lavaanModel, 
#' #                   what = "est",
#' #                   fade = FALSE)
#' 
#' regsem <- regularizeSEM(# pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized paramters:
#'   regularizedParameterLabels = paste0("l", 6:15),
#'   # which penalty should be used?
#'   penalty = "lasso",
#'   # in case of lasso and adaptive lasso, we can specify the number of lambda
#'   # values to use. aCV4SEM will automatically find lambda_max and fit
#'   # models for nLambda values between 0 and lambda_max. For the other
#'   # penalty functions, lambdas must be specified explicitly
#'   nLambdas = 5)
#' 
#' # use the plot-function to plot the regularized parameters:
#' plot(regsem)
#' 
#' # elements of regsem can be accessed with the @ operator:
#' regsem@parameters[1,]
#' 
#' # AIC and BIC:
#' AIC(regsem)
#' BIC(regsem)
#' 
#' # The best parameters can also be extracted with:
#' coef(regsem, criterion = "AIC")
#' coef(regsem, criterion = "BIC")
#' 
#' ## The fitted model can then be used as basis for an approximate cross-validation
#' # (see ?aCV4SEM::aCV4regularizedSEM) or approximate influence functions
#' # (see ?aCV4SEM::aI4regularizedSEM)
#' @export
regularizeSEM <- function(lavaanModel, 
                          regularizedParameterLabels,
                          penalty, 
                          lambdas = NULL,
                          nLambdas = NULL,
                          alphas = NULL, 
                          adaptiveLassoWeights = NULL,
                          raw = TRUE,
                          control = aCV4SEM::controlGLMNET()){
  
  inputArguments <- as.list(environment())
  
  if(!is(control, "controlGLMNET")){
    stop("control must be of class controlGLMNET. These objects can be generated with the controlGLMNET function. See ?aCV4SEM::controlGLMNET")
  }
  
  if(!penalty %in% c("lasso", "ridge", "adaptiveLasso", "elasticNet")) stop("Currently supported penalty functions are: lasso, ridge, adaptiveLasso, and elasticNet")
  
  if(is.null(lambdas) && is.null(nLambdas)) stop("Function requires either lambdas or nLambdas to be specified.")
  if(!is.null(lambdas) && !is.null(nLambdas)) stop("Please specify either lambdas or nLambdas, but not both.")
  if(is.null(lambdas) && any(!is.null(alphas), penalty %in% c("ridge", "elasticNet"))) stop("nLambdas is currently not supported for ridge and elasticNet penalty. Please specify lambdas (e.g., lambdas = seq(0,1,.1)).")
  
  if(!is(lavaanModel, "lavaan")){
    stop("lavaanModel must be of class lavaan")
  }
  
  if(lavaanModel@Options$estimator != "ML") stop("lavaanModel must be fit with ml estimator.")
  
  rawData <- try(lavaan::lavInspect(lavaanModel, "data"))
  if(is(rawData, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  
  ### initialize model ####
  startingValues <- control$startingValues
  if(any(startingValues == "est")){
    SEM <- aCV4SEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                   transformVariances = TRUE,
                                   whichPars = "est",
                                   addMeans = control$addMeans, 
                                   activeSet = control$activeSet)
  }else if(any(startingValues == "start")){
    SEM <- aCV4SEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                   transformVariances = TRUE,
                                   whichPars = "start",
                                   addMeans = control$addMeans, 
                                   activeSet = control$activeSet)
  }else if(is.numeric(startingValues)){
    
    if(!all(names(startingValues) %in% names(aCV4SEM::getLavaanParameters(lavaanModel)))) stop("Parameter names of startingValues do not match those of the lavaan object. See aCV4SEM::getLavaanParameters(lavaanModel).")
    SEM <- aCV4SEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                   transformVariances = TRUE,
                                   whichPars = "start", 
                                   fit = FALSE,
                                   addMeans = control$addMeans, 
                                   activeSet = control$activeSet)
    SEM <- aCV4SEM:::setParameters(SEM = SEM, labels = names(startingValues), value = startingValues, raw = FALSE)
    SEM <- try(aCV4SEM:::fit(SEM))
    if(is(SEM, "try-error") || !is.finite(SEM$m2LL)) stop("Infeasible starting values.")
    
  }else{
    stop("Invalid startingValues passed to GLMNET. See ?controlGLMNET for more information.")
  }
  
  # get parameters
  parameters <- aCV4SEM:::getParameters(SEM, raw = raw)
  aCV4SEM:::checkRegularizedParameters(parameters = parameters, 
                                       regularizedParameterLabels = regularizedParameterLabels,
                                       SEM$getParameters(), 
                                       raw = raw)
  
  initialHessian <- control$initialHessian
  if(is.matrix(initialHessian) && nrow(initialHessian) == length(parameters) && ncol(initialHessian) == length(parameters)){
    
    if(!all(rownames(initialHessian) %in% names(aCV4SEM::getLavaanParameters(lavaanModel))) ||
       !all(colnames(initialHessian) %in% names(aCV4SEM::getLavaanParameters(lavaanModel)))
    ) stop("initialHessian must have the parameter names as rownames and colnames. See aCV4SEM::getLavaanParameters(lavaanModel).")
    
  }else if(any(initialHessian == "compute")){
    
    initialHessian <- aCV4SEM:::getHessian(SEM = SEM, raw = TRUE)
    
  }else if(any(initialHessian == "scoreBased")){
    
    scores <- aCV4SEM:::getScores(SEM = SEM, raw = TRUE)
    FisherInformation <- matrix(0, nrow = ncol(scores), ncol = ncol(scores))
    rownames(FisherInformation) <- colnames(FisherInformation) <- colnames(scores)
    for(score in 1:nrow(scores)){
      FisherInformation <- FisherInformation + t(-.5*scores[score,, drop = FALSE]) %*%(-.5*scores[score,, drop = FALSE]) # we are using the -2 log-Likelihood
    }
    
    initialHessian <- -2*(-FisherInformation) # negative log-likelihood
    # make symmetric; just in case...
    initialHessian <- .5*(initialHessian + t(initialHessian))
    while(any(eigen(initialHessian, only.values = TRUE)$values < 0)){
      diag(initialHessian) <- diag(initialHessian) + 1e-4
    }
    
  }else if(length(initialHessian) == 1 && is.numeric(initialHessian)){
    initialHessian <- diag(initialHessian,length(parameters))
    rownames(initialHessian) <- names(parameters)
    colnames(initialHessian) <- names(parameters)
  }else{
    stop("Invalid initialHessian passed to GLMNET. See ?controlGLMNET for more information.")
  }
  
  initialHessian4Optimizer <- initialHessian
  
  if(penalty == "adaptiveLasso" && is.null(adaptiveLassoWeights)){
    message("adaptiveLasso selected, but no adaptiveLassoWeights provided. Using the default abs(parameters)^(-1).")
    adaptiveLassoWeights <- abs(parameters)^(-1)
  }else if(is.null(adaptiveLassoWeights)){
    adaptiveLassoWeights <- rep(1, length(parameters)) # we are using the 
    names(adaptiveLassoWeights) <- names(parameters)
  }else{
    message("Using user specified adaptive lasso weights.")
    
    if(any(alphas != 1)) warning("Combining ridge and elastic net with adaptive lasso weights is unexplored territory.")
  }
  inputArguments$adaptiveLassoWeights <- adaptiveLassoWeights
  
  if(!is.null(alphas) && (!all(alphas %in% c(0,1))) && penalty != "elasticNet") {
    stop("non-null alpha parameter only valid for elasticNet")
  }
  if(is.null(alphas)){
    if(penalty == "elasticNet") stop("elasticNet requires specification of alpha parameter.")
  }
  if(penalty %in% c("lasso", "adaptiveLasso")){ 
    alphas <- 1
  }else if(penalty == "ridge"){
    alphas <- 0
  }
  
  if(!is.null(nLambdas)){
    message("Automatically selecting the maximal lambda value. Note: This may not work properly if a model with all regularized parameters set to zero is not identified.")
    maxLambda <- aCV4SEM:::getMaxLambda(SEM = SEM, 
                                        regularizedParameterLabels = regularizedParameterLabels, 
                                        adaptiveLassoWeights = adaptiveLassoWeights,
                                        initialHessian4Optimizer = initialHessian4Optimizer, 
                                        control = control,
                                        N = nrow(rawData))
    lambdas <- seq(0, maxLambda, length.out = nLambdas)
  }
  
  tuningGrid <- expand.grid("lambda" = lambdas, "alpha" = alphas)
  
  fits <- data.frame(
    tuningGrid,
    "m2LL" = NA,
    "regM2LL"= NA,
    "nonZeroParameters" = NA,
    "convergence" = NA
  )
  
  parameterEstimates <- as.data.frame(matrix(NA,nrow = nrow(tuningGrid), ncol = length(parameters)))
  colnames(parameterEstimates) <- names(parameters)
  parameterEstimates <- cbind(
    tuningGrid,
    parameterEstimates
  )
  
  if(control$saveHessian){
    Hessians <- list(
      "lambda" = tuningGrid$lambda,
      "alpha" = tuningGrid$alpha,
      "Hessian" = lapply(1:nrow(tuningGrid), 
                         matrix, 
                         data= NA, 
                         nrow=nrow(initialHessian), 
                         ncol=ncol(initialHessian))
    )
  }else{
    Hessians <- list(NULL)
  }
  
  if(control$verbose == 0){
    progressbar = txtProgressBar(min = 0, 
                                 max = nrow(tuningGrid), 
                                 initial = 0, 
                                 style = 3)
  }
  
  for(it in 1:nrow(tuningGrid)){
    if(control$verbose == 0){
      setTxtProgressBar(progressbar,it)
    }else{
      cat(paste0("\nIteration [", it, "/", nrow(tuningGrid),"]\n"))
    }
    
    lambda <- tuningGrid$lambda[it]
    alpha <- tuningGrid$alpha[it]
    
    result <- try(aCV4SEM:::GLMNET(SEM = SEM, 
                                   regularizedParameterLabels = regularizedParameterLabels, 
                                   lambda = lambda, 
                                   alpha = alpha,
                                   adaptiveLassoWeights = adaptiveLassoWeights,
                                   initialHessian = initialHessian4Optimizer,
                                   stepSize = control$stepSize,
                                   sig = control$sig,
                                   gam = control$gam,
                                   maxIterOut = control$maxIterOut,
                                   maxIterIn = control$maxIterIn,
                                   maxIterLine = control$maxIterLine,
                                   epsOut = control$epsOut,
                                   epsIn = control$epsIn,
                                   convergenceCriterion = control$convergenceCriterion,
                                   verbose = control$verbose))
    if(is(result, "try-error")) next
    
    fits$m2LL[it] <- result$m2LL
    fits$regM2LL[it] <- result$regM2LL
    fits$nonZeroParameters[it] <- result$nonZeroParameters
    fits$convergence[it] <- result$convergence
    newParameters <- aCV4SEM:::getParameters(result$SEM, raw = FALSE)
    
    parameterEstimates[it, names(parameters)] <- newParameters[names(parameters)]
    
    if(control$saveHessian) Hessians$Hessian[[it]] <- result$Hessian
    
    # set initial values for next iteration
    SEM <- aCV4SEM:::setParameters(SEM = SEM, labels = names(newParameters), value = newParameters, raw = FALSE)
    SEM <- try(aCV4SEM:::fit(SEM))
    if(is(SEM, "try-Error")){
      # reset
      warning("Fit for lambda = ",lambda, "alpha = ", alpha, " resulted in Error!")
      SEM <- aCV4SEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                     transformVariances = TRUE,
                                     whichPars = "start",
                                     addMeans = control$addMeans)
      initialHessian4Optimizer <- NULL
    }else{
      initialHessian4Optimizer <- result$Hessian
    }
    
  }
  
  internalOptimization <- list(
    "HessiansOfDifferentiablePart" = Hessians
  )
  
  results <- new("regularizedSEM",
                 parameters = parameterEstimates,
                 fits = fits,
                 parameterLabels = names(parameters),
                 regularized = regularizedParameterLabels,
                 internalOptimization = internalOptimization,
                 inputArguments = inputArguments)
  
  return(results)
  
}

#' optimizeSEM
#' 
#' optimize a SEM 
#' 
#' @param SEM model of class Rcpp_SEMCpp. 
#' @param raw controls if the internal transformations of aCV4SEM is used.
#' @param method optimizer method. See ?optim
#' @param control control optimizer. See ?optim
#' @param hessian Should the optimizer return the hessian?
optimizeSEM <- function(SEM, 
                        raw = TRUE, 
                        method = "BFGS", 
                        control = list(), 
                        hessian = FALSE){
  
  
  parameters <- aCV4SEM:::getParameters(SEM, raw = raw)
  
  fitFun <- function(par, SEM, raw){
    m2LL <- aCV4SEM:::fitFunction(par = par, SEM = SEM, raw = raw)
    if(m2LL == 9999999999999) return(m2LL)
    return(m2LL)
  }
  
  gradFun <- function(par, SEM, raw){
    gradients <- aCV4SEM:::derivativeFunction(par = par, SEM = SEM, raw = raw)
    if(all(gradients == 9999999999999)) return(gradients)
    return(gradients)
  }
  
  
  optimized <- optim(par = parameters, 
                     fn = fitFun, 
                     gr = gradFun, 
                     SEM = SEM,
                     raw = raw,
                     method = method,
                     control = control,
                     hessian = hessian)
  
  SEM <- aCV4SEM:::setParameters(SEM, names(optimized$par), optimized$par, raw = raw)
  SEM <- try(aCV4SEM:::fit(SEM), silent = TRUE)
  
  return(list("SEM" = SEM,
              "optimizer" = optimized))
}

#' checkRegularizedParameters
#' 
#' internal function to check if the regularized parameters can be found in the model
#' 
#' @param parameters labeled vector with parameter values
#' @param regularizedParameterLabels labels of regularized parameters
#' @param parameterTable table with all parameters
#' @param raw controls if the internal transformations of aCV4SEM is used.
checkRegularizedParameters <- function(parameters, regularizedParameterLabels, parameterTable, raw){
  
  
  # check regularizedParameterLabels
  if(any(!regularizedParameterLabels %in% names(parameters))) stop(paste0("Could not find parameter ", 
                                                                          paste0(regularizedParameterLabels[!regularizedParameterLabels %in% names(parameters)], collapse = ", "))
  )
  
  # check if variances are regularized and the parameters are raw
  if(raw){
    for(regularizedParmeter in regularizedParameterLabels){
      if(all(parameterTable$location[parameterTable$label == regularizedParmeter] != "Smatrix")) next
      if(any(parameterTable$row[parameterTable$label == regularizedParmeter] == parameterTable$col[parameterTable$label == regularizedParmeter])){
        warning("Regularizing raw_value of variances (raw = TRUE). The actual variances are given by var = exp(raw_value).")
      }
    }
  }
}

#' getMaxLambda
#' 
#' generates a the first lambda which sets all regularized parameters to zero
#' @param SEM model of class Rcpp_SEMCpp
#' @param regularizedParameterLabels labels of regularized parameters
#' @param adaptiveLassoWeights vector with weights for adaptive LASSO
#' @param initialHessian4Optimizer initial Hessian used by the optimizer
#' @param control additional arguments passed to GLMNET. See ?controlGLMNET
#' @param N sample size
getMaxLambda <- function(SEM, regularizedParameterLabels, adaptiveLassoWeights, initialHessian4Optimizer, control, N){
  initialParameters <- aCV4SEM:::getParameters(SEM = SEM, raw = TRUE)
  lambda <- .Machine$double.xmax^(.05)
  result <- aCV4SEM:::GLMNET(SEM = SEM, 
                             regularizedParameterLabels = regularizedParameterLabels, 
                             lambda = lambda, 
                             alpha = 1,
                             adaptiveLassoWeights = adaptiveLassoWeights,
                             initialHessian = initialHessian4Optimizer,
                             stepSize = control$stepSize,
                             sig = control$sig,
                             gam = control$gam,
                             maxIterOut = control$maxIterOut,
                             maxIterIn = control$maxIterIn,
                             maxIterLine = control$maxIterLine,
                             epsOut = control$epsOut,
                             epsIn = control$epsIn,
                             convergenceCriterion = control$convergenceCriterion,
                             verbose = control$verbose)
  sparseParameters <- result$parameters
  SEM <- aCV4SEM:::setParameters(SEM = SEM, labels = names(sparseParameters), values = sparseParameters, raw = TRUE)
  SEM <- aCV4SEM:::fit(SEM = SEM)
  gradients <- aCV4SEM:::getGradients(SEM = SEM, raw = TRUE)
  
  # define maxLambda as the maximal gradient of the regularized parameters
  maxLambda <- max(abs(gradients[regularizedParameterLabels]) * adaptiveLassoWeights[regularizedParameterLabels]^(-1))
  
  SEM <- aCV4SEM:::setParameters(SEM = SEM, labels = names(initialParameters), values = initialParameters, raw = TRUE)
  SEM <- aCV4SEM:::fit(SEM = SEM)
  
  return((1/N)*(maxLambda+.1*maxLambda)) # adding some wiggle room as well
}