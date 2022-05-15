#' regularizeSEMWithCustomPenalty
#' 
#' optimize a SEM with custom penalty functions.
#' 
#' @param lavaanModel model of class lavaan 
#' @param individualPenaltyFunction penalty function which takes the current parameter values as first argument, the tuning parameters as second, and the penaltyFunctionArguments as third argument and 
#' returns a single value - the value of the penalty function for a single person. If the true penalty function is non-differentiable (e.g., lasso) a smooth
#' approximation of this function should be provided.
#' @param individualPenaltyFunctionGradient gradients of the penalty function. Function should take the current parameter values as first argument, the tuning parameters as second, and the penaltyFunctionArguments as third argument and 
#' return a vector of the same length as parameters. If the true penalty function is non-differentiable (e.g., lasso) a smooth
#' approximation of this function should be provided. If set to NULL, numDeriv will be used.
#' @param individualPenaltyFunctionHessian Hessian of the penalty function. Function should take the current parameter values as first argument, the tuning parameters as second, and the penaltyFunctionArguments as third argument and 
#' return a matrix with (length as parameters)^2 number of elements. If the true penalty function is non-differentiable (e.g., lasso) a smooth
#' approximation of this function should be provided. If set to NULL, numDeriv will be used.
#' @param tuningParameters data.frame with tuning parameter values. Important: The function will iterate over the rows of these tuning parameters and pass them to your penalty function
#' @param penaltyFunctionArguments arguments passed to individualPenaltyFunction, individualPenaltyFunctionGradient, and individualPenaltyFunctionHessian
#' @param control option to set parameters of the optimizer
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
#' ## Defining a custom penalty function is a bit more complicated than
#' # using the default ones provided in regularizeSEM (see ?aCV4SEM::regularizeSEM).
#' # We start with the definition of the penalty function. Make sure that the derivatives
#' # of this function are well defined (i.e., the function is smooth)!
#' # We will use the lasso penalty as an example.
#' 
#' # The penalty function MUST accept three arguments; how you use them is up to you.
#' 
#' # The first argument are the parameters. aCV4SEM will pass a vector with the current
#' # parameter values and their names to your function. The parameter labels will be
#' # exactly the same as those used by lavaan! So you can check them before with:
#' 
#' aCV4SEM::getLavaanParameters(lavaanModel)
#' 
#' # The vector passed to your function will look like the vector above.
#' 
#' # The second argument is called tuningParameters MUST be a data.frame with the
#' # tuning-parameters. aCV4SEM will automatically iterate over the rows of this object.
#' # In case of LASSO regularization there is only one tuning parameter: lambda.
#' # Therefore, we specify the tuningParameters object as:
#' 
#' tuningParameters <- data.frame(lambda = seq(0,1,.1)) # we will test 11 lambdas here
#' print(tuningParameters)
#' 
#' # The third argument is called penaltyFunctionArguments and we can pass anything
#' # we want here. For the lasso penalty, we need two additional things:
#' # 1) We need the smoothing-parameter epsilon, which makes sure that our
#' # penalty is differentiable
#' # 2) We need the labels of the regularized parameters, so that we can only
#' # penalize those while all others remain unpenalized.
#' 
#' penaltyFunctionArguments <- list(
#'   eps = 1e-10,
#'   regularizedParameterLabels = paste0("l", 6:15)
#' )
#' 
#' # Now, it is time to specify our custom penalty function:
#' 
#' smoothLASSO <- function(
    #'   # here are our three arguments:
#'   parameters,
#'   tuningParameters,
#'   penaltyFunctionArguments
#' ){
#'   # to make it easier to see what is going on:
#'   lambda <- tuningParameters$lambda # tuningParameters will be ONE ROW OF the
#'   # tuningParameters object we created before -> it will contain just
#'   # one lambda in this case!
#'   eps <- penaltyFunctionArguments$eps
#'   regularizedParameterLabels <- penaltyFunctionArguments$regularizedParameterLabels
#' 
#'   regularizedParameters <- parameters[regularizedParameterLabels]
#' 
#'   # now, let's define our penalty function:
#'   penaltyLasso <- lambda*sum(sqrt(regularizedParameters^2 + eps))
#' 
#'   return(penaltyLasso)
#' }
#' # Important: This penalty function is assumed to be for a single individual only.
#' # aCV4SEM will multiply it with sample size N to get the penalty value of the
#' # full sample!
#' 
#' #### Now we are ready to optimize! ####
#' regsemApprox <- regularizeSEMWithCustomPenalty(lavaanModel = lavaanModel,
#'                                                individualPenaltyFunction = smoothLASSO,
#'                                                tuningParameters = tuningParameters,
#'                                                penaltyFunctionArguments = penaltyFunctionArguments)
#' 
#' # let's compare the results to an exact optimization:
#' 
#' regsemExact <- regularizeSEM(
#'   lavaanModel = lavaanModel,
#'   regularizedParameterLabels = paste0("l", 6:15),
#'   penalty = "lasso",
#'   lambdas = tuningParameters$lambda)
#' 
#' head(regsemExact@parameters[,regsemExact@parameterLabels] -
#'        regsemApprox@parameters[,regsemExact@parameterLabels])
#' # Note that the parameter estimates are basically identical.
#' 
#' ## To select a model, we used the approximate cross-validation function:
#' aCV <- aCV4regularizedSEMWithCustomPenalty(regularizedSEMWithCustomPenalty = regsemApprox,
#'                                            k = nrow(dataset))
#' plot(aCV)
#' # To extract the best parameter estimates:
#' coef(aCV)
#' 
#' #### More advanced features ####
#' # You may have seen that regularizeSEMWithCustomPenalty returned warnings that :
#' # 1. The gradients of the penalty function are approximated.
#' # 2. The Hessian of the penalty function is approximated.
#' # Especially (1) can slow down the optimization. We can take care of this by
#' # defining a custom gradient function. This function MUST take the same three arguments
#' # as the penalty function
#' 
#' smoothLASSOGradient <- function(parameters,
#'                                 tuningParameters,
#'                                 penaltyFunctionArguments){
#'   gradients <- rep(0, length(parameters))
#'   names(gradients) <- names(parameters)
#' 
#'   # to make it easier to see what is going on:
#'   lambda <- tuningParameters$lambda
#'   eps <- penaltyFunctionArguments$eps
#'   regularizedParameterLabels <- penaltyFunctionArguments$regularizedParameterLabels
#' 
#'   regularizedParameters <- parameters[regularizedParameterLabels]
#' 
#'   # now, lets compute the gradients of the penalty function:
#'   gradients[regularizedParameterLabels] <- lambda*
#'     regularizedParameters*
#'     (1/sqrt((regularizedParameters)^2 + eps))
#' 
#'   return(gradients)
#' }
#' 
#' # Now, let's optimize again:
#' regsemApprox2 <- regularizeSEMWithCustomPenalty(lavaanModel = lavaanModel,
#'                                                individualPenaltyFunction = smoothLASSO,
#'                                                individualPenaltyFunctionGradient = smoothLASSOGradient,
#'                                                tuningParameters = tuningParameters,
#'                                                penaltyFunctionArguments = penaltyFunctionArguments)
#' 
#' # let's compare the results to the first approximate optimization:
#' 
#' head(regsemApprox2@parameters[,regsemApprox2@parameterLabels] -
#'        regsemApprox@parameters[,regsemApprox2@parameterLabels])
#' # In this case, the additional efforts of defining a custom gradient function
#' # we quite unnecessary. However, there are other cases where this may play a
#' # role.
#' 
#' # In our experience, defining an additional Hessian function is rarely necessary,
#' # but for completeness sake:
#' smoothLASSOHessian <- function(parameters,
#'                                tuningParameters,
#'                                penaltyFunctionArguments){
#'   hessian <- matrix(0,
#'                     length(parameters),
#'                     length(parameters),
#'                     dimnames = list(names(parameters),
#'                                     names(parameters)))
#' 
#'   # to make it easier to see what is going on:
#'   lambda <- tuningParameters$lambda
#'   eps <- penaltyFunctionArguments$eps
#'   regularizedParameterLabels <- penaltyFunctionArguments$regularizedParameterLabels
#' 
#'   regularizedParameters <- parameters[regularizedParameterLabels]
#' 
#'   diag(hessian[regularizedParameterLabels,
#'                regularizedParameterLabels]) <- -lambda*
#'     regularizedParameters^2*
#'     (regularizedParameters^2 + eps)^(-3/2) +
#'     lambda*
#'     (1/sqrt(regularizedParameters^2 + eps))
#' 
#'   return(hessian)
#' }
#' 
#' # Now, let's optimize again:
#' regsemApprox3 <- regularizeSEMWithCustomPenalty(lavaanModel = lavaanModel,
#'                                                 individualPenaltyFunction = smoothLASSO,
#'                                                 individualPenaltyFunctionGradient = smoothLASSOGradient,
#'                                                 individualPenaltyFunctionHessian = smoothLASSOHessian,
#'                                                 tuningParameters = tuningParameters,
#'                                                 penaltyFunctionArguments = penaltyFunctionArguments)
#' 
#' # let's compare the results to the second approximate optimization:
#' 
#' head(regsemApprox3@parameters[,regsemApprox3@parameterLabels] -
#'        regsemApprox2@parameters[,regsemApprox3@parameterLabels])
#' 
#' # Again, the parameter estimates are basically identical.
#' @export
regularizeSEMWithCustomPenalty <- function(lavaanModel, 
                                           individualPenaltyFunction,
                                           individualPenaltyFunctionGradient = NULL,
                                           individualPenaltyFunctionHessian = NULL,
                                           tuningParameters,
                                           penaltyFunctionArguments,
                                           control = aCV4SEM::controlQuasiNewtonBFGS()){
  
  if(is.null(individualPenaltyFunctionGradient)){
    message("Using numDeriv to approximate the gradient of the individualPenaltyFunction. This may result in slow optimization. We highly recommend that you pass an analytic solution.")
    individualPenaltyFunctionGradient <- aCV4SEM::genericGradientApproximiation
    penaltyFunctionArguments <- c(penaltyFunctionArguments, individualPenaltyFunction = individualPenaltyFunction)
  }
  if(is.null(individualPenaltyFunctionHessian)){
    message("Using numDeriv to approximate the Hessian of the individualPenaltyFunction. This may result in slow optimization. We highly recommend that you pass an analytic solution.")
    individualPenaltyFunctionHessian <- aCV4SEM::genericHessianApproximiation
    penaltyFunctionArguments <- c(penaltyFunctionArguments, individualPenaltyFunction = individualPenaltyFunction)
  }
  
  inputArguments <- as.list(environment())
  
  if(!is.data.frame(tuningParameters)) stop("tuningParameters must be a data frame (e.g., data.frame(lambda = seq(0,1,.1))).")
  
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
    stop("Invalid startingValues passed to regularizeSEMWithCustomPenalty. See ?controlQuasiNewtonBFGS for more information.")
  }
  
  # get parameters
  parameters <- aCV4SEM:::getParameters(SEM, raw = TRUE)
  
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
      diag(initialHessian) <- diag(initialHessian) + 1e-2
    }
    
  }else if(length(initialHessian) == 1 && is.numeric(initialHessian)){
    initialHessian <- diag(initialHessian,length(parameters))
    rownames(initialHessian) <- names(parameters)
    colnames(initialHessian) <- names(parameters)
  }else{
    stop("Invalid initialHessian passed to regularizeSEMWithCustomPenalty. See ?controlQuasiNewtonBFGS for more information.")
  }
  
  initialHessian4Optimizer <- initialHessian
  
  fits <- data.frame(
    "m2LL" = rep(NA,nrow(tuningParameters)),
    "regM2LL"= rep(NA,nrow(tuningParameters)),
    "convergence" = rep(NA,nrow(tuningParameters))
  )
  
  fits <- cbind(tuningParameters, fits)
  
  parameterEstimates <- as.data.frame(matrix(NA,nrow = nrow(tuningParameters), ncol = length(parameters)))
  colnames(parameterEstimates) <- names(parameters)
  parameterEstimates <- cbind(
    tuningParameters,
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
                                 max = nrow(tuningParameters), 
                                 initial = 0, 
                                 style = 3)
  }
  
  for(i in 1:nrow(tuningParameters)){
    
    if(control$verbose == 0){
      setTxtProgressBar(progressbar,i)
    }else{
      cat(paste0("\nIteration [", it, "/", nrow(tuningParameters),"]\n"))
    }
    
    currentTuningParameters <- tuningParameters[i,,drop = FALSE]
    
    result <- aCV4SEM:::quasiNewtonBFGS(SEM = SEM, 
                                        individualPenaltyFunction = individualPenaltyFunction, 
                                        individualPenaltyFunctionGradient = individualPenaltyFunctionGradient, 
                                        individualPenaltyFunctionHessian = individualPenaltyFunctionHessian,
                                        currentTuningParameters = currentTuningParameters,
                                        penaltyFunctionArguments = penaltyFunctionArguments,
                                        initialHessian = initialHessian4Optimizer,
                                        stepSize = control$stepSize,
                                        sig = control$sig,
                                        gam = control$gam,
                                        maxIterOut = control$maxIterOut,
                                        maxIterIn = control$maxIterIn,
                                        maxIterLine = control$maxIterLine,
                                        epsOut = control$epsOut,
                                        epsIn = control$epsIn,
                                        useMultipleConvergenceCriteria = control$useMultipleConvergenceCriteria,
                                        regM2LLChangeEps = control$regM2LLChangeEps,
                                        verbose = control$verbose)
    
    SEM <- result$SEM
    
    parameterEstimates[i, names(parameters)] <- aCV4SEM:::getParameters(result$SEM, raw = FALSE)[names(parameters)]
    fits$m2LL[i] <- result$m2LL
    fits$regM2LL[i] <- result$regM2LL
    fits$convergence[i] <- result$convergence
    
    initialHessian4Optimizer <- result$Hessian - individualPenaltyFunctionHessian(aCV4SEM:::getParameters(result$SEM, raw = FALSE)[names(parameters)], 
                                                                                  currentTuningParameters,
                                                                                  penaltyFunctionArguments)
    
    if(control$saveHessian) Hessians$Hessian[[it]] <- initialHessian4Optimizer # note: we are returning the Hessian only for the likelihood, not for the penalty
    
  }
  
  internalOptimization <- list(
    "HessiansOfDifferentiablePart" = Hessians
  )
  
  results <- new("regularizedSEMWithCustomPenalty",
                 parameters = parameterEstimates,
                 fits = fits,
                 parameterLabels = names(parameters),
                 internalOptimization = internalOptimization,
                 inputArguments = inputArguments)
  
  return(results)
  
}

#' regularizeSEMWithCustomPenaltyRsolnp
#' 
#' Optimize a SEM with custom penalty function using the Rsolnp optimizer (see ?Rsolnp::solnp). This optimizer is the default in regsem (see ?regsem::cv_regsem).
#' 
#' @param lavaanModel model of class lavaan 
#' @param individualPenaltyFunction penalty function which takes the current parameter values as first argument, the tuning parameters as second, and the penaltyFunctionArguments as third argument and 
#' returns a single value - the value of the penalty function for a single person. If the true penalty function is non-differentiable (e.g., lasso) a smooth
#' approximation of this function should be provided.
#' @param tuningParameters data.frame with tuning parameter values. Important: The function will iterate over the rows of these tuning parameters and pass them to your penalty function
#' @param penaltyFunctionArguments arguments passed to individualPenaltyFunction, individualPenaltyFunctionGradient, and individualPenaltyFunctionHessian
#' @param startingValues option to provide initial starting values. Only used for the first lambda. Three options are supported. Setting to "est" will use the estimates
#' from the lavaan model object. Setting to "start" will use the starting values of the lavaan model. Finally, a labeled vector with parameter
#' values can be passed to the function which will then be used as starting values.
#' @param saveHessian should the Hessian of the optimizer be saved? This will take a lot of space!
#' @param control option to set parameters of the optimizer; see ?Rsolnp::solnp
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
#' ## Defining a custom penalty function is a bit more complicated than
#' # using the default ones provided in regularizeSEM (see ?aCV4SEM::regularizeSEM).
#' # We start with the definition of the penalty function. Make sure that the derivatives
#' # of this function are well defined (i.e., the function is smooth)!
#' # We will use the lasso penalty as an example.
#' 
#' # The penalty function MUST accept three arguments; how you use them is up to you.
#' 
#' # The first argument are the parameters. aCV4SEM will pass a vector with the current
#' # parameter values and their names to your function. The parameter labels will be
#' # exactly the same as those used by lavaan! So you can check them before with:
#' 
#' aCV4SEM::getLavaanParameters(lavaanModel)
#' 
#' # The vector passed to your function will look like the vector above.
#' 
#' # The second argument is called tuningParameters MUST be a data.frame with the
#' # tuning-parameters. aCV4SEM will automatically iterate over the rows of this object.
#' # In case of LASSO regularization there is only one tuning parameter: lambda.
#' # Therefore, we specify the tuningParameters object as:
#' 
#' tuningParameters <- data.frame(lambda = seq(0,1,.1)) # we will test 11 lambdas here
#' print(tuningParameters)
#' 
#' # The third argument is called penaltyFunctionArguments and we can pass anything
#' # we want here. For the lasso penalty, we need two additional things:
#' # 1) We need the smoothing-parameter epsilon, which makes sure that our
#' # penalty is differentiable
#' # 2) We need the labels of the regularized parameters, so that we can only
#' # penalize those while all others remain unpenalized.
#' 
#' penaltyFunctionArguments <- list(
#'   eps = 1e-10,
#'   regularizedParameterLabels = paste0("l", 6:15)
#' )
#' 
#' # Now, it is time to specify our custom penalty function:
#' 
#' smoothLASSO <- function(
    #'   # here are our three arguments:
#'   parameters,
#'   tuningParameters,
#'   penaltyFunctionArguments
#' ){
#'   # to make it easier to see what is going on:
#'   lambda <- tuningParameters$lambda # tuningParameters will be ONE ROW OF the
#'   # tuningParameters object we created before -> it will contain just
#'   # one lambda in this case!
#'   eps <- penaltyFunctionArguments$eps
#'   regularizedParameterLabels <- penaltyFunctionArguments$regularizedParameterLabels
#' 
#'   regularizedParameters <- parameters[regularizedParameterLabels]
#' 
#'   # now, let's define our penalty function:
#'   penaltyLasso <- lambda*sum(sqrt(regularizedParameters^2 + eps))
#' 
#'   return(penaltyLasso)
#' }
#' # Important: This penalty function is assumed to be for a single individual only.
#' # aCV4SEM will multiply it with sample size N to get the penalty value of the
#' # full sample!
#' 
#' #### Now we are ready to optimize! ####
#' regsemApprox <- regularizeSEMWithCustomPenaltyRsolnp(lavaanModel = lavaanModel,
#'                                                individualPenaltyFunction = smoothLASSO,
#'                                                tuningParameters = tuningParameters,
#'                                                penaltyFunctionArguments = penaltyFunctionArguments)
#' 
#' # let's compare the results to an exact optimization:
#' 
#' regsemExact <- regularizeSEM(
#'   lavaanModel = lavaanModel,
#'   regularizedParameterLabels = paste0("l", 6:15),
#'   penalty = "lasso",
#'   lambdas = tuningParameters$lambda)
#' 
#' head(regsemExact@parameters[,regsemExact@parameterLabels] -
#'        regsemApprox@parameters[,regsemExact@parameterLabels])
#' # Note that the parameter estimates are basically identical.
#' 
#' ## To select a model, we used the approximate cross-validation function:
#' aCV <- aCV4regularizedSEMWithCustomPenalty(regularizedSEMWithCustomPenalty = regsemApprox,
#'                                            k = nrow(dataset))
#' plot(aCV)
#' # To extract the best parameter estimates:
#' coef(aCV)
#' @export
regularizeSEMWithCustomPenaltyRsolnp <- function(lavaanModel, 
                                                 individualPenaltyFunction,
                                                 tuningParameters,
                                                 penaltyFunctionArguments,
                                                 startingValues = "est",
                                                 saveHessian = FALSE,
                                                 control = list("trace" = 0)){
  
  inputArguments <- as.list(environment())
  
  if(!is.data.frame(tuningParameters)) stop("tuningParameters must be a data frame (e.g., data.frame(lambda = seq(0,1,.1))).")
  
  if(!is(lavaanModel, "lavaan")){
    stop("lavaanModel must be of class lavaan")
  }
  
  if(lavaanModel@Options$estimator != "ML") stop("lavaanModel must be fit with ml estimator.")
  
  rawData <- try(lavaan::lavInspect(lavaanModel, "data"))
  if(is(rawData, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  
  sampleSize <- nrow(rawData)
  
  ### initialize model ####
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
    stop("Invalid startingValues passed to regularizeSEMWithCustomPenaltyRsolnp.")
  }
  
  # get parameters
  parameters <- aCV4SEM:::getParameters(SEM, raw = TRUE)
  
  # define fitting function
  fitfun <- function(parameters, 
                     SEM, 
                     sampleSize,
                     individualPenaltyFunction, 
                     tuningParameters,
                     penaltyFunctionArguments){
    SEM <- try(aCV4SEM:::setParameters(SEM = SEM, labels = names(parameters), values = parameters, raw = TRUE), silent = TRUE)
    if(is(SEM, "try-error")){
      return(99999999999999999)
    }
    SEM <- try(aCV4SEM:::fit(SEM), silent = TRUE)
    if(is(SEM, "try-error") || !is.finite(SEM$m2LL)){
      return(99999999999999999)
    }
    
    penalty <- sampleSize*individualPenaltyFunction(parameters, tuningParameters, penaltyFunctionArguments)
    if(is(penalty, "try-error") || !is.finite(penalty)){
      return(99999999999999999)
    }
    
    return(SEM$m2LL + penalty)
  }
  
  fits <- data.frame(
    "m2LL" = rep(NA,nrow(tuningParameters)),
    "regM2LL"= rep(NA,nrow(tuningParameters)),
    "convergence" = rep(NA,nrow(tuningParameters))
  )
  
  parameterEstimates <- as.data.frame(matrix(NA,nrow = nrow(tuningParameters), ncol = length(parameters)))
  colnames(parameterEstimates) <- names(parameters)
  parameterEstimates <- cbind(
    tuningParameters,
    parameterEstimates
  )
  
  if(saveHessian){
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
  
  progressbar = txtProgressBar(min = 0, 
                               max = nrow(tuningParameters), 
                               initial = 0, 
                               style = 3)
  
  parametersInit <- parameters
  
  for(i in 1:nrow(tuningParameters)){
    
    setTxtProgressBar(progressbar,i)
    
    currentTuningParameters <- tuningParameters[i,,drop = FALSE]
    
    result <- Rsolnp::solnp(pars = parametersInit, fun = fitfun, SEM = SEM, 
                            sampleSize = sampleSize,
                            individualPenaltyFunction = individualPenaltyFunction, 
                            tuningParameters = currentTuningParameters,
                            penaltyFunctionArguments = penaltyFunctionArguments,
                            control = control)
    parametersInit <- result$pars
    
    SEM <- try(aCV4SEM:::setParameters(SEM = SEM, labels = names(parametersInit), values = parametersInit, raw = TRUE), silent = TRUE)
    SEM <- try(aCV4SEM:::fit(SEM), silent = TRUE)
    
    parameterEstimates[i, names(parameters)] <- aCV4SEM:::getParameters(SEM, raw = FALSE)[names(parameters)]
    fits$m2LL[i] <- SEM$m2LL
    fits$regM2LL[i] <- result$value[length(result$value)]
    fits$convergence[i] <- result$convergence == 0
    
    if(saveHessian) Hessians$Hessian[[it]] <- result$hessian # note: we are returning the Hessian for the likelihood AND the penalty
    
    if(result$convergence != 0) warning(paste0("Optimizer did not converge for ", paste0(names(currentTuningParameters), " = ",currentTuningParameters) , collapse = ", "))
  }
  
  internalOptimization <- list(
    "HessiansOfDifferentiablePart" = Hessians
  )
  
  results <- new("regularizedSEMWithCustomPenalty",
                 parameters = parameterEstimates,
                 fits = fits,
                 parameterLabels = names(parameters),
                 internalOptimization = internalOptimization,
                 inputArguments = inputArguments)
  
  return(results)
  
}