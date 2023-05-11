#' .gpOptimizationInternal
#' 
#' Internal function: This function computes the regularized models
#' for all penaltiy functions which are implemented for glmnet and gist.
#' Use the dedicated penalty functions (e.g., lessSEM::gpLasso) to penalize
#' the model.
#' 
#' @param par labeled vector with starting values
#' @param weights labeled vector with weights for each of the parameters in the 
#' model.
#' @param fn R function which takes the parameters AND their labels 
#' as input and returns the fit value (a single value)
#' @param gr R function which takes the parameters AND their labels
#' as input and returns the gradients of the objective function. 
#' If set to NULL, numDeriv will be used to approximate the gradients 
#' @param additionalArguments additional argument passed to fn and gr
#' @param isCpp boolean: are fn and gr C++ function pointers?
#' @param penalty string: name of the penalty used in the model
#' @param tuningParameters data.frame with tuning parameter values
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related procedures
#' (currently gist).
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta() and controlGlmnet() functions.
#' @returns Object of class gpRegularized
#' @keywords internal
.gpOptimizationInternal <- function(par,
                                    weights,
                                    fn,
                                    gr = NULL,
                                    additionalArguments,
                                    isCpp = FALSE,
                                    penalty,
                                    tuningParameters,
                                    method,
                                    control){
  
  inputArguments <- as.list(environment())
  
  notes <- c("Notes:")
  
  if(! method %in% c("ista", "glmnet")) 
    stop("Currently ony methods = 'ista' and methods = 'glmnet' are supported")
  
  if(method == "ista" && !is(control, "controlIsta")) 
    stop("control must be of class controlIsta. See ?controlIsta.")
  if(method == "glmnet" && !is(control, "controlGlmnet")) 
    stop("control must be of class controlGlmnet See ?controlGlmnet")
  
  if(method == "glmnet" && any(control$initialHessian == "lavaan")) {
    notes <- c(notes,
               "Changing initialHessian from 'lavaan' to 'compute'")
    control$initialHessian <- "compute"
  }
  
  # gradient function
  if(is.null(gr)){
    
    # The following is necessary to ensure that the arguments are called 
    # correctly even if the user does not use the naming of our examples:
    additionalArguments$fn <- fn
    
    # gradient function using numDeriv:
    gr <- function(par, parameterLabels, additionalArguments){
      grad <- numDeriv::grad(func = additionalArguments$fn, 
                             x = par, 
                             parameterLabels = parameterLabels,
                             additionalArguments = additionalArguments)
      names(grad) <- names(par)
      return(grad)
    }
    
  }
  
  # check functions
  parameterLabels <- names(par)
  
  if(!isCpp){
    
    check_fn <- fn(par, parameterLabels, additionalArguments)
    
    if(length(check_fn) != 1) stop("fn returns more than one element!")
    check_fnIsColVec <- try(
      fn(matrix(par, ncol = 1), parameterLabels, additionalArguments),
      silent = TRUE
    )
    # we can also check if the function may expect a column vector instead of 
    # a row vector:
    
    check_fnIsRowVec <- try(
      fn(matrix(par, nrow = 1), parameterLabels, additionalArguments),
      silent = TRUE
    )
    if(is(check_fnIsRowVec, "try-error") & !is(check_fnIsColVec, "try-error")){
      fn_usr <- fn
      fn <- function(par, parameterLabels, additionalArguments
      ){
        return(fn_usr(matrix(par, ncol = 1), parameterLabels, additionalArguments))
      }
    }
  }
  
  if(!isCpp){
    
    check_fn <- gr(par, parameterLabels, additionalArguments)
    
    if(length(check_fn) != length(par)) stop("gr has different length than par!")
    
    # we can also check if the function may expect a column vector instead of 
    # a row vector:
    
    check_fnIsColVec <- try(
      gr(matrix(par, ncol = 1), parameterLabels, additionalArguments),
      silent = TRUE
    )
    check_fnIsRowVec <- try(
      gr(matrix(par, nrow = 1), parameterLabels, additionalArguments),
      silent = TRUE
    )
    if(is(check_fnIsRowVec, "try-error") & !is(check_fnIsColVec, "try-error")){
      gr_usr <- gr
      gr <- function(par, parameterLabels, additionalArguments
      ){
        return(gr_usr(matrix(par, ncol = 1), parameterLabels, additionalArguments))
      }
    }
  }
  
  # make sure that the weights are in the correct order
  par <- par
  # set weights
  if(!is.numeric(weights)){
    regularized <- weights
    if(penalty == "adaptiveLasso"){
      warning("Using starting values to construct adaptive lasso weights. Make sure that this is correct!")
      weights <- 1/abs(par)
      weights[!names(weights) %in% regularized] <- 0
    }else{
      weights <- par
      weights[] <- 0
      weights[regularized] <- 1
    }
  }
  
  if(is.null(names(weights))) stop("weights must have the same names as the parameters")
  if(length(weights) != length(par)) stop("weights must be of the same length as the parameter vector.")
  if(any(!is.numeric(weights))) stop("weights must be numeric")
  weights <- weights[names(par)]
  
  #### glmnet requires an initial Hessian ####
  if(method == "glmnet"){
    initialHessian <- control$initialHessian
    if(is.matrix(initialHessian) && 
       nrow(initialHessian) == length(par) && 
       ncol(initialHessian) == length(par)){
      
    }else if(any(initialHessian == "compute")){
      
      if(isCpp){
        initialHessian <- numDeriv::hessian(func = 
                                              function(par, parameterLabels, additionalArguments
                                              ){
                                                names(par) <- parameterLabels
                                                return(callFitFunction(fitFunctionSEXP = fn, 
                                                                       parameters = par, 
                                                                       userSuppliedElements = additionalArguments))
                                              }, 
                                            x = par,
                                            parameterLabels = parameterLabels,
                                            additionalArguments = additionalArguments)
      }else{
        initialHessian <- numDeriv::hessian(func = 
                                              function(par, parameterLabels, additionalArguments
                                              ){
                                                return(fn(par, parameterLabels, additionalArguments))
                                              }, 
                                            x = par,
                                            parameterLabels = parameterLabels,
                                            additionalArguments = additionalArguments)
      }
      
    }else if(length(initialHessian) == 1 && is.numeric(initialHessian)){
      initialHessian <- diag(initialHessian,length(par))
      rownames(initialHessian) <- names(par)
      colnames(initialHessian) <- names(par)
    }else{
      stop("Invalid initialHessian passed to GLMNET. See ?controlGLMNET for more information.")
    }
    
    control$initialHessian <- initialHessian
  }
  
  
  #### prepare regularized model object ####
  if(method == "glmnet"){
    
    controlIntern <- list(
      initialHessian = control$initialHessian,
      stepSize = control$stepSize,
      sigma = control$sigma,
      gamma = control$gamma,
      maxIterOut = control$maxIterOut,
      maxIterIn = control$maxIterIn,
      maxIterLine = control$maxIterLine,
      breakOuter = control$breakOuter,
      breakInner = control$breakInner,
      convergenceCriterion = control$convergenceCriterion, 
      verbose = control$verbose
    )
    
    if(isCpp){
      if(penalty %in% c("lasso", "elasticNet", "adaptiveLasso", "ridge")){
        regularizedModel <- new(glmnetEnetGeneralPurposeCpp, 
                                weights, 
                                controlIntern)
      }else{
        penaltyType <- rep("none", length(weights))
        penaltyType[weights != 0] <- penalty
        regularizedModel <- new(glmnetMixedPenaltyGeneralPurposeCpp, 
                                weights, 
                                penaltyType,
                                controlIntern)
      }
      
    }else{
      
      if(penalty %in% c("lasso", "elasticNet", "adaptiveLasso", "ridge")){
        regularizedModel <- new(glmnetEnetGeneralPurpose, 
                                weights, 
                                controlIntern)
      }else{
        penaltyType <- rep("none", length(weights))
        penaltyType[weights != 0] <- penalty
        regularizedModel <- new(glmnetMixedPenaltyGeneralPurpose, 
                                weights, 
                                penaltyType,
                                controlIntern)
      }
    }
    
  }else if(method == "ista"){
    
    controlIntern <- list(
      L0 = control$L0,
      eta = control$eta,
      accelerate = control$accelerate,
      maxIterOut = control$maxIterOut,
      maxIterIn = control$maxIterIn,
      breakOuter = control$breakOuter,
      convCritInner = control$convCritInner,
      sigma = control$sigma,
      stepSizeInheritance = control$stepSizeInheritance,
      verbose = control$verbose
    )
    
    if(penalty %in% c("ridge", "lasso", "adaptiveLasso", "elasticNet")){
      
      if(isCpp){
        regularizedModel <- new(istaEnetGeneralPurposeCpp, 
                                weights, 
                                controlIntern)
      }else{
        regularizedModel <- new(istaEnetGeneralPurpose, 
                                weights, 
                                controlIntern)
      }
      
      
    }else{
      
      
      if(isCpp){
        penaltyType <- rep("none", length(weights))
        penaltyType[weights != 0] <- penalty
        regularizedModel <- new(istaMixedPenaltyGeneralPurposeCpp, 
                                weights, 
                                penaltyType,
                                controlIntern)
      }else{
        penaltyType <- rep("none", length(weights))
        penaltyType[weights != 0] <- penalty
        regularizedModel <- new(istaMixedPenaltyGeneralPurpose, 
                                weights, 
                                penaltyType,
                                controlIntern)
      }
      
    }
  }
  
  #### define tuning parameters and prepare fit results ####
  ## get max lambda ##
  if(!is.null(tuningParameters$nLambdas)){
    # for lasso type penalties, the maximal lambda value can be determined
    # automatically
    notes <- c(notes,
               paste0(
                 "Automatically selecting the maximal lambda value. ",
                 "This may fail if a model with all regularized parameters set to zero is not identified.")
    )
    
    maxLambda <- .gpGetMaxLambda(regularizedModel,
                                 par,
                                 fn,
                                 gr,
                                 additionalArguments,
                                 weights)
    
    if(tuningParameters$reverse){
      tuningParameters <- data.frame(
        lambda = rev(curveLambda(maxLambda = maxLambda, 
                                 lambdasAutoCurve = tuningParameters$curve, 
                                 tuningParameters$nLambdas)),
        alpha = 1,
        theta = 0
      )
    }else{
      tuningParameters <- data.frame(
        lambda = curveLambda(maxLambda = maxLambda, 
                             lambdasAutoCurve = tuningParameters$curve, 
                             tuningParameters$nLambdas),
        alpha = 1,
        theta = 0
      )
    }
    
    inputArguments$tuningParameters = tuningParameters
    
  }
  
  fits <- data.frame(
    tuningParameters,
    "objectiveValue" = NA,
    "regObjectiveValue"= NA,
    "nonZeroParameters" = NA,
    "convergence" = NA
  )
  
  parameterEstimates <- as.data.frame(matrix(NA,
                                             nrow = nrow(tuningParameters), 
                                             ncol = length(par)))
  colnames(parameterEstimates) <- names(par)
  parameterEstimates <- cbind(
    tuningParameters,
    parameterEstimates
  )
  
  if(method == "glmnet" && control$saveDetails){
    Hessians <- list(
      "lambda" = tuningParameters$lambda,
      "alpha" = tuningParameters$alpha,
      "theta" = tuningParameters$theta,
      "Hessian" = lapply(1:nrow(tuningParameters), 
                         matrix, 
                         data= NA, 
                         nrow=nrow(initialHessian), 
                         ncol=ncol(initialHessian))
    )
  }else{
    Hessians <- list(NULL)
  }
  
  #### print progress ####
  if(control$verbose == 0){
    progressbar = utils::txtProgressBar(min = 0, 
                                        max = nrow(tuningParameters), 
                                        initial = 0, 
                                        style = 3)
  }
  
  #### Iterate over all tuning parameter combinations and fit models ####
  
  for(it in 1:nrow(tuningParameters)){
    if(control$verbose == 0){
      utils::setTxtProgressBar(progressbar,it)
    }else{
      cat(paste0("\nIteration [", it, "/", nrow(tuningParameters),"]\n"))
    }
    
    
    if(penalty %in% c("ridge", "lasso", "adaptiveLasso", "elasticNet")){
      
      result <- try(regularizedModel$optimize(par,
                                              fn,
                                              gr,
                                              additionalArguments,
                                              tuningParameters$lambda[it],
                                              tuningParameters$alpha[it])
      )
      
    }else{
      
      result <- try(regularizedModel$optimize( par,
                                               fn,
                                               gr,
                                               additionalArguments,
                                               rep(tuningParameters$lambda[it], length(par)),
                                               rep(tuningParameters$theta[it], length(par)),
                                               rep(tuningParameters$alpha[it], length(par)))
      )
      
    }
    
    if(is(result, "try-error")) next
    
    rawParameters <- result$rawParameters
    parameterEstimates[it, names(rawParameters)] <- rawParameters[names(rawParameters)]
    
    fits$nonZeroParameters[it] <- length(rawParameters) - 
      sum(rawParameters[weights[names(rawParameters)] != 0] == 0)
    fits$regObjectiveValue[it] <- result$fit
    fits$convergence[it] <- result$convergence
    if(!isCpp) fits$objectiveValue[it] <- fn(rawParameters, parameterLabels, additionalArguments)
    
    if(method == "glmnet" && control$saveDetails) 
      Hessians$Hessian[[it]] <- result$Hessian
    
    # set initial values for next iteration
    
    
    if(method == "glmnet"){
      if(control$saveDetails) Hessians$Hessian[[it]] <- result$Hessian
      
      # set Hessian for next iteration
      regularizedModel$setHessian(result$Hessian)
      
    } 
    
    
  }
  
  internalOptimization <- list(
    "HessiansOfDifferentiablePart" = Hessians
  )
  
  notes <- unique(notes)
  if(length(notes) > 1){
    cat("\n")
    rlang::inform(notes)
  }
  
  results <- new("gpRegularized",
                 penalty = penalty,
                 parameters = parameterEstimates,
                 fits = fits,
                 parameterLabels = names(rawParameters),
                 weights = weights,
                 regularized = names(weights)[weights!=0],
                 internalOptimization = internalOptimization,
                 inputArguments = inputArguments)
  
  return(results)
  
}