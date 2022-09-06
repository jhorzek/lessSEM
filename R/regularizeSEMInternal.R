#' .regularizeSEMInternal
#' 
#' Internal function: This function computes the regularized models
#' for all penaltiy functions which are implemented for glmnet and gist.
#' Use the dedicated penalty functions (e.g., lessSEM::lasso) to penalize
#' the model.
#' 
#' @param lavaanModel model of class lavaan 
#' @param penalty string: name of the penalty used in the model
#' @param weights labeled vector with weights for each of the parameters in the 
#' model.
#' @param tuningParameters data.frame with tuning parameter values
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related procedures
#' (currently gist).
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta() and controlGlmnet() functions.
#' @keywords internal
.regularizeSEMInternal <- function(lavaanModel,
                                   penalty,
                                   weights,
                                   tuningParameters,
                                   method, 
                                   modifyModel,
                                   control){
  
  inputArguments <- as.list(environment())
  
  if(! method %in% c("ista", "glmnet")) 
    stop("Currently ony methods = 'ista' and methods = 'glmnet' are supported")
  if(method == "glmnet" & !penalty %in% c("ridge", "lasso", "adaptiveLasso", "elasticNet")) 
    stop(paste0(
      "glmnet only supports the following penalty functions: ",
      paste0(c("ridge", "lasso", "adaptiveLasso", "elasticNet"), collapse = ", ")
    )
    )
  
  if(method == "ista" && !is(control, "controlIsta")) 
    stop("control must be of class controlIsta. See ?controlIsta.")
  if(method == "glmnet" && !is(control, "controlGlmnet")) 
    stop("control must be of class controlGlmnet See ?controlGlmnet")
  
  if(!is(lavaanModel, "lavaan"))
    stop("lavaanModel must be of class lavaan")
  
  if(lavaanModel@Options$estimator != "ML") 
    stop("lavaanModel must be fit with ml estimator.")
  
  rawData <- try(lavaan::lavInspect(lavaanModel, "data"))
  if(is(rawData, "try-error")) 
    stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  N <- nrow(rawData)
  
  ### initialize model ####
  startingValues <- control$startingValues
  
  if(!any(startingValues == "est") & penalty == "adaptiveLasso" & !is.numeric(weights)){
    createAdaptiveLassoWeights <- TRUE
  }else{
    createAdaptiveLassoWeights <- FALSE
  }
  
  if(any(startingValues == "est")){
    SEM <- .SEMFromLavaan(lavaanModel = lavaanModel,
                          whichPars = "est",
                          addMeans = modifyModel$addMeans, 
                          activeSet = modifyModel$activeSet,
                          dataSet = modifyModel$dataSet,
                          transformations = modifyModel$transformations)
  }else if(any(startingValues == "start")){
    SEM <- .SEMFromLavaan(lavaanModel = lavaanModel,
                          whichPars = "start",
                          addMeans = modifyModel$addMeans, 
                          activeSet = modifyModel$activeSet,
                          dataSet = modifyModel$dataSet,
                          transformations = modifyModel$transformations)
  }else if(is.numeric(startingValues)){
    
    if(!all(names(startingValues) %in% names(getLavaanParameters(lavaanModel))))
      stop("Parameter names of startingValues do not match those of the lavaan object. See lessSEM::getLavaanParameters(lavaanModel).")
    SEM <- .SEMFromLavaan(lavaanModel = lavaanModel,
                          whichPars = "start", 
                          fit = FALSE,
                          addMeans = modifyModel$addMeans, 
                          activeSet = modifyModel$activeSet,
                          dataSet = modifyModel$dataSet,
                          transformations = modifyModel$transformations)
    SEM <- .setParameters(SEM = SEM, labels = names(startingValues), 
                          values = startingValues, 
                          raw = FALSE)
    SEM <- try(.fit(SEM))
    if(is(SEM, "try-error") || !is.finite(SEM$m2LL)) 
      stop("Infeasible starting values.")
    
  }else{
    stop("Invalid startingValues passed to elasticNet. See e.g., ?controlIsta for more information.")
  }
  
  # get parameters in raw form
  startingValues <- .getParameters(SEM, raw = TRUE)
  rawParameters <- .getParameters(SEM, raw = TRUE)
  
  # set weights
  if(!is.numeric(weights)){
    regularized <- weights
    if(penalty == "adaptiveLasso"){
      
      if(createAdaptiveLassoWeights){
        control_s <- control
        control_s$verbose <- 0
        # optimize model: We set lambda = 0, so we get the MLE
        cat("Computing MLE for adaptive lasso...\n")
        MLE <- lasso(
          lavaanModel = lavaanModel,
          lambdas = 0, 
          regularized = regularized,
          method = method, 
          modifyModel = modifyModel,
          control = control_s
        )
        
        # MLE:
        MLEs <- unlist(MLE@parameters[,MLE@parameterLabels])
        weights <- 1/abs(MLEs)
        weights[!names(weights) %in% regularized] <- 0
      }else{
        
        weights <- 1/abs(startingValues)
        weights[!names(weights) %in% regularized] <- 0
        
      }
    }else{
      weights <- startingValues
      weights[] <- 0
      weights[regularized] <- 1
    }
  }
  
  # make sure that all regularized parameters are actually in the model
  if(any(! regularized %in% names(startingValues))){
    stop(
      paste0("The following parameter was specified as regularized, but could not be found in the model:\n",
                paste0(regularized[! regularized %in% names(startingValues)], collapse = ", ")),
         "\nThe model parameters are:\n",
         paste0(names(startingValues), collapse = ", ")
    )
  }
  
  if(any(grepl("~~", names(weights)) & weights != 0)) 
    warning("Be careful when regularizing variances. These are implemented with a log-transform in lessSEM and you may not get the results you expected.")
  
  # make sure that the weights are in the correct order
  if(is.null(names(weights))) stop("weights must have the same names as the parameters")
  if(length(weights) != length(rawParameters)) stop("weights must be of the same length as the parameter vector.")
  if(any(!is.numeric(weights))) stop("weights must be numeric")
  weights <- weights[names(rawParameters)]
  
  #### glmnet requires an initial Hessian ####
  if(method == "glmnet"){
    initialHessian <- control$initialHessian
    if(is.matrix(initialHessian) && nrow(initialHessian) == length(rawParameters) && ncol(initialHessian) == length(rawParameters)){
      
      if(!all(rownames(initialHessian) %in% names(getLavaanParameters(lavaanModel))) ||
         !all(colnames(initialHessian) %in% names(getLavaanParameters(lavaanModel)))
      ) stop("initialHessian must have the parameter names as rownames and colnames. See lessSEM::getLavaanParameters(lavaanModel).")
      
    }else if(any(initialHessian == "compute")){
      
      initialHessian <- .getHessian(SEM = SEM, raw = TRUE)
      
    }else if(any(initialHessian == "scoreBased")){
      
      scores <- .getScores(SEM = SEM, raw = TRUE)
      FisherInformation <- matrix(0, nrow = ncol(scores), ncol = ncol(scores))
      rownames(FisherInformation) <- colnames(FisherInformation) <- colnames(scores)
      for(score in 1:nrow(scores)){
        FisherInformation <- FisherInformation + t(-.5*scores[score,, drop = FALSE]) %*%(-.5*scores[score,, drop = FALSE]) # we are using the -2 log-Likelihood
      }
      
      initialHessian <- -2*(-FisherInformation) # negative log-likelihood
      # make symmetric; just in case...
      initialHessian <- .5*(initialHessian + t(initialHessian))
      
    }else if(length(initialHessian) == 1 && is.numeric(initialHessian)){
      initialHessian <- diag(initialHessian,length(rawParameters))
      rownames(initialHessian) <- names(rawParameters)
      colnames(initialHessian) <- names(rawParameters)
    }else{
      stop("Invalid initialHessian passed to glmnet See ?controlGlmnet for more information.")
    }
    
    if(any(eigen(initialHessian, only.values = )$values < 0)){
      # make positive definite
      # see https://nhigham.com/2021/02/16/diagonally-perturbing-a-symmetric-matrix-to-make-it-positive-definite/
      eigenValues = eigen(initialHessian, only.values = )$values
      diagMat = diag(-1.1*min(eigenValues), nrow(initialHessian), ncol(initialHessian))
      initialHessian = initialHessian +  diagMat
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
      breakOuter = N*control$breakOuter,
      breakInner = N*control$breakInner,
      convergenceCriterion = control$convergenceCriterion, 
      verbose = control$verbose
    )
    
    regularizedModel <- new(glmnetEnet, 
                            weights, 
                            controlIntern)
    
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
      
      regularizedModel <- new(istaEnet, 
                              weights, 
                              controlIntern)
      
    }else if(penalty == "cappedL1"){
      
      regularizedModel <- new(istaCappedL1, 
                              weights, 
                              controlIntern)
      
    }else if(penalty == "lsp"){
      
      regularizedModel <- new(istaLSP, 
                              weights, 
                              controlIntern)
      
    }else if(penalty  == "scad"){
      
      regularizedModel <- new(istaScad, 
                              weights, 
                              controlIntern)
      
    }else if(penalty == "mcp"){
      
      regularizedModel <- new(istaMcp, 
                              weights, 
                              controlIntern)
      
    }else{
      stop("Unknow penalty selected.")
    }
    
  }
  
  #### define tuning parameters and prepare fit results ####
  ## get max lambda ##
  if(!is.null(tuningParameters$nLambdas)){
    # for lasso type penalties, the maximal lambda value can be determined
    # automatically
    message(paste0(
      "Automatically selecting the maximal lambda value.\n",
      "Note: This may fail if a model with all regularized parameters set to zero is not identified.")
    )
    
    maxLambda <- .getMaxLambda_C(regularizedModel = regularizedModel,
                                 SEM = SEM,
                                 rawParameters = rawParameters,
                                 weights = weights,
                                 N = N)
    if(tuningParameters$reverse){
      tuningParameters <- data.frame(
        lambda = rev(seq(0,
                         maxLambda,
                         length.out = tuningParameters$nLambdas)),
        alpha = 1
      )
    }else{
      tuningParameters <- data.frame(
        lambda = seq(0,
                     maxLambda,
                     length.out = tuningParameters$nLambdas),
        alpha = 1
      )
    }
    
    inputArguments$tuningParameters = tuningParameters
    
  }
  
  fits <- data.frame(
    tuningParameters,
    "m2LL" = NA,
    "regM2LL"= NA,
    "nonZeroParameters" = NA,
    "convergence" = NA
  )
  
  parameterEstimates <- as.data.frame(matrix(NA,
                                             nrow = nrow(tuningParameters), 
                                             ncol = length(rawParameters)))
  colnames(parameterEstimates) <- names(rawParameters)
  parameterEstimates <- cbind(
    tuningParameters,
    parameterEstimates
  )
  
  if(!is.null(modifyModel$transformations)){
    transformationsAre <- .getParameters(SEM,
                                         raw = FALSE,
                                         transformations = TRUE)
    transformationsAre <- transformationsAre[!names(transformationsAre)%in%names(rawParameters)]
    
    transformations <- as.data.frame(matrix(NA,
                                            nrow = nrow(tuningParameters), 
                                            ncol = length(transformationsAre)))
    colnames(transformations) <- names(transformationsAre)
    transformations <- cbind(
      tuningParameters,
      transformations
    )
  }else{
    transformations <- data.frame()
  }
  
  if(method == "glmnet" && control$saveHessian){
    Hessians <- list(
      "lambda" = tuningParameters$lambda,
      "alpha" = tuningParameters$alpha,
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
      
      result <- try(regularizedModel$optimize(rawParameters,
                                              SEM,
                                              tuningParameters$lambda[it],
                                              tuningParameters$alpha[it])
      )
      
    }else if(penalty %in% c("lsp", "scad", "mcp")){
      
      result <- try(regularizedModel$optimize(rawParameters,
                                              SEM,
                                              tuningParameters$theta[it],
                                              tuningParameters$lambda[it])
      )
      
    }else if(penalty == "cappedL1"){
      
      result <- try(regularizedModel$optimize(rawParameters,
                                              SEM,
                                              tuningParameters$theta[it],
                                              tuningParameters$lambda[it],
                                              tuningParameters$alpha[it])
      )
      
    }
    
    if(is(result, "try-error")) next
    
    rawParameters <- result$rawParameters
    fits$nonZeroParameters[it] <- length(rawParameters) - 
      sum(rawParameters[weights[names(rawParameters)] != 0] == 0)
    fits$regM2LL[it] <- result$fit
    fits$convergence[it] <- result$convergence
    
    # get unregularized fit:
    SEM <- .setParameters(SEM, 
                          names(rawParameters), 
                          values = rawParameters, 
                          raw = TRUE)
    fits$m2LL[it] <- SEM$fit()
    # transform internal parameter representation to "natural" form
    transformedParameters <- .getParameters(SEM,
                                            raw = FALSE)
    parameterEstimates[it, 
                       names(rawParameters)] <- transformedParameters[names(rawParameters)]
    
    if(!is.null(modifyModel$transformations)){
      transformationsAre <- .getParameters(SEM,
                                           raw = FALSE,
                                           transformations = TRUE)
      transformationsAre <- transformationsAre[!names(transformationsAre)%in%names(rawParameters)]
      transformations[it,
                      names(transformationsAre)] <- transformationsAre
    }
    
    if(method == "glmnet" && control$saveHessian) 
      Hessians$Hessian[[it]] <- result$Hessian
    
    # set initial values for next iteration
    if(is(SEM, "try-Error")){
      # reset
      warning(paste0("Fit for ",
                     paste0(names(tuningParameters),
                            tuningParameters[it,], 
                            sep = " = "),
                     " resulted in Error!"))
      
      SEM <- .SEMFromLavaan(lavaanModel = lavaanModel,
                            whichPars = startingValues,
                            addMeans = control$addMeans)
      
      if(method == "glmnet"){
        regularizedModel$setHessian(controlIntern$initialHessian)
      }
      
    }else{
      
      if(method == "glmnet"){
        if(control$saveHessian) Hessians$Hessian[[it]] <- result$Hessian
        
        # set Hessian for next iteration
        regularizedModel$setHessian(result$Hessian)
        
      }
      
    }
    
  }
  
  internalOptimization <- list(
    "HessiansOfDifferentiablePart" = Hessians,
    "functionCalls" = SEM$functionCalls,
    "gradientCalls" = SEM$gradientCalls
  )
  
  results <- new("regularizedSEM",
                 penalty = penalty,
                 parameters = parameterEstimates,
                 fits = fits,
                 parameterLabels = names(rawParameters),
                 weights = weights,
                 regularized = names(weights)[weights!=0],
                 transformations = transformations,
                 internalOptimization = internalOptimization,
                 inputArguments = inputArguments)
  
  return(results)
  
}