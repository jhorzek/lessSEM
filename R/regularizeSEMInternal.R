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
#' @param notes option to pass a notes to function. All notes of the current
#' function will be added
#' @return regularized SEM 
#' @keywords internal
.regularizeSEMInternal <- function(lavaanModel,
                                   penalty,
                                   weights,
                                   tuningParameters,
                                   method, 
                                   modifyModel,
                                   control,
                                   notes = NULL){
  
  if(is.null(notes)){
    printNotes <- TRUE
    notes <- c("Notes:")
  }else{
    # notes already exists and we only append the new ones
    printNotes <- FALSE
  }
  
  inputArguments <- as.list(environment())
  
  if(! method %in% c("ista", "glmnet")) 
    stop("Currently ony methods = 'ista' and methods = 'glmnet' are supported")
  if(method == "glmnet" & !penalty %in% c("ridge", "lasso", "adaptiveLasso", "elasticNet", "scad", "cappedL1", "mcp", "lsp")) 
    stop(
      paste0(
        "glmnet only supports the following penalty functions: ",
        paste0(c("ridge", "lasso", "adaptiveLasso", "elasticNet", "scad", "cappedL1", "mcp", "lsp"), 
               collapse = ", ")
      )
    )
  
  if(method == "ista" && !is(control, "controlIsta")) 
    stop("control must be of class controlIsta. See ?controlIsta.")
  if(method == "glmnet" && !is(control, "controlGlmnet")) 
    stop("control must be of class controlGlmnet See ?controlGlmnet")
  
  control$breakOuter <- .adaptBreakingForWls(lavaanModel = lavaanModel, 
                                             currentBreaking = control$breakOuter,
                                             selectedDefault = ifelse(method == "ista",
                                                                      control$breakOuter == controlIsta()$breakOuter,
                                                                      control$breakOuter == controlGlmnet()$breakOuter
                                             ))
  
  
  if(method == "glmnet" && is.vector(lavaanModel)){
    if(any(control$initialHessian == "lavaan")){
      notes <- c(notes,
                 "You specified a multi-group model. Switching initialHessian from 'lavaan' to 'compute'.")
      control$initialHessian <- "compute"
    }
  }
  if(!is.null(modifyModel$transformations)){
    
    notes <- c(notes,
               paste("Your model has transformations. If you transform variances",
                     "the variance estimates returned by lessSEM may not be the true variances",
                     "but transformations thereof. Check model@transformations to find the actual variance",
                     "estimates for your regularized variances")
    )
    if(method == "glmnet"){
      if(any(control$initialHessian == "lavaan")){
        cat(" - Switching initialHessian from 'lavaan' to 'compute'.\n")
        control$initialHessian <- "compute"
      }
    }
  }
  
  .checkLavaanModel(lavaanModel = lavaanModel)
  
  ## Setup Multi-Core ##
  
  .setupMulticore(control)
  
  ### initialize model ####
  startingValues <- control$startingValues
  if(!any(startingValues == "est") & penalty == "adaptiveLasso" & !is.numeric(weights)){
    createAdaptiveLassoWeights <- TRUE
  }else{
    createAdaptiveLassoWeights <- FALSE
  }
  if(penalty == "adaptiveLasso" & !is.numeric(weights) & !is.null(modifyModel$transformations)){
    createAdaptiveLassoWeights <- TRUE
  }
  
  ### initialize model ####
  if(is(lavaanModel, "lavaan")){
    SEM <- .initializeSEMForRegularization(lavaanModel = lavaanModel,
                                           startingValues = startingValues,
                                           modifyModel = modifyModel)
  }else if(is.vector(lavaanModel)){
    SEM <- .initializeMultiGroupSEMForRegularization(lavaanModels = lavaanModel,
                                                     startingValues = startingValues,
                                                     modifyModel = modifyModel)
  }
  
  # check if we have a likelihood objective function:
  usesLikelihood <- all(SEM$getEstimator() == "fiml")
  
  N <- SEM$sampleSize
  
  # get parameters in raw form
  startingValues <- .getParameters(SEM, raw = TRUE)
  rawParameters <- .getParameters(SEM, raw = TRUE)
  
  # set weights
  weights <- .initializeWeights(weights = weights, 
                                penalty = penalty, 
                                method = method,
                                createAdaptiveLassoWeights = createAdaptiveLassoWeights, 
                                control = control, 
                                lavaanModel = lavaanModel, 
                                modifyModel = modifyModel, 
                                startingValues = startingValues, 
                                rawParameters = rawParameters)
  
  #### glmnet requires an initial Hessian ####
  if(method == "glmnet"){
    initialHessian <- .computeInitialHessian(initialHessian = control$initialHessian, 
                                             rawParameters = rawParameters, 
                                             lavaanModel = lavaanModel, 
                                             SEM = SEM,
                                             addMeans = modifyModel$addMeans,
                                             stepSize = control$stepSize,
                                             notes = notes)
    notes <- initialHessian$notes
    control$initialHessian <- initialHessian$initialHessian
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
    
    if(penalty %in% c("ridge", "lasso", "adaptiveLasso", "elasticNet")){
      if(is(SEM, "Rcpp_SEMCpp")){
        regularizedModel <- new(glmnetEnetSEM, 
                                weights, 
                                controlIntern)
      }else if(is(SEM, "Rcpp_mgSEM")){
        regularizedModel <- new(glmnetEnetMgSEM, 
                                weights, 
                                controlIntern)
      }
    }else if(penalty == "scad"){
      if(is(SEM, "Rcpp_SEMCpp")){
        regularizedModel <- new(glmnetScadSEM,
                                weights,
                                controlIntern)
      }else if(is(SEM, "Rcpp_mgSEM")){
        regularizedModel <- new(glmnetScadMgSEM,
                                weights,
                                controlIntern)
      }
    }else if(penalty == "cappedL1"){
      if(is(SEM, "Rcpp_SEMCpp")){
        regularizedModel <- new(glmnetCappedL1SEM,
                                weights,
                                controlIntern)
      }else if(is(SEM, "Rcpp_mgSEM")){
        regularizedModel <- new(glmnetCappedL1MgSEM,
                                weights,
                                controlIntern)
      }
    }else if(penalty == "mcp"){
      if(is(SEM, "Rcpp_SEMCpp")){
        regularizedModel <- new(glmnetMcpSEM,
                                weights,
                                controlIntern)
      }else if(is(SEM, "Rcpp_mgSEM")){
        regularizedModel <- new(glmnetMcpMgSEM,
                                weights,
                                controlIntern)
      }
    }else if(penalty == "lsp"){
      if(is(SEM, "Rcpp_SEMCpp")){
        regularizedModel <- new(glmnetLspSEM,
                                weights,
                                controlIntern)
      }else if(is(SEM, "Rcpp_mgSEM")){
        regularizedModel <- new(glmnetLspMgSEM,
                                weights,
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
      
      if(is(SEM, "Rcpp_SEMCpp")){
        regularizedModel <- new(istaEnetSEM, 
                                weights, 
                                controlIntern)
      }else if(is(SEM, "Rcpp_mgSEM")){
        regularizedModel <- new(istaEnetMgSEM, 
                                weights, 
                                controlIntern)
      }
      
    }else if(penalty == "cappedL1"){
      
      if(is(SEM, "Rcpp_SEMCpp")){
        regularizedModel <- new(istaCappedL1SEM, 
                                weights, 
                                controlIntern)
      }else if(is(SEM, "Rcpp_mgSEM")){
        regularizedModel <- new(istaCappedL1mgSEM, 
                                weights, 
                                controlIntern)
      }
      
    }else if(penalty == "lsp"){
      
      if(is(SEM, "Rcpp_SEMCpp")){
        regularizedModel <- new(istaLSPSEM, 
                                weights, 
                                controlIntern)
      }else if(is(SEM, "Rcpp_mgSEM")){
        regularizedModel <- new(istaLSPMgSEM, 
                                weights, 
                                controlIntern)
      }
      
    }else if(penalty  == "scad"){
      
      if(is(SEM, "Rcpp_SEMCpp")){
        regularizedModel <- new(istaScadSEM, 
                                weights, 
                                controlIntern)
      }else if(is(SEM, "Rcpp_mgSEM")){
        regularizedModel <- new(istaScadMgSEM, 
                                weights, 
                                controlIntern)
      }
      
    }else if(penalty == "mcp"){
      
      if(is(SEM, "Rcpp_SEMCpp")){
        regularizedModel <- new(istaMcpSEM, 
                                weights, 
                                controlIntern)
      }else if(is(SEM, "Rcpp_mgSEM")){
        regularizedModel <- new(istaMcpMgSEM, 
                                weights, 
                                controlIntern)
      }
      
    }else{
      stop("Unknow penalty selected.")
    }
    
  }
  
  #### define tuning parameters and prepare fit results ####
  ## get max lambda ##
  if(!is.null(tuningParameters$nLambdas)){
    # for lasso type penalties, the maximal lambda value can be determined
    # automatically
    notes <- c(
      notes,
      paste0(
        "Automatically selecting the maximal lambda value.",
        " This may fail if a model with all regularized parameters set to zero is not identified."
      )
    )
    
    maxLambda <- .getMaxLambda_C(regularizedModel = regularizedModel,
                                 SEM = SEM,
                                 rawParameters = rawParameters,
                                 weights = weights,
                                 N = N)
    if(tuningParameters$reverse){
      tuningParameters <- data.frame(
        lambda = rev(curveLambda(maxLambda = maxLambda, 
                                 lambdasAutoCurve = tuningParameters$curve, 
                                 tuningParameters$nLambdas)),
        alpha = 1
      )
    }else{
      tuningParameters <- data.frame(
        lambda = curveLambda(maxLambda = maxLambda, 
                             lambdasAutoCurve = tuningParameters$curve, 
                             tuningParameters$nLambdas),
        alpha = 1
      )
    }
    
    inputArguments$tuningParameters = tuningParameters
    
  }
  
  fits <- data.frame(
    tuningParameters,
    "objectiveValue" = NA,
    "regObjectiveValue" = NA,
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
  
  if(method == "glmnet" && control$saveDetails){
    Hessians <- list(
      "lambda" = tuningParameters$lambda,
      "alpha" = tuningParameters$alpha,
      "Hessian" = lapply(1:nrow(tuningParameters), 
                         matrix, 
                         data= NA, 
                         nrow=nrow(control$initialHessian), 
                         ncol=ncol(control$initialHessian))
    )
  }else{
    Hessians <- list(NULL)
  }
  
  # save model implied matrices
  implied <- list(means = vector("list", nrow(tuningParameters)),
                  covariances = vector("list", nrow(tuningParameters)))
  
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
      if(method == "ista")
        result <- try(regularizedModel$optimize(rawParameters,
                                                SEM,
                                                tuningParameters$theta[it],
                                                tuningParameters$lambda[it],
                                                tuningParameters$alpha[it])
        )
      if(method == "glmnet")
        result <- try(regularizedModel$optimize(rawParameters,
                                                SEM,
                                                tuningParameters$theta[it],
                                                tuningParameters$lambda[it])
        )
    }
    
    if(is(result, "try-error")) next
    
    rawParameters <- result$rawParameters
    fits$nonZeroParameters[it] <- length(rawParameters) - 
      sum(rawParameters[weights[names(rawParameters)] != 0] == 0)
    
    fits$regObjectiveValue[it] <- result$fit
    if(usesLikelihood)
      fits$regM2LL[it] <- fits$regObjectiveValue[it]
    
    fits$convergence[it] <- result$convergence
    
    # get unregularized fit:
    SEM <- .setParameters(SEM, 
                          names(rawParameters), 
                          values = rawParameters, 
                          raw = TRUE)
    
    fits$objectiveValue[it] <- SEM$fit()
    if(usesLikelihood)
      fits$m2LL[it] <- fits$objectiveValue[it]
    
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
    
    if(method == "glmnet" && control$saveDetails) 
      Hessians$Hessian[[it]] <- result$Hessian
    
    # save implied
    if(is(SEM, "Rcpp_SEMCpp") & control$saveDetails){
      implied$means[[it]] <- SEM$impliedMeans
      implied$covariances[[it]] <- SEM$impliedCovariance
      
      rownames(implied$means[[it]]) <- SEM$manifestNames
      dimnames(implied$covariances[[it]]) <- list(SEM$manifestNames,
                                                  SEM$manifestNames)
    }else if(is(SEM, "Rcpp_mgSEM")){
      implied$means[[it]] <- NULL
      implied$covariances[[it]] <- NULL
    }
    
    # set initial values for next iteration
    if(is(SEM, "try-Error")){
      # reset
      warning(paste0("Fit for ",
                     paste0(names(tuningParameters),
                            tuningParameters[it,], 
                            sep = " = "),
                     " resulted in Error!"))
      
      if(is(lavaanModel, "lavaan")){
        SEM <- .initializeSEMForRegularization(lavaanModel = lavaanModel,
                                               startingValues = startingValues,
                                               modifyModel = modifyModel)
      }else if(is.vector(lavaanModel)){
        SEM <- .initializeMultiGroupSEMForRegularization(lavaanModels = lavaanModel,
                                                         startingValues = startingValues,
                                                         modifyModel = modifyModel)
      }
      
      if(method == "glmnet"){
        regularizedModel$setHessian(controlIntern$initialHessian)
      }
      
    }else{
      
      if(method == "glmnet"){
        if(control$saveDetails) Hessians$Hessian[[it]] <- result$Hessian
        
        # set Hessian for next iteration
        regularizedModel$setHessian(result$Hessian)
        
      }
      
    }
    
  }
  
  if(is(SEM, "Rcpp_SEMCpp")) internalOptimization <- list(
    "isMultiGroup" = FALSE,
    "implied" = implied,
    "HessiansOfDifferentiablePart" = Hessians,
    "functionCalls" = SEM$functionCalls,
    "gradientCalls" = SEM$gradientCalls,
    "N" = SEM$sampleSize,
    "estimator"= SEM$getEstimator()
  )
  if(is(SEM, "Rcpp_mgSEM")) internalOptimization <- list(
    "isMultiGroup" = TRUE,
    "implied" = implied,
    "HessiansOfDifferentiablePart" = Hessians,
    "functionCalls" = NA,
    "gradientCalls" = NA,
    "N" = SEM$sampleSize,
    "estimator"= SEM$getEstimator()
  )
  
  if(any(internalOptimization$estimator == "wls")){
    notes <- c(notes, "WLS (and variants) is a very new feature and not yet thoroughly tested. Please be wary of bugs!")
  }
  
  results <- new("regularizedSEM",
                 penalty = penalty,
                 parameters = parameterEstimates,
                 fits = fits,
                 parameterLabels = names(rawParameters),
                 weights = weights,
                 regularized = names(weights)[weights!=0],
                 transformations = transformations,
                 internalOptimization = internalOptimization,
                 inputArguments = inputArguments,
                 notes = notes)
  
  notes <- unique(notes)
  
  if(printNotes & (length(notes) > 1)){
    cat("\n")
    rlang::inform(
      notes
    )
  }
  
  return(results)
  
}