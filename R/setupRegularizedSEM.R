#' .initializeSEMForRegularization
#' 
#' initializes the internal C++ SEM for regularization functions
#' @param lavaanModel model of class lavaan
#' @param startingValues either set to est, start, or labeled vector with starting values
#' @param modifyModel user supplied model modifications
#' @return model to be used by the regularization procedure
#' @keywords internal
.initializeSEMForRegularization <- function(lavaanModel,
                                            startingValues,
                                            modifyModel){
  if(any(startingValues == "est")){
    SEM <- .SEMFromLavaan(lavaanModel = lavaanModel,
                          whichPars = "est",
                          addMeans = modifyModel$addMeans, 
                          activeSet = modifyModel$activeSet,
                          dataSet = modifyModel$dataSet,
                          transformations = modifyModel$transformations,
                          transformationList = modifyModel$transformationList,
                          transformationGradientStepSize = modifyModel$transformationGradientStepSize)
  }else if(any(startingValues == "start")){
    SEM <- .SEMFromLavaan(lavaanModel = lavaanModel,
                          whichPars = "start",
                          addMeans = modifyModel$addMeans, 
                          activeSet = modifyModel$activeSet,
                          dataSet = modifyModel$dataSet,
                          transformations = modifyModel$transformations,
                          transformationList = modifyModel$transformationList,
                          transformationGradientStepSize = modifyModel$transformationGradientStepSize)
  }else if(is.numeric(startingValues)){
    
    if(!all(names(startingValues) %in% names(getLavaanParameters(lavaanModel))))
      stop("Parameter names of startingValues do not match those of the lavaan object. See lessSEM::getLavaanParameters(lavaanModel).")
    SEM <- .SEMFromLavaan(lavaanModel = lavaanModel,
                          whichPars = "start", 
                          fit = FALSE,
                          addMeans = modifyModel$addMeans, 
                          activeSet = modifyModel$activeSet,
                          dataSet = modifyModel$dataSet,
                          transformations = modifyModel$transformations,
                          transformationList = modifyModel$transformationList,
                          transformationGradientStepSize = modifyModel$transformationGradientStepSize)
    SEM <- .setParameters(SEM = SEM, 
                          labels = names(startingValues), 
                          values = startingValues, 
                          raw = FALSE)
    SEM <- try(.fit(SEM))
    if(is(SEM, "try-error") || !is.finite(SEM$m2LL)) 
      stop("Infeasible starting values.")
    
  }else{
    stop("Invalid startingValues passed to elasticNet. See e.g., ?controlIsta for more information.")
  }
  
  return(SEM)
}

#' .initializeMultiGroupSEMForRegularization
#' 
#' initializes the internal C++ SEM for regularization functions
#' @param lavaanModels vector with models of class lavaan
#' @param startingValues either set to est, start, or labeled vector with starting values
#' @param modifyModel user supplied model modifications
#' @return model to be used by the regularization procedure
#' @keywords internal
.initializeMultiGroupSEMForRegularization <- function(lavaanModels,
                                                      startingValues,
                                                      modifyModel){
  
  if(any(startingValues == "est")){
    SEM <- .multiGroupSEMFromLavaan(lavaanModels = lavaanModels,
                                    whichPars = "est",
                                    addMeans = modifyModel$addMeans,
                                    transformations = modifyModel$transformations,
                                    transformationList = modifyModel$transformationList,
                                    transformationGradientStepSize = modifyModel$transformationGradientStepSize
    )
  }else if(any(startingValues == "start")){
    SEM <- .multiGroupSEMFromLavaan(lavaanModels = lavaanModels,
                                    whichPars = "start",
                                    addMeans = modifyModel$addMeans,
                                    transformations = modifyModel$transformations,
                                    transformationList = modifyModel$transformationList,
                                    transformationGradientStepSize = modifyModel$transformationGradientStepSize)
  }else if(is.numeric(startingValues)){
    
    SEM <- .multiGroupSEMFromLavaan(lavaanModels = lavaanModel,
                                    whichPars = "start", 
                                    fit = FALSE,
                                    addMeans = modifyModel$addMeans,
                                    transformations = modifyModel$transformations,
                                    transformationList = modifyModel$transformationList,
                                    transformationGradientStepSize = modifyModel$transformationGradientStepSize)
    SEM <- .setParameters(SEM = SEM, 
                          labels = names(startingValues), 
                          values = startingValues, 
                          raw = FALSE)
    SEM <- try(.fit(SEM))
    if(is(SEM, "try-error") || !is.finite(SEM$m2LL)) 
      stop("Infeasible starting values.")
    
  }else{
    stop("Invalid startingValues passed to elasticNet. See e.g., ?controlIsta for more information.")
  }
  
  return(SEM)
}

#' .initializeWeights
#' 
#' initialize the adaptive lasso weights
#' @param weights weight argument passed to function
#' @param penalty penalty used
#' @param createAdaptiveLassoWeights should adaptive lasso weights be created?
#' @param control list with control elements for optimizer
#' @param lavaanModel model of type lavaan
#' @param modifyModel list with model modifications
#' @param startingValues either set to est, start, or labeled vector with starting values
#' @param rawParameters raw parameters
#' @keywords internal
.initializeWeights <- function(weights, 
                               penalty, 
                               method, 
                               createAdaptiveLassoWeights, 
                               control, 
                               lavaanModel, 
                               modifyModel, 
                               startingValues, 
                               rawParameters){
  if(is.null(weights)){
    # using bfgs without regularization -> no weights necessary
    weights <- startingValues
    weights[] <- 0
    regularized <- ""
  }else if(!is.numeric(weights)){
    regularized <- weights
    
    if(any(!regularized %in% names(rawParameters))){
      stop("The parameter(s) ", paste0(regularized[!regularized %in% names(rawParameters)], sep = ", "),
           " were specified as regularized but could not be found in the model. The model parameters are ", 
           paste0(names(rawParameters), sep = ", "))
    }
    
    if(penalty == "adaptiveLasso"){
      
      if(createAdaptiveLassoWeights){
        
        if(method == "bfgs") {
          method <- "glmnet" # bfgs is used by smooth lasso
          control <- controlGlmnet()
        }
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
  
  if(any(grepl("~~", names(weights)) & weights != 0)) 
    warning("Be careful when regularizing variances. These are implemented with a log-transform in lessSEM and you may not get the results you expected.")
  
  # make sure that the weights are in the correct order
  if(is.null(names(weights))) stop("weights must have the same names as the parameters")
  if(length(weights) != length(rawParameters)) stop("weights must be of the same length as the parameter vector.")
  if(any(!is.numeric(weights))) stop("weights must be numeric")
  weights <- weights[names(rawParameters)]
  
  return(weights)
}


#' .computeInitialHessian
#' 
#' computes the initial Hessian used in the optimization. Because we use the parameter
#' estimates from lavaan as starting values, it typcially makes sense to just use the
#' Hessian of the lavaan model as initial Hessian
#' @param initialHessian option to provide an initial Hessian to the optimizer. 
#' Must have row and column names corresponding to the parameter labels. use 
#' getLavaanParameters(lavaanModel) to 
#' see those labels. If set to "scoreBased", the outer product of the scores 
#' will be used as an approximation 
#' (see https://en.wikipedia.org/wiki/Berndt%E2%80%93Hall%E2%80%93Hall%E2%80%93Hausman_algorithm).
#' If set to "compute", the initial hessian will be computed. If set to a single 
#' value, a diagonal matrix with the single value along the diagonal will be used.
#' The default is "lavaan" which extracts the Hessian from the lavaanModel. This Hessian
#' will typically deviate from that of the internal SEM represenation of lessSEM (due to
#' the transformation of the variances), but works quite well in practice.
#' @param rawParameters vector with raw parameters
#' @param lavaanModel lavaan model object
#' @param SEM internal SEM representation
#' @param addMeans should a mean structure be added to the model?
#' @return Hessian matrix
#' @keywords internal
.computeInitialHessian <- function(initialHessian, rawParameters, lavaanModel, SEM, addMeans){
  
  if(is.matrix(initialHessian) && nrow(initialHessian) == length(rawParameters) && ncol(initialHessian) == length(rawParameters)){
    
    if(!all(rownames(initialHessian) %in% names(getLavaanParameters(lavaanModel))) ||
       !all(colnames(initialHessian) %in% names(getLavaanParameters(lavaanModel)))
    ) stop("initialHessian must have the parameter names as rownames and colnames. See lessSEM::getLavaanParameters(lavaanModel).")
    
  }
  
  if(any(initialHessian == "lavaan")){
    
    lavaanParameters <- getLavaanParameters(lavaanModel) 
    if((!lavaanModel@Options$meanstructure) && addMeans){
      .printNote("Your lavaan model has no mean structure. Switching initialHessian from 'lavaan' to 'compute'.")
      initialHessian <- "compute"
    }else if(!lavaanModel@Options$do.fit){
      .printNote("Your lavaan model was not optimized. Switching initialHessian from 'lavaan' to 'compute'.")
      initialHessian <- "compute"
      
    }else if(any(!names(lavaanParameters) %in% names(rawParameters)) |
             any(!names(rawParameters) %in% names(lavaanParameters))
    ){
      .printNote("Your model seems to have transformations. Switching initialHessian from 'lavaan' to 'compute'.")
      initialHessian <- "compute"
    }else{
      lavaanVcov <- try(lavaan:::vcov(lavaanModel),
                        silent = TRUE)
      if(is(lavaanVcov, "try-error")){
        warning("Could not extract initial Hessian from lavaan. Switching to ",
                "initialHessian = 'compute'.")
        initialHessian <- "compute"
        
      }else{
        
        lavaanVcov <- lavaanVcov[!duplicated(rownames(lavaanVcov)), 
                                 !duplicated(colnames(lavaanVcov))][names(rawParameters),
                                                                    names(rawParameters)]
        initialHessian <- 2*solve(lavaanVcov)
        
        if(any(eigen(initialHessian, only.values = TRUE)$values <= 0)){
          # make positive definite
          # see https://nhigham.com/2021/02/16/diagonally-perturbing-a-symmetric-matrix-to-make-it-positive-definite/
          eigenValues = eigen(initialHessian, only.values = )$values
          diagMat = diag(-1.1*min(eigenValues), nrow(initialHessian), ncol(initialHessian))
          initialHessian = initialHessian +  diagMat
        }
        
        return(initialHessian)
      }
    }
    
  }
  
  
  if(any(initialHessian == "compute")){
    
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
  
  if(any(eigen(initialHessian, only.values = TRUE)$values <= 0)){
    # make positive definite
    # see https://nhigham.com/2021/02/16/diagonally-perturbing-a-symmetric-matrix-to-make-it-positive-definite/
    eigenValues = eigen(initialHessian, only.values = TRUE)$values
    diagMat = diag(-1.1*min(min(eigenValues),-2), nrow(initialHessian), ncol(initialHessian))
    initialHessian = initialHessian +  diagMat
  }
  
  return(initialHessian)
}

#' .createTransformations
#' 
#' compiles the transformation function and adapts the parameterTable
#' @param transformations string with transformations
#' @param parameterLabels labels of parameteres already in the model
#' @param parameterTable existing parameter table
#' @return list with parameterTable and transformation function pointer
.createTransformations <- function(transformations, parameterLabels, parameterTable){
  if(!is(transformations, "character")) stop("transformations must be a string.")
  
  transformationFunctions <- .compileTransformations(syntax = transformations,
                                                     parameterLabels = parameterLabels)
  
  transformationFunctionPointer <- transformationFunctions$getPtr()
  
  # let's check for new parameters which have to be added to the model:
  
  addParameters <- transformationFunctions$parameters[!transformationFunctions$parameters %in% parameterTable$label]
  
  if(length(addParameters) != 0)
    parameterTable <- rbind(
      parameterTable,
      data.frame(
        "label" = addParameters,
        "location" = "transformation", 
        "row" = 1, 
        "col" = 1,
        "value" = 0,
        "rawValue" = 0,
        "isTransformation" = FALSE
      )
    )
  parameterTable$isTransformation[parameterTable$label %in% transformationFunctions$isTransformation] <- TRUE
  
  if(all(!is.na(transformationFunctions$startingValues))){
    
    for(i in 1:length(transformationFunctions$startingValues)){
      
      whichPar <- parameterTable$label == names(transformationFunctions$startingValues)[i]
      
      parameterTable$rawValue[whichPar] <- transformationFunctions$startingValues[i]
      
      parameterTable$value[whichPar] <- transformationFunctions$startingValues[i]
      
    }
  }
  
  return(
    list(transformationFunctionPointer = transformationFunctionPointer,
         parameterTable = parameterTable)
  )
}

#' .createMultiGroupTransformations
#' 
#' compiles the transformation function and adapts the parameter vector
#' @param transformations string with transformations
#' @param parameterValues values of parameters already in the model
#' @return list with extended parameter vector and transformation function pointer
.createMultiGroupTransformations <- function(transformations, parameterValues){
  if(!is(transformations, "character")) stop("transformations must be a string.")
  
  isTransformation <- rep(FALSE, length(parameterValues))
  names(isTransformation) <- names(parameterValues)
  
  transformationFunctions <- .compileTransformations(syntax = transformations,
                                                     parameterLabels = names(parameterValues))
  
  transformationFunctionPointer <- transformationFunctions$getPtr()
  
  # let's check for new parameters which have to be added to the model:
  
  addParameters <- transformationFunctions$parameters[!transformationFunctions$parameters %in% names(parameterValues)]
  
  if(length(addParameters) != 0){
    addParametersVec <- rep(FALSE, length(addParameters))
    names(addParametersVec) <- addParameters
    isTransformation <- c(isTransformation, addParametersVec)
    
    addParametersVec <- rep(0, length(addParameters))
    names(addParametersVec) <- addParameters
    parameterValues <- c(parameterValues, addParametersVec)
  }
  
  isTransformation[names(isTransformation) %in% transformationFunctions$isTransformation] <- TRUE
  
  if(all(!is.na(transformationFunctions$startingValues))){
    
    for(i in 1:length(transformationFunctions$startingValues)){
      
      whichPar <- names(parameterValues) == names(transformationFunctions$startingValues)[i]
      
      parameterValues[whichPar] <- transformationFunctions$startingValues[i]
      
    }
  }
  
  return(
    list(transformationFunctionPointer = transformationFunctionPointer,
         isTransformation = isTransformation,
         parameterValues = parameterValues)
  )
}

#' .checkLavaanModel
#' 
#' checks model of type lavaan
#' 
#' @param lavaanModel m0del of type lavaan
#' @return nothing
.checkLavaanModel <- function(lavaanModel){
  
  if(is.vector(lavaanModel)){
    for(m in lavaanModel){
      .checkLavaanModel(m)
    }
  }else{
    
    if(is(lavaanModel, "lavaan")){
      
      if(!is(lavaanModel, "lavaan"))
        stop("lavaanModel must be of class lavaan")
      
      if(lavaanModel@Options$estimator != "ML") 
        stop("lavaanModel must be fit with ml estimator.")
      
    }
  }
}