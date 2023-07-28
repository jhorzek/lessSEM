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
    if(is(SEM, "try-error") || !is.finite(SEM$objectiveValue)) 
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
    
    SEM <- .multiGroupSEMFromLavaan(lavaanModels = lavaanModels,
                                    whichPars = "start", 
                                    fit = FALSE,
                                    addMeans = modifyModel$addMeans,
                                    transformations = modifyModel$transformations,
                                    transformationList = modifyModel$transformationList,
                                    transformationGradientStepSize = modifyModel$transformationGradientStepSize)
    SEM <- .setParameters(SEM = SEM, 
                          labels = names(startingValues), 
                          values = startingValues, 
                          raw = TRUE)
    SEM <- try(.fit(SEM))
    if(is(SEM, "try-error") || !is.finite(SEM$objectiveValue)) 
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
#' @return vector with weights
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
#' @param stepSize initial step size
#' @param notes option to pass a notes to function. All notes of the current
#' function will be added
#' @return Hessian matrix and notes
#' @keywords internal
.computeInitialHessian <- function(initialHessian, rawParameters, lavaanModel, SEM, addMeans, stepSize, notes = NULL){
  
  if(is.null(notes)){
    printNotes <- TRUE
    notes <- c("Notes:")
  }else{
    # notes already exists and we only append the new ones
    printNotes <- FALSE
  }
  
  if(is.matrix(initialHessian) && nrow(initialHessian) == length(rawParameters) && ncol(initialHessian) == length(rawParameters)){
    
    if(!all(rownames(initialHessian) %in% names(rawParameters)) ||
       !all(colnames(initialHessian) %in% names(rawParameters))
    ){
      if(!all(rownames(initialHessian) %in% names(getLavaanParameters(lavaanModel))) ||
         !all(colnames(initialHessian) %in% names(getLavaanParameters(lavaanModel)))
      ) stop("initialHessian must have the parameter names as rownames and colnames. See lessSEM::getLavaanParameters(lavaanModel).")
    }
    
    if(any(eigen(initialHessian)$values < 0))
      stop("Provided initial Hessian is not positive definite.")
    
    return(list(initialHessian = initialHessian, notes = notes))
  }
  
  if(any(initialHessian == "lavaan")){
    
    lavaanParameters <- getLavaanParameters(lavaanModel) 
    if((!lavaanModel@Options$meanstructure) && addMeans){
      notes <- c(notes,
                 "Your lavaan model has no mean structure. Switching initialHessian from 'lavaan' to 'compute'.")
      initialHessian <- "compute"
    }else if(!lavaanModel@Options$do.fit){
      notes <- c(notes,
                 "Your lavaan model was not optimized. Switching initialHessian from 'lavaan' to 'compute'.")
      initialHessian <- "compute"
      
    }else if(any(!names(lavaanParameters) %in% names(rawParameters)) |
             any(!names(rawParameters) %in% names(lavaanParameters))
    ){
      notes <- c(notes,
                 "Your model seems to have transformations. Switching initialHessian from 'lavaan' to 'compute'.")
      initialHessian <- "compute"
    }else{
      lavaanVcov <- suppressWarnings(try(lavaan::vcov(lavaanModel),
                                         silent = TRUE))
      if(is(lavaanVcov, "try-error") | is.null(lavaanVcov)){
        notes <- c(notes,
                   "Could not extract initial Hessian from lavaan. Switching to initialHessian = 'compute'.")
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
        
        return(list(initialHessian = initialHessian, notes = notes))
      }
    }
    
  }
  
  
  if(any(initialHessian == "compute")){
    
    initialHessian <- .getHessian(SEM = SEM, raw = TRUE)
    
  }else if(any(initialHessian == "gradNorm")){
    
    # Adapted from Optim.jl
    # see https://github.com/JuliaNLSolvers/Optim.jl/blob/f43e6084aacf2dabb2b142952acd3fbb0e268439/src/multivariate/solvers/first_order/bfgs.jl#L104
    SEM$fit()
    
    gr <- .getGradients(SEM = SEM, raw = TRUE)
    
    initialHessian <- max(abs(gr)) * diag(length(gr)) * stepSize
    
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
  notes <- unique(notes)
  if(printNotes & length(notes) > 1){
    cat("\n")
    rlang::inform(notes)
  }
  
  return(
    list(
      initialHessian = initialHessian,
      notes = notes)
  )
}

#' .createTransformations
#' 
#' compiles the transformation function and adapts the parameterTable
#' @param transformations string with transformations
#' @param parameterLabels labels of parameteres already in the model
#' @param parameterTable existing parameter table
#' @return list with parameterTable and transformation function pointer
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
.checkLavaanModel <- function(lavaanModel){
  
  if(is.vector(lavaanModel)){
    for(m in lavaanModel){
      .checkLavaanModel(m)
    }
  }else{
    
    if(!is(lavaanModel, "lavaan"))
      stop("lavaanModel must be of class lavaan")
    if(lavaanModel@Options$categorical)
      stop("Categorical data is currently not supported.")
    if(lavaanModel@Options$missing %in% c("two.stage", "robust.two.stage"))
      stop(paste0("Missing = ", lavaanModel@Options$missing, " is currently not supported."))
    if(lavaanModel@Options$missing %in% c("fiml", "direct"))
      stop(paste0("Missing = ", lavaanModel@Options$missing, " was selected. Please change to missing = 'ml'."))
    if(lavaanModel@Options$missing %in% c("ml.x"))
      warning(paste0("Missing = ", lavaanModel@Options$missing, " was selected. We did not test this option yet; please check if it worked as expected."))
    if(lavaanModel@Options$effect.coding != "")
      stop(paste0("effect.coding = ", lavaanModel@Options$effect.coding, " is currently not supported."))
    if(lavaanModel@Options$ceq.simple)
      stop(paste0("ceq.simple = ", lavaanModel@Options$ceq.simple, " is currently not supported."))
    if(!is.null(lavaanModel@Options$group.label))
      stop("Lavaans multi-group models are currently not supported. See https://jhorzek.github.io/lessSEM/articles/Definition-Variables-and-Multi-Group-SEM.html for specifying multi-group models in lessSEM")
    if(!lavaanModel@Options$estimator %in% c("ML", "GLS", "WLS", "ULS", "MLM", "MLV", "MLMVS", "MLF", "MLR", "WLSM", "DWLS", "ULSMV"))
      stop(paste0("estimator = ", lavaanModel@Options$estimator, " is currently not supported."))
    if(lavaanModel@Options$likelihood %in% c("wishart"))
      stop(paste0("likelihood = ", lavaanModel@Options$likelihood, " is currently not supported."))
    if(lavaanModel@Options$bounds != "none")
      warning(paste0("bounds = ", lavaanModel@Options$bounds, " is currently not supported."))
    if(length(unlist(lavaanModel@Options$optim.bounds)) != 0)
      warning("optim.bounds is currently not supported.")
  }
}