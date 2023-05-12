#' .regularizeSmoothSEMInternal
#' 
#' Internal function: This function computes the regularized models
#' for all smooth penalty functions which are implemented for bfgs.
#' Use the dedicated penalty functions (e.g., lessSEM::smoothLasso) to penalize
#' the model.
#' 
#' @param lavaanModel model of class lavaan 
#' @param penalty string: name of the penalty used in the model
#' @param weights labeled vector with weights for each of the parameters in the 
#' model.
#' @param tuningParameters data.frame with tuning parameter values
#' @param epsilon epsilon > 0; controls the smoothness of the approximation. Larger values = smoother 
#' @param tau parameters below threshold tau will be seen as zeroed
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param method optimizer used. Currently only "bfgs" is supported.
#' @param control used to control the optimizer. This element is generated with 
#' the controlBFGS function. See ?controlBFGS for more details.
#' @param notes option to pass a notes to function. All notes of the current
#' function will be added
#' @returns regularizedSEM
#' @keywords internal
.regularizeSmoothSEMInternal <- function(lavaanModel,
                                         penalty,
                                         weights,
                                         tuningParameters,
                                         epsilon,
                                         tau,
                                         method = "bfgs", 
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
  
  if(!penalty %in% c("ridge", "lasso", "adaptiveLasso", "elasticNet")) 
    stop(paste0(
      "bfgs only supports the following penalty functions: ",
      paste0(c("ridge", "lasso", "adaptiveLasso", "elasticNet"), collapse = ", ")
    )
    )
  
  if(!is(control, "controlBFGS")) 
    stop("control must be of class controlBfgs See ?controlBfgs.")
  
  control$breakOuter <- .adaptBreakingForWls(lavaanModel = lavaanModel, 
                                             currentBreaking = control$breakOuter,
                                             selectedDefault = control$breakOuter == controlBFGS()$breakOuter
  )
  
  if(!is.null(modifyModel$transformations)){
    if(any(control$initialHessian == "lavaan")){
      notes <- c(notes,
                 "Your model has transformations. Switching initialHessian from 'lavaan' to 'compute'.")
      control$initialHessian <- "compute"
    }
  }
  if(is.vector(lavaanModel)){
    if(any(control$initialHessian == "lavaan")){
      notes <- c(notes,
                 "You specified a multi-group model. Switching initialHessian from 'lavaan' to 'compute'.")
      control$initialHessian <- "compute"
    }
  }
  
  .checkLavaanModel(lavaanModel = lavaanModel)
  
  ## Setup Multi-Core ##
  
  .setupMulticore(control)
  
  ## check starting values ##
  
  startingValues <- control$startingValues
  if(!any(startingValues == "est") & penalty == "adaptiveLasso" & !is.numeric(weights)){
    createAdaptiveLassoWeights <- TRUE
  }else{
    createAdaptiveLassoWeights <- FALSE
  }
  
  ### initialize model ####
  if(is(lavaanModel, "lavaan")){
    SEM <- .initializeSEMForRegularization(lavaanModel = lavaanModel,
                                           startingValues = startingValues,
                                           modifyModel = modifyModel)
  }else{
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
  initialHessian <- .computeInitialHessian(initialHessian = control$initialHessian, 
                                           rawParameters = rawParameters, 
                                           lavaanModel = lavaanModel, 
                                           SEM = SEM,
                                           addMeans = modifyModel$addMeans,
                                           stepSize = control$stepSize,
                                           notes = notes)
  notes <- initialHessian$notes
  control$initialHessian <- initialHessian$initialHessian
  
  #### prepare regularized model object ####
  
  controlIntern <- list(
    epsilon = epsilon,
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
  
  if(is(SEM, "Rcpp_SEMCpp")){
    regularizedModel <- new(bfgsEnetSEM, 
                            weights, 
                            controlIntern)
  }else if(is(SEM, "Rcpp_mgSEM")){
    regularizedModel <- new(bfgsEnetMgSEM, 
                            weights, 
                            controlIntern)
  }
  
  #### define tuning parameters and prepare fit results ####
  
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
  
  if(control$saveDetails){
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
    
    result <- try(regularizedModel$optimize(rawParameters,
                                            SEM,
                                            tuningParameters$lambda[it],
                                            tuningParameters$alpha[it])
    )
    
    if(is(result, "try-error")) next
    
    rawParameters <- result$rawParameters
    fits$nonZeroParameters[it] <- length(rawParameters) - 
      sum(abs(rawParameters[weights[names(rawParameters)] != 0]) <= tau)
    
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
    
    if(control$saveDetails) 
      Hessians$Hessian[[it]] <- result$Hessian
    
    # set initial values for next iteration
    if(is(SEM, "try-Error")){
      # reset
      stop(paste0("Fit for ",
                  paste0(names(tuningParameters),
                         tuningParameters[it,], 
                         sep = " = "),
                  " resulted in Error!"))
      
      SEM <- .SEMFromLavaan(lavaanModel = lavaanModel,
                            whichPars = startingValues,
                            addMeans = control$addMeans)
      
      
      regularizedModel$setHessian(controlIntern$initialHessian)
      
      
    }else{
      
      
      if(control$saveDetails) Hessians$Hessian[[it]] <- result$Hessian
      
      # set Hessian for next iteration
      regularizedModel$setHessian(result$Hessian)
      
    }
    
  }
  
  internalOptimization <- list(
    "isMultiGroup" = is(SEM, "Rcpp_mgSEM"),
    "HessiansOfDifferentiablePart" = Hessians,
    "N" = SEM$sampleSize,
    "estimator"= SEM$getEstimator()
  )
  
  notes <- unique(notes)
  
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
  
  if(printNotes & length(notes) > 1){
    cat("\n")
    rlang::inform(
      notes
    )
  }
  
  return(results)
  
}


#' newTau
#' 
#' assign new value to parameter tau used by approximate optimization. Any regularized
#' value below tau will be evaluated as zeroed which directly impacts the AIC, BIC, etc.
#' 
#' @param regularizedSEM object fitted with approximate optimization
#' @param tau new tau value
#' @examples 
#' library(lessSEM)
#' 
#' # Identical to regsem, lessSEM builds on the lavaan
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
#' # Regularization:
#' 
#' lsem <- smoothLasso(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   epsilon = 1e-10,
#'   tau = 1e-4,
#'   lambdas = seq(0,1,length.out = 50))
#' newTau(regularizedSEM = lsem, tau = .1)
#' @returns regularizedSEM, but with new regularizedSEM@fits$nonZeroParameters
#' @export
newTau <- function(regularizedSEM, tau){
  if(!is(regularizedSEM,"regularizedSEM")) stop("regularizedSEM must be of class regularizedSEM")
  if(! regularizedSEM@penalty %in% c("lasso", "adaptiveLasso", "elasticNet")) stop("penalty must be of type lasso, adaptiveLasso, or elasticNet")
  if(is.null(regularizedSEM@inputArguments$tau)) stop("Could not find tau in regularizedSEM. Did you use a smoothed penalty?")
  
  nParameters <- length(regularizedSEM@parameterLabels)
  
  regularized <- names(regularizedSEM@weights)[regularizedSEM@weights != 0]
  
  if(any(!regularized %in% regularizedSEM@regularized)) stop("Error while identifying regularized parameters.")
  
  regularizedSEM@fits$nonZeroParameters <- 
    nParameters - 
    apply(regularizedSEM@parameters[,regularized],
          1,
          function(x) sum(abs(x) <= tau)
    )
  regularizedSEM@fits$nonZeroParameters[regularizedSEM@fits$lambda == 0] <- nParameters
  return(regularizedSEM)
}