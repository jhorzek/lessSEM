elasticNet <- function(lavaanModel,
                       weights,
                       lambdas,
                       alphas,
                       method = "ista", 
                       control = controlIsta()){
  
  inputArguments <- as.list(environment())
  
  if(method == "ista" && !is(control, "controlIsta")) 
    stop("control must be of class controlIsta. See ?controlIsta.")
  
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
  if(any(startingValues == "est")){
    SEM <- linr:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                transformVariances = TRUE,
                                whichPars = "est",
                                addMeans = control$addMeans, 
                                activeSet = control$activeSet)
  }else if(any(startingValues == "start")){
    SEM <- linr:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                transformVariances = TRUE,
                                whichPars = "start",
                                addMeans = control$addMeans, 
                                activeSet = control$activeSet)
  }else if(is.numeric(startingValues)){
    
    if(!all(names(startingValues) %in% names(linr::getLavaanParameters(lavaanModel)))) stop("Parameter names of startingValues do not match those of the lavaan object. See linr::getLavaanParameters(lavaanModel).")
    SEM <- linr:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                transformVariances = TRUE,
                                whichPars = "start", 
                                fit = FALSE,
                                addMeans = control$addMeans, 
                                activeSet = control$activeSet)
    SEM <- linr:::setParameters(SEM = SEM, labels = names(startingValues), value = startingValues, raw = FALSE)
    SEM <- try(linr:::fit(SEM))
    if(is(SEM, "try-error") || !is.finite(SEM$m2LL)) 
      stop("Infeasible starting values.")
    
  }else{
    stop("Invalid startingValues passed to elasticNet. See e.g., ?controlIsta for more information.")
  }
  
  # get parameters
  startingValues <- linr:::getParameters(SEM, raw = TRUE)
  rawParameters <- linr:::getParameters(SEM, raw = TRUE)
  
  # make sure that the weights are in the correct order
  if(is.null(names(weights))) stop("weights must have the same names as the parameters")
  if(length(weights) != length(rawParameters)) stop("weights must be of the same length as the parameter vector.")
  if(any(!is.numeric(weights))) stop("weights must be numeric")
  weights <- weights[names(rawParameters)]
  
  #### glmnet requires an initial Hessian ####
  if(method == "glmnet"){
    initialHessian <- control$initialHessian
    if(is.matrix(initialHessian) && nrow(initialHessian) == length(rawParameters) && ncol(initialHessian) == length(rawParameters)){
      
      if(!all(rownames(initialHessian) %in% names(linr::getLavaanParameters(lavaanModel))) ||
         !all(colnames(initialHessian) %in% names(linr::getLavaanParameters(lavaanModel)))
      ) stop("initialHessian must have the parameter names as rownames and colnames. See linr::getLavaanParameters(lavaanModel).")
      
    }else if(any(initialHessian == "compute")){
      
      initialHessian <- linr:::getHessian(SEM = SEM, raw = TRUE)
      
    }else if(any(initialHessian == "scoreBased")){
      
      scores <- linr:::getScores(SEM = SEM, raw = TRUE)
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
      initialHessian <- diag(initialHessian,length(rawParameters))
      rownames(initialHessian) <- names(rawParameters)
      colnames(initialHessian) <- names(rawParameters)
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
      maxIterOut = control$maxIterOut,
      maxIterIn = control$maxIterIn,
      breakOuter = control$breakOuter,
      convCritInner = control$convCritInner,
      sigma = control$sigma,
      stepSizeInheritance = control$stepSizeInheritance,
      verbose = control$verbose
    )
    
    regularizedModel <- new(istaEnet, 
                            weights, 
                            controlIntern)
  }
  
  #### define tuning parameters and prepare fit results ####
  tuningGrid <- expand.grid("lambda" = lambdas, 
                            "alpha" = alphas)
  
  fits <- data.frame(
    tuningGrid,
    "m2LL" = NA,
    "regM2LL"= NA,
    "nonZeroParameters" = NA,
    "convergence" = NA
  )
  
  parameterEstimates <- as.data.frame(matrix(NA,
                                             nrow = nrow(tuningGrid), 
                                             ncol = length(rawParameters)))
  colnames(parameterEstimates) <- names(rawParameters)
  parameterEstimates <- cbind(
    tuningGrid,
    parameterEstimates
  )
  
  if(method == "glmnet" && control$saveHessian){
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
  
  #### print progress ####
  if(control$verbose == 0){
    progressbar = txtProgressBar(min = 0, 
                                 max = nrow(tuningGrid), 
                                 initial = 0, 
                                 style = 3)
  }
  
  #### Iterate over all tuning parameter combinations and fit models ####
  
  for(it in 1:nrow(tuningGrid)){
    if(control$verbose == 0){
      setTxtProgressBar(progressbar,it)
    }else{
      cat(paste0("\nIteration [", it, "/", nrow(tuningGrid),"]\n"))
    }
    
    lambda <- tuningGrid$lambda[it]
    alpha <- tuningGrid$alpha[it]
    
    result <- try(regularizedModel$optimize(rawParameters,
                                            SEM,
                                            lambda,
                                            alpha)
    )
    if(is(result, "try-error")) next
    
    rawParameters <- result$rawParameters
    fits$regM2LL[it] <- result$fit
    fits$convergence[it] <- result$convergence
    
    SEM <- linr::setParameters(SEM, 
                               names(rawParameters), 
                               values = rawParameters, 
                               raw = TRUE)
    fits$m2LL[it] <- SEM$fit()
    transformedParameters <- linr:::getParameters(SEM,
                                                  raw = FALSE)
    
    parameterEstimates[it, names(rawParameters)] <- transformedParameters[names(rawParameters)]
    
    if(method == "glmnet" && control$saveHessian) 
      Hessians$Hessian[[it]] <- result$Hessian
    
    # set initial values for next iteration
    if(is(SEM, "try-Error")){
      # reset
      warning("Fit for lambda = ",
              lambda, "alpha = ", 
              alpha,
              " resulted in Error!")
      
      SEM <- linr:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                  transformVariances = TRUE,
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
    "HessiansOfDifferentiablePart" = Hessians
  )
  
  results <- new("regularizedSEM",
                 parameters = parameterEstimates,
                 fits = fits,
                 parameterLabels = names(rawParameters),
                 weights = weights,
                 regularized = names(weights)[weights!=0],
                 internalOptimization = internalOptimization,
                 inputArguments = inputArguments)
  
  return(results)
  
}