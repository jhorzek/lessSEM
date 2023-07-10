#' Class for regularized SEM
#' @slot penalty penalty used (e.g., "lasso")
#' @slot parameters data.frame with parameter estimates
#' @slot fits data.frame with all fit results
#' @slot parameterLabels character vector with names of all parameters
#' @slot weights vector with weights given to each of the parameters in the penalty
#' @slot regularized character vector with names of regularized parameters
#' @slot transformations if the model has transformations, the transformed parameters 
#' are returned
#' @slot internalOptimization list of elements used internally
#' @slot inputArguments list with elements passed by the user to the general
#' @slot notes internal notes that have come up when fitting the model
#' @export
setClass(Class = "regularizedSEM",
         representation = representation(
           penalty = "character",
           parameters="data.frame",
           fits="data.frame", 
           parameterLabels = "character",
           weights = "numeric",
           regularized = "character",
           transformations = "data.frame",
           internalOptimization="list", 
           inputArguments="list",
           notes = "character"
         )
)

#' show
#' @param object object of class regularizedSEM
#' @return No return value, just prints estimates
#' @export
setMethod("show", "regularizedSEM", function (object) {
  #modelName <-deparse(substitute(object)) # get the name of the object
  isMultiGroup <- object@internalOptimization$isMultiGroup
  method <- object@inputArguments$method
  penalty <- object@penalty
  estimator <- object@internalOptimization$estimator
  
  basics <- c()
  
  if(!isMultiGroup){
    basics <- c(basics, "Model of class regularizedSEM")
    basics <- c(basics, paste0("Estimator: ", estimator))
  }else{
    basics <- c(basics, "Multi-group model of class regularizedSEM")
    basics <- c(basics, paste0("Estimators in sub-groups: ", paste0(estimator, collapse = ", ")))
  }
  basics <- c(basics, paste0("Penalty: ", penalty))
  basics <- c(basics, paste0("Method: ", method))
  
  basics <- c(basics, paste0("Regularized parameters: ", 
                             paste0(object@regularized, collapse = ", "))
  )
  
  nextSteps <- c("Next steps:")
  nextSteps <- c(nextSteps,
                 c(
                   "Use coef(object) to get the parameter estimates of the model. With coef(object, criterion = 'BIC') parameters estimates at the lowest BIC can be extracted.",
                   "Use plot(object) to plot the parameter estimates of the model.",
                   "Use fitIndices(object) to get the fit indices."
                 )
  )
  
  rlang::inform(basics)
  rlang::inform(nextSteps)
})

#' summary
#' @param object object of class regularizedSEM
#' @param ... not used
#' @return No return value, just prints estimates
#' @export
setMethod("summary", "regularizedSEM", function (object, ...) {
  show(object)
})

#' coef
#' 
#' Returns the parameter estimates of a regularizedSEM
#' 
#' @param object object of class regularizedSEM
#' @param ... criterion can be one of the ones returned by fitIndices. If set to NULL, all parameters will be returned
#' @return parameters of the model as data.frame
#' @export
setMethod("coef", "regularizedSEM", function (object, ...) {
  dotdotdot <- list(...)
  if("criterion" %in% names(dotdotdot)){
    criterion <- dotdotdot$criterion
  }else{
    criterion <- NULL
  }
  
  tuningParameters <- object@parameters[, !colnames(object@parameters) %in% object@parameterLabels,
                                        drop=FALSE] 
  estimates <- as.matrix(object@parameters[,object@parameterLabels,
                                           drop=FALSE])
  
  if(ncol(object@transformations) != 0){
    transformations <- as.matrix(object@transformations[,
                                                        !colnames(object@transformations) %in% colnames(tuningParameters), 
                                                        drop = FALSE])
  }else{
    transformations <- matrix(nrow = 0, ncol = 0)
  }
  
  if(!is.null(criterion)){
    
    fits <- fitIndices(object)
    
    if(!criterion %in% colnames(fits))
      stop("Could not find", criterion, "in fitIndices(object).")
    
    bestFit <- which(fits[,criterion] == min(fits[,criterion]))[1]
    
    coefs <- new("lessSEMCoef")
    coefs@tuningParameters <- tuningParameters[bestFit,,drop = FALSE]
    coefs@estimates <- estimates[bestFit,,drop = FALSE]
    if(ncol(object@transformations) != 0){
      coefs@transformations <- transformations[bestFit,,drop = FALSE]
    }else{
      coefs@transformations <- transformations
    }
    
    return(coefs) 
  }
  
  coefs <- new("lessSEMCoef")
  coefs@tuningParameters <- tuningParameters
  coefs@estimates <- estimates
  coefs@transformations <- transformations
  
  return(coefs)
})

#' AIC
#' 
#' returns the AIC
#' 
#' @param object object of class regularizedSEM
#' @param ... not used
#' @param k multiplier for number of parameters
#' @return AIC values
#' @export
setMethod("AIC", "regularizedSEM", function (object, ..., k = 2) {
  if(object@penalty == "ridge" & !all(object@inputArguments$tuningParameters == 0))
    stop("AIC not supported for this penalty.")
  
  fits <- object@fits
  fits$AIC <- fits$m2LL + k*fits$nonZeroParameters
  
  return(fits)
})

#' BIC
#' 
#' returns the BIC
#' 
#' @param object object of class regularizedSEM
#' @param ... not used
#' @return BIC values
#' @export
setMethod("BIC", "regularizedSEM", function (object, ...) {
  N <- object@internalOptimization$N
  fits <- object@fits
  
  if(object@penalty == "ridge" & !all(object@inputArguments$tuningParameters == 0))
    stop("BIC not supported for this penalty.")
  
  fits <- object@fits
  fits$BIC <- fits$m2LL + log(N)*fits$nonZeroParameters
  
  return(fits)
  
})

#' plots the regularized and unregularized parameters for all levels of lambda
#' 
#' @param x object of class gpRegularized
#' @param y not used
#' @param ... use regularizedOnly=FALSE to plot all parameters
#' @return either an object of ggplot2 or of plotly
#' @export
setMethod("plot", 
          c(x = "regularizedSEM", y = "missing"), 
          function (x, y, ...) {
            if("regularizedOnly" %in% names(list(...))){
              regularizedOnly <- list(...)$regularizedOnly
            }else{
              regularizedOnly <- TRUE
            }
            parameters <- x@parameters
            tuningParameters <- x@parameters[,!colnames(x@parameters)%in%x@parameterLabels,drop=FALSE]
            tuningParameters <- tuningParameters[,apply(tuningParameters,2,function(x) length(unique(x)) > 1),drop=FALSE]
            
            nTuning <- ncol(tuningParameters)
            
            if(nTuning > 2) 
              stop("Plotting currently only supported for up to 2 tuning parameters")
            if(nTuning == 2 & !requireNamespace("plotly", quietly = TRUE))
              stop("Plotting more than one tuning parameter requires the package plotly")
            
            
            if(regularizedOnly){
              
              parameters <- cbind(
                tuningParameters,
                parameters[,x@regularized, drop = FALSE]
              )
              parametersLong <- tidyr::pivot_longer(data = parameters, cols = x@regularized)
              
            }else{
              
              parameters <- cbind(
                tuningParameters,
                parameters
              )
              parametersLong <- tidyr::pivot_longer(data = parameters, cols = x@parameterLabels)
              
            }
            
            if(nTuning == 1){
              
              return(
                ggplot2::ggplot(data = parametersLong,
                                mapping = ggplot2::aes(
                                  x = .data[[colnames(tuningParameters)]], 
                                  y = .data[["value"]],
                                  group = .data[["name"]])) +
                  ggplot2::geom_line(colour = "#008080")+
                  ggplot2::ggtitle("Regularized Parameters")
              )
              
            }else{
              parametersLong$name <- paste0(parametersLong$name, 
                                            "_", 
                                            unlist(parametersLong[,colnames(tuningParameters)[2]]))
              parametersLong$tp1 <- unlist(parametersLong[,colnames(tuningParameters)[1]])
              parametersLong$tp2 <- unlist(parametersLong[,colnames(tuningParameters)[2]])
              plt <- plotly::layout(
                plotly::plot_ly(parametersLong, 
                                x = ~tp1, y = ~tp2, z = ~value, 
                                type = 'scatter3d',
                                mode = 'lines',
                                opacity = 1,
                                color = ~name,
                                split = ~tp2,
                                line = list(width = 6, 
                                            reverscale = FALSE)
                ), 
                scene = list(xaxis = list(title = colnames(tuningParameters)[1]), 
                             yaxis = list(title = colnames(tuningParameters)[2]))
              )
              return(plt)
              
            }
            
          })

#' estimates
#' 
#' @param object object of class regularizedSEM
#' @param criterion fit index (e.g., AIC) used to select the parameters
#' @param transformations boolean: Should transformations be returned?
#' @return returns a matrix with estimates
#' @export
setMethod("estimates", "regularizedSEM", function(object, criterion = NULL, transformations = FALSE) {
  
  if(transformations)
    return(cbind(
      coef(object, criterion = criterion)@estimates,
      coef(object, criterion = criterion)@transformations)
    )
  
  return(coef(object, criterion = criterion)@estimates)
  
})

#' fitIndices
#' 
#' @param object object of class regularizedSEM
#' @return returns a data.frame with fit indices
#' @export
setMethod("fitIndices", "regularizedSEM", function(object) {
  
  fits <- object@fits
  
  usesLikelihood <- any(!is.na(fits$m2LL))
  
  if(usesLikelihood){
    multiGroup <- !is(object@inputArguments$lavaanModel, "lavaan")
    
    if(!multiGroup){
      dataset <- lavInspect(object@inputArguments$lavaanModel, "data")
      # remove empty rows:
      if(any(apply(dataset,1,function(x) all(is.na(x))))){
        warning("Your data set has rows where all observations are missing. lessSEM will",
                "remove those rows, but it is recommended to do so before fitting the models.")
        dataset <- dataset[!apply(dataset,1,function(x) all(is.na(x))),,drop = FALSE]
      }
      sampstats <- lavInspect(object@inputArguments$lavaanModel, "sampstat")
      N <- nrow(dataset)
    }
    
    # fit indices
    fits$AIC <- AIC(object)$AIC
    fits$BIC <- BIC(object)$BIC
    
    # The following variants of the AIC are adapted from here:
    # https://search.r-project.org/CRAN/refmans/AICcmodavg/html/AICc.html
    if(!multiGroup){
      fits$AICc <- fits$m2LL + 2*fits$nonZeroParameters * (N/(N - fits$nonZeroParameters - 1))
      
      # Chi^2
      
      if(is.null(sampstats$mean)){
        sampstats$mean <- apply(dataset, 2, mean, na.rm = TRUE)
        satPar <- nrow(sampstats$cov)*(ncol(sampstats$cov)+1)/2
      }else{
        satPar <- nrow(sampstats$cov)*(ncol(sampstats$cov)+1)/2 + length(sampstats$mean)
      }
      
      saturatedFit <- -2*sum(apply(dataset, 1, function(x) mvtnorm::dmvnorm(x = x[!is.na(x), drop = FALSE], 
                                                                            mean = sampstats$mean[!is.na(x), drop = FALSE], 
                                                                            sigma = sampstats$cov[!is.na(x), !is.na(x), drop = FALSE], 
                                                                            log = TRUE))
      )
      
      fits$chisq <- fits$m2LL - saturatedFit
      fits$df <- satPar - fits$nonZeroParameters
      
      # RMSEA
      # degrees of freedom
      lambda <- fits$chisq - fits$df
      # Note: lavaan uses df*N instead of df*(N-1)!
      N <- nrow(dataset)
      
      fits$rmsea <- 0
      fits$rmsea[lambda >= 0] <- sqrt(lambda[lambda>=0] / (fits$df[lambda>=0] * N))
    }
  }
  
  return(fits)
  
})