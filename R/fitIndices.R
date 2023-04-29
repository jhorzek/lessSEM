#' S4 method to compute fit indices (e.g., AIC, BIC, ...)
#'
#' @param object a model fitted with lessSEM
#' @return returns a data.frame with fit indices
#' @export
setGeneric("fitIndices", function(object) {
  standardGeneric("fitIndices")
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
      
      saturatedFit <- -2*sum(apply(dataset, 1, function(x) mvtnorm::dmvnorm(x = x[!is.na(x)], 
                                                                            mean = sampstats$mean[!is.na(x)], 
                                                                            sigma = sampstats$cov[!is.na(x), !is.na(x)], 
                                                                            log = TRUE))
      )
      
      fits$chisq <- fits$m2LL - saturatedFit
      fits$df <- satPar - fits$nonZeroParameters
      
      # RMSEA
      # degrees of freedom
      lambda <- fits$chisq - fits$df
      # Note: lavaan uses df*N instead of df*(N-1)!
      N <- nrow(dataset)
      
      fits$rmsea <- sqrt(lambda / (fits$df * N))
    }
  }
  
  return(fits)
  
})

#' fitIndices
#' 
#' @param object object of class regularizedSEMMixedPenalty
#' @return returns a data.frame with fit indices
#' @export
setMethod("fitIndices", "regularizedSEMMixedPenalty", function(object) {
  
  fits <- object@fits
  
  usesLikelihood <- any(!is.na(fits$m2LL))
  
  if(usesLikelihood){
    multiGroup <- !is(object@inputArguments$lavaanModel, "lavaan")
    
    if(!multiGroup){
      dataset <- lavInspect(object@inputArguments$lavaanModel, "data")
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
      
      saturatedFit <- -2*sum(apply(dataset, 1, function(x) mvtnorm::dmvnorm(x = x[!is.na(x)], 
                                                                            mean = sampstats$mean[!is.na(x)], 
                                                                            sigma = sampstats$cov[!is.na(x), !is.na(x)], 
                                                                            log = TRUE))
      )
      
      fits$chisq <- fits$m2LL - saturatedFit
      fits$df <- satPar - fits$nonZeroParameters
      
      # RMSEA
      # degrees of freedom
      lambda <- fits$chisq - fits$df
      # Note: lavaan uses df*N instead of df*(N-1)!
      N <- nrow(dataset)
      
      fits$rmsea <- sqrt(lambda / (fits$df * N))
    }
  }
  
  return(fits)
  
})

#' fitIndices
#' 
#' @param object object of class cvRegularizedSEM
#' @return returns a data.frame with fit indices
#' @export
setMethod("fitIndices", "cvRegularizedSEM", function(object) {
  # In case of cross-validated models, we do not compute any fit measures
  # but just return the cv-fit
  return(object@cvfits)
})