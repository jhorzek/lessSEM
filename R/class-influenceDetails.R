# Note: we define a custom logLik - Function because the generic one is 
# using df = number of parameters which might be confusing.
setClass("approximateInfluence",
         representation = representation(
           subsets="matrix",
           tuningParameters = "data.frame",
           subsetParameters="array",
           subsetFits = "array",
           parameterLabels = "character",
           regularized  = "character"
         ))

#' AIC
#' 
#' returns the AIC
#' 
#' @param object object of class approximateInfluence
#' @export
setMethod("AIC", "approximateInfluence", function (object) {
  AICs <- as.data.frame(matrix(NA, 
                               nrow = nrow(object@tuningParameters)*ncol(object@subsets), 
                               ncol = ncol(object@tuningParameters)+3,
                               dimnames = list(NULL, c(colnames(object@tuningParameters), "subset", "nonZeroParameters", "AIC")))
  )
  AICs$lambda <- rep(object@tuningParameters$lambda, ncol(object@subsets))
  AICs$alpha <- rep(object@tuningParameters$alpha, ncol(object@subsets))
  AICs$subset <- rep(1:ncol(object@subsets), each = nrow(object@tuningParameters))
  notZeroed <- object@subsetParameters != 0
  
  for(i in 1:nrow(object@tuningParameters)){
    select <- AICs$lambda == object@tuningParameters$lambda[i] & AICs$alpha == object@tuningParameters$alpha[i]
    AICs$nonZeroParameters[select] <- apply(notZeroed[,,i],1,sum)
    AICs$AIC[select] <- object@subsetFits[,,i] + 2*AICs$nonZeroParameters[select]
  }
  
  return(AICs)
})