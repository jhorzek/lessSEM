getParameters.lavaan <- function(model, removeDuplicates = TRUE){
  parameters <- coef(model)
  if(!removeDuplicates) return(parameters)
  return(parameters[unique(names(parameters))])
}

# setParameters.lavaan <- function(par, model){
#   
#   parameters <- getParameters.lavaan(model, removeDuplicates = FALSE)
#   
#   if(any(!names(par) %in% names(parameters))) stop("Unkown parameters passed to setParameters.lavaan")
#   
#   for(parlabel in names(par)) parameters[names(parameters) == parlabel] <- par[parlabel]
#   
#   # extract some model contents for convenience
#   lavmodel = model@Model
#   
#   # change parameters; lavaan expects a vector with duplicates for parameters with constraints
#   lavmodel <- lav_model_set_parameters(lavmodel, parameters)
#   
#   # also: update the model object
#   model@Model <- lavmodel
#   model@implied <- lavaan:::lav_model_implied(lavmodel = lavmodel)
# 
#   return(model)
# }