#' getLavaanParameters
#' 
#' helper function: returns a labeled vector with parameters from lavaan
#'
#' @param lavaanModel model of class lavaan
#' @param removeDuplicates should duplicated parameters be removed?
#' @returns returns a labeled vector with parameters from lavaan
#' @export
getLavaanParameters <- function(lavaanModel, removeDuplicates = TRUE){
  if(!is(lavaanModel, "lavaan")) stop("lavaanModel must be of class lavaan.")
  parameters <- coef(lavaanModel)
  if(!removeDuplicates) return(parameters)
  return(parameters[unique(names(parameters))])
}

#' .lavaan2regsemLabels
#' 
#' helper function: regsem and lavaan use slightly different parameter labels. This function
#' can be used to get both sets of labels.
#'
#' @param lavaanModel model of class lavaan
#' @returns a list with lavaan and regsem labels
#' @keywords internal
.lavaan2regsemLabels <- function(lavaanModel){
  if(!is(lavaanModel, "lavaan")) stop("lavaanModel must be of class lavaan.")
  # extract parameters
  parameterIDs <- lavaanModel@ParTable$id[lavaanModel@ParTable$free != 0]
  ops <- lavaanModel@ParTable$op
  lavaanLabels <- paste0(lavaanModel@ParTable$lhs, ops, lavaanModel@ParTable$rhs)
  lavaanLabels[lavaanModel@ParTable$label!=""] <- lavaanModel@ParTable$label[lavaanModel@ParTable$label!=""]
  lavaanLabels <- lavaanLabels[lavaanModel@ParTable$free != 0]
  
  opsRegsem <- lavaanModel@ParTable$op
  opsRegsem[opsRegsem == "=~"] <- "->"
  regsemLabels <- paste(lavaanModel@ParTable$lhs, opsRegsem, lavaanModel@ParTable$rhs)
  regsemLabels[opsRegsem == "~1"] <- paste("1 ->", lavaanModel@ParTable$lhs[opsRegsem == "~1"])
  #regsemLabels[lavaanModel@ParTable$label!=""] <- lavaanModel@ParTable$label[lavaanModel@ParTable$label!=""]
  regsemLabels <- regsemLabels[lavaanModel@ParTable$free != 0]
  
  return( 
    list("lavaanLabels" = lavaanLabels,
         "regsemLabels" = regsemLabels)
  )
}

#' .cvregsem2LavaanParameters
#' 
#' helper function: regsem and lavaan use slightly different parameter labels. This function
#' can be used to translate the parameter labels of a cv_regsem object to lavaan labels
#' 
#' @param cvregsemModel model of class cvregsem
#' @param lavaanModel model of class lavaan
#' @returns regsem parameters with lavaan labels
#' @keywords internal
.cvregsem2LavaanParameters <- function(cvregsemModel, lavaanModel){
  if(!is(cvregsemModel, "cvregsem")) stop("cvregsemModel must be of class cvregsem.")
  if(!is(lavaanModel, "lavaan")) stop("lavaanModel must be of class lavaan.")
  parameters <- cvregsemModel$parameters
  
  parLabels <- .lavaan2regsemLabels(lavaanModel = lavaanModel)
  
  # resort parameters
  parameters <- parameters[,parLabels$regsemLabels]
  
  colnames(parameters) <- parLabels$lavaanLabels
  
  return(parameters)
}

#' regsem2LavaanParameters
#' 
#' helper function: regsem and lavaan use slightly different parameter labels. This function
#' can be used to translate the parameter labels of a cv_regsem object to lavaan labels
#' 
#' @param regsemModel model of class regsem
#' @param lavaanModel model of class lavaan
#' @returns regsem parameters with lavaan labels
#' @export
regsem2LavaanParameters <- function(regsemModel, lavaanModel){
  if(!is(regsemModel, "regsem")) stop("regsemModel must be of class cvregsem.")
  if(!is(lavaanModel, "lavaan")) stop("lavaanModel must be of class lavaan.")
  
  parameters <- unlist(regsemModel$out$pars)
  
  parLabels <- .lavaan2regsemLabels(lavaanModel = lavaanModel)
  
  # resort parameters
  parameters <- parameters[parLabels$regsemLabels]
  
  names(parameters) <- parLabels$lavaanLabels
  
  return(parameters)
}

#' lavaan2lslxLabels
#' 
#' helper function: lslx and lavaan use slightly different parameter labels. This function
#' can be used to get both sets of labels.
#'
#' @param lavaanModel model of class lavaan
#' @returns list with lavaan labels and lslx labels
#' @export
lavaan2lslxLabels <- function(lavaanModel){
  if(!is(lavaanModel, "lavaan")) stop("lavaanModel must be of class lavaan.")
  # extract parameters
  parameterIDs <- lavaanModel@ParTable$id[lavaanModel@ParTable$free != 0]
  ops <- lavaanModel@ParTable$op
  lavaanLabels <- paste0(lavaanModel@ParTable$lhs, ops, lavaanModel@ParTable$rhs)
  lavaanLabels[lavaanModel@ParTable$label!=""] <- lavaanModel@ParTable$label[lavaanModel@ParTable$label!=""]
  lavaanLabels <- lavaanLabels[lavaanModel@ParTable$free != 0]
  
  parFrameLslx <- data.frame(lhs = lavaanModel@ParTable$rhs, ops = ops, rhs = lavaanModel@ParTable$lhs, g = "/g")
  parFrameLslx <- parFrameLslx[lavaanModel@ParTable$free != 0,]
  parFrameLslx$ops[parFrameLslx$ops == "=~"] <- "<-"
  parFrameLslx$ops[parFrameLslx$ops == "~~"] <- "<->"
  parFrameLslx$lhs[parFrameLslx$ops == "~1"] <- parFrameLslx$rhs[parFrameLslx$ops == "~1"]
  parFrameLslx$rhs[parFrameLslx$ops == "~1"] <- ""
  parFrameLslx$ops[parFrameLslx$ops == "~1"] <- "<-1"
  
  lslxLabels <- apply(parFrameLslx,1,paste0, collapse = "")
  names(lslxLabels) <- NULL
  
  return( 
    list("lavaanLabels" = lavaanLabels,
         "lslxLabels" = lslxLabels)
  )
}

#' lessSEM2Lavaan
#' 
#' Creates a lavaan model object from lessSEM (only if possible).
#' 
#' @param regularizedSEM object created with lessSEM
#' @param lambda value for tuning parameter lambda
#' @param theta value for tuning parameter theta
#' @return lavaan model
lessSEM2Lavaan <- function(regularizedSEM, lambda, theta = NULL){
  if("theta" %in% colnames(regularizedSEM@fits) &
     is.null(theta))
    stop("Your model uses tuning parameter theta, but no theta value was specified")
  
  if(!any(regularizedSEM@fits$lambda == lambda))
    stop("Could not find the specified lambda in your model.")
  
  if(!is.null(theta) &&
     !any(regularizedSEM@fits$theta == theta))
    stop("Could not find the specified theta in your model.")
  
  # extract coefficients and remove tuning parameters:
  whichRow <- regularizedSEM@fits$lambda == lambda
  if(!is.null(theta)){
    whichRow <- whichRow & regularizedSEM@fits$theta == theta
  }
  
  if(sum(whichRow) != 1) 
    stop("Error while selecting parameters: Instead of returning parameters for a single model, multiple model parameters have been returned")
  
  # extract lesssem parameters
  if(ncol(regularizedSEM@transformations) == 0){
    # model without transformations
    lessSEMEstimates <- regularizedSEM@parameters
  }else{
    lessSEMEstimates <- cbind(regularizedSEM@parameters, regularizedSEM@transformations)
  }
  
  # check if it is a multi-group model
  if(is(object = regularizedSEM@inputArguments$lavaanModel, class2 = "list")){
    nModels <- length(regularizedSEM@inputArguments$lavaanModel)
    lavaanModels <- vector("list", nModels)
    
    for(m in 1:length(lavaanModels)){
      # We only need those parameters that are also in the lavaan model:
      lavaanModel <- regularizedSEM@inputArguments$lavaanModel[[m]]
      expectedParameters <- names(getLavaanParameters(lavaanModel))
      
      lessSEMEstimates <- unlist(cbind(regularizedSEM@parameters, regularizedSEM@transformations)[whichRow,expectedParameters])
      
      # Now we can change the parameters of the lavaan model to match ours
      lavaanParTable <- lavaan::parTable(object = lavaanModel)
      lavaanParTable$se <- NA
      
      for(i in 1:nrow(lavaanParTable)){
        if(! "label" %in% colnames(lavaanParTable) ||
           lavaanParTable$label[i] == ""){
          label <- paste0(lavaanParTable$lhs[i], lavaanParTable$op[i], lavaanParTable$rhs[i])
        }else{
          label <- lavaanParTable$label[i]
        }
        
        if(label %in% names(lessSEMEstimates)){
          lavaanParTable$est[i] <- lessSEMEstimates[label]
        }
      }
      
      lavaanModels[[m]] <- suppressWarnings(lavaan::sem(model = lavaanParTable,
                                                        data = lavInspect(lavaanModel, "data"),
                                                        do.fit= FALSE))
      
      
    }
    return(lavaanModels)
    
  }
  
  # We only need those parameters that are also in the lavaan model:
  lavaanModel <- regularizedSEM@inputArguments$lavaanModel
  expectedParameters <- names(getLavaanParameters(lavaanModel))
  
  lessSEMEstimates <- unlist(lessSEMEstimates[whichRow,expectedParameters])
  
  # Now we can change the parameters of the lavaan model to match ours
  lavaanParTable <- lavaan::parTable(object = lavaanModel)
  lavaanParTable$se <- NA
  
  for(i in 1:nrow(lavaanParTable)){
    if(! "label" %in% colnames(lavaanParTable) ||
       lavaanParTable$label[i] == ""){
      label <- paste0(lavaanParTable$lhs[i], lavaanParTable$op[i], lavaanParTable$rhs[i])
    }else{
      label <- lavaanParTable$label[i]
    }
    
    if(label %in% names(lessSEMEstimates)){
      lavaanParTable$est[i] <- lessSEMEstimates[label]
    }
  }
  
  updatedLavaanModel <- suppressWarnings(lavaan::sem(model = lavaanParTable,
                                                     data = lavInspect(lavaanModel, "data"),
                                                     do.fit= FALSE))
  return(updatedLavaanModel)
}