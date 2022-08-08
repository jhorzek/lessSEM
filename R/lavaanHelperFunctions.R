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
#' @export
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
#' @export
.cvregsem2LavaanParameters <- function(cvregsemModel, lavaanModel){
  if(!is(cvregsemModel, "cvregsem")) stop("cvregsemModel must be of class cvregsem.")
  if(!is(lavaanModel, "lavaan")) stop("lavaanModel must be of class lavaan.")
  parameters <- cvregsemModel$parameters
  
  parLabels <- lessSEM:::.lavaan2regsemLabels(lavaanModel = lavaanModel)
  
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
  
  parLabels <- lessSEM:::.lavaan2regsemLabels(lavaanModel = lavaanModel)
  
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