#' getLavaanParameters
#' 
#' helper function: returns a labeled vector with parameters from lavaan
#'
#' @param lavaanModel model of class lavaan
#' @param removeDuplicates should duplicated parameters be removed?
#' @returns returns a labeled vector with parameters from lavaan
#' @examples 
#' library(lessSEM)
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
#' getLavaanParameters(lavaanModel)                       
#' @export
getLavaanParameters <- function(lavaanModel, removeDuplicates = TRUE){
  if(!is(lavaanModel, "lavaan")) stop("lavaanModel must be of class lavaan.")
  parameters <- lavaan::coef(lavaanModel)
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
#' @examples
#' ## The following is adapted from ?regsem::regsem.
#' #library(lessSEM)
#' #library(regsem)
#' ## put variables on same scale for regsem
#' #HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
#' #
#' #mod <- '
#' #f =~ 1*x1 + l1*x2 + l2*x3 + l3*x4 + l4*x5 + l5*x6 + l6*x7 + l7*x8 + l8*x9
#' #'
#' ## Recommended to specify meanstructure in lavaan
#' #lavaanModel <- cfa(mod, HS, meanstructure=TRUE)
#' #
#' #regsemModel <- regsem(lavaanModel, 
#' #                lambda = 0.3, 
#' #                gradFun = "ram",
#' #                type="lasso",
#' #                pars_pen=c("l1", "l2", "l6", "l7", "l8"))
#' # regsem2LavaanParameters(regsemModel = regsemModel,
#' #                         lavaanModel = lavaanModel)
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
#' lavaan2lslxLabels(lavaanModel)
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
#' Creates a lavaan model object from lessSEM (only if possible). Pass either
#' a criterion or a combination of lambda, alpha, and theta.
#' 
#' @param regularizedSEM object created with lessSEM
#' @param criterion criterion used for model selection. Currently supported are
#' "AIC" or "BIC"
#' @param lambda value for tuning parameter lambda
#' @param alpha value for tuning parameter alpha
#' @param theta value for tuning parameter theta
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
#' regularized <- lasso(lavaanModel,
#'                      regularized = paste0("l", 11:15), 
#'                      lambdas = seq(0,1,.1))
#' 
#' # using criterion
#' lessSEM2Lavaan(regularizedSEM = regularized, 
#'                criterion = "AIC")
#'                
#' # using tuning parameters (note: we only have to specify the tuning
#' # parameters that are actually used by the penalty function. In case
#' # of lasso, this is lambda):
#' lessSEM2Lavaan(regularizedSEM = regularized, 
#'                lambda = 1)
#' @return lavaan model
#' @export
lessSEM2Lavaan <- function(regularizedSEM, criterion = NULL, lambda = NULL, alpha = NULL, theta = NULL){
  
  if(is.null(criterion) & is.null(lambda) & is.null(alpha) & is.null(theta))
    stop("Either specify a criterion or lambda, alpha, and theta.")
  
  if(!is.null(criterion)){
    
    bestEstimates <- coef(regularizedSEM, criterion = criterion)
    lambda <- bestEstimates@tuningParameters$lambda
    alpha <- bestEstimates@tuningParameters$alpha
    theta <- bestEstimates@tuningParameters$theta
    
  }
  
  if("theta" %in% colnames(regularizedSEM@fits) &
     is.null(theta)){
    if(length(unique(regularizedSEM@fits$theta)) != 1)
      stop("Your model uses tuning parameter theta, but no theta value was specified")
  }
  
  if("alpha" %in% colnames(regularizedSEM@fits) &
     is.null(theta)){
    if(length(unique(regularizedSEM@fits$alpha)) != 1)
      stop("Your model uses tuning parameter alpha, but no alpha value was specified")
  }
  
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
                                                        data = lavaan::lavInspect(lavaanModel, "data"),
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
                                                     data = lavaan::lavInspect(lavaanModel, "data"),
                                                     do.fit= FALSE))
  return(updatedLavaanModel)
}

#' .updateLavaan
#' 
#' updates a lavaan model. lavaan has an update function that does exactly that,
#' but it seems to not work with testthat. This is an attempt to hack around the
#' issue...
#' 
#' @param lavaanModel fitted lavaan model
#' @param key label of the element that should be updated
#' @param value new value for the updated element
#' @return lavaan model 
.updateLavaan <- function(lavaanModel, key, value){
  
  callArgs <- names(lavaanModel@call)[names(lavaanModel@call) != ""]
  reconstructCall <- list(
    model = parTable(lavaanModel),
    data = lavInspect(lavaanModel, "data")
  )
  for(i in callArgs){
    if(i %in% names(lavaanModel@Options))
      reconstructCall <- c(reconstructCall,
                           lavaanModel@Options[i])
  }
  
  reconstructCall[[key]] <- value
  
  return(do.call("sem", reconstructCall))
  
}