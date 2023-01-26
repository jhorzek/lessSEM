regularize <- function(lavaanModel,
                       method = "ista", 
                       modifyModel = lessSEM::modifyModel(),
                       control = lessSEM::controlIsta()){
  regularizedSEMBase <- list(
    lavaanModel = lavaanModel,
    method = method,
    modifyModel = modifyModel,
    control = control,
    penalties = list()
  )
  
  class(regularizedSEMBase) <- "regularizedSEMBase"
  
  return(regularizedSEMBase)
}

addLasso <- function(regularizedSEMBase,
                     regularized,
                     lambdas){
  
  if(!is(regularizedSEMBase, "regularizedSEMBase"))
    stop("regularizedSEMBase must be of class regularizedSEMBase. ",
         "These models can be created with the regularize() function.")
  
  tps <- expand.grid(lambda = lambdas,
                     theta = 0)
  
  regularizedSEMBase$penalties <- c(
    regularizedSEMBase$penalties,
    list(regularized = regularized, 
         tps = tps,
         type = "lasso")
  )
  
  return(regularizedSEMBase)
}

addScad <- function(regularizedSEMBase,
                    regularized,
                    lambdas,
                    thetas){
  
  if(!is(regularizedSEMBase, "regularizedSEMBase"))
    stop("regularizedSEMBase must be of class regularizedSEMBase. ",
         "These models can be created with the regularize() function.")
  
  tps <- expand.grid(lambda = lambdas,
                     theta = thetas)
  
  regularizedSEMBase$penalties <- c(
    regularizedSEMBase$penalties,
    list(regularized = regularized, 
         tps = tps,
         type = "scad")
  )
  
  return(regularizedSEMBase)
}

fit <- function(regularizedSEMBase){
  
  if(length(regularizedSEMBase$penalties) == 0)
    stop("Penaties are missing")
  
}