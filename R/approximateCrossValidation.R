approximateCrossValidation <- function(lavaanModel, SEM = NULL, k, individualPenaltyFunction = NULL, raw = FALSE, penalty = NULL, ...){
  
  if(is.null(SEM)){
    return(approximateCrossValidationLavaan(lavaanModel = lavaanModel, 
                                            k = k,
                                            individualPenaltyFunction = individualPenaltyFunction, 
                                            raw = raw, 
                                            ... = ...))
  }
  
  if(is(SEM, "Rcpp_SEMCpp")){
    return(approximateCrossValidationRcpp_SEMCpp(SEM = SEM, 
                                                 k = k,
                                                 individualPenaltyFunction = individualPenaltyFunction, 
                                                 raw = raw, 
                                                 ... = ...))
  }
  
  if(!is(lavaanModel, "lavaan")){
    stop("lavaanModel must be of class lavaan")
  }
  
  if(lavaanModel@Options$estimator != "ML") stop("lavaanModel must be fit with ml estimator.")
  
  data <- try(lavaan::lavInspect(lavaanModel, "data"))
  if(is(data, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  
  
  if(is(SEM, "regsem")){
    if(is.null(penalty)){stop("Requires specification of penalty")}
    return(approximateCrossValidationRegsem(SEM = SEM, 
                                            k = k,
                                            individualPenaltyFunction = individualPenaltyFunction, 
                                            ... = ...))
  }
  
  if(is(SEM, "cvregsem")){
    if(is.null(penalty)){stop("Requires specification of penalty")}
    return(approximateCrossValidationCvregsem(cvregsemModel = SEM, 
                                              lavaanModel = lavaanModel,
                                              k = k,
                                              penalty = penalty, 
                                              ... = ...))
  }
  
  if(is(SEM, "lslx")){
    return(approximateCrossValidationLslx(SEM = SEM, 
                                          k = k,
                                          individualPenaltyFunction = individualPenaltyFunction, 
                                          raw = raw, 
                                          ... = ...))
  }
  
}

approximateCrossValidationRcpp_SEMCpp <- function(SEM, k, individualPenaltyFunction = NULL, raw = FALSE, ...){
  if(!is(SEM, "Rcpp_SEMCpp")){
    stop("SEM must be of class Rcpp_SEMCpp")
  }
  
  additionalArguments <- list(...)
  
  parameters <- getParameters(SEM = SEM, raw = raw)
  dataSet <- SEM$rawData
  N <- nrow(dataSet)
  
  # step 1: compute scores
  scores <- getScores(SEM = SEM, raw = raw)
  
  # compute penalty scores
  
  if(!is.null(individualPenaltyFunction)){
    for(i in 1:N){
      
      scores[i,] <- scores[i,] + do.call(what = numDeriv::grad, 
                                         args = c(list(
                                           "func" = individualPenaltyFunction,
                                           "x" = parameters[colnames(scores)]),
                                           additionalArguments
                                         ))
      
    }
  }
  
  hessian <- getHessian(SEM = SEM, raw = raw)
  
  if(!is.null(individualPenaltyFunction)){
    hessian + N*do.call(what = numDeriv::hessian, 
                        args = c(list(
                          "func" = individualPenaltyFunction,
                          "x" = parameters[rownames(hessian)]),
                          additionalArguments
                        ))
  }
  
  # step 2: subgroup-parameters
  if(k < N){
    
    randomCases <- sample(1:N,N)
    subsets <- split(randomCases, sort(randomCases%%k)) # https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks
    
  }else if(k == N){
    
    subsets <- vector("list",N)
    names(subsets) <- 1:N
    for(i in 1:N) subsets[[i]] <- i
    
  }else{
    stop(paste0("k must be <= ", N))
  }
  
  subsetParameters <- vector("list",k)
  leaveOutFits <- rep(0, k)
  names(leaveOutFits) <- paste0("sample", 1:k)
  
  inverseHessian <- (N/(N-1))*solve(hessian)
  for(s in 1:k){
    
    direction <- -t(inverseHessian%*%apply(scores[-subsets[[s]],,drop=FALSE],2,sum))
    
    # Taking a step size of 1 can result in nonsensical parameters. Therefore, a rudimentary line search is used here:
    stepLength <- 1
    steps <- 1
    while(steps<100){
      parameters_s <- parameters + stepLength*direction[,names(parameters)]
      SEM <- setParameters(SEM = SEM, labels = names(parameters), values = parameters_s, raw = raw)
      SEM <- try(fit(SEM), silent = TRUE)
      if(any(class(SEM) == "try-error")) {
        stepLength <- stepLength*.9
        steps <- steps + 1
        next
      }
      # check positive definiteness
      tryChol <- try(chol(SEM$impliedCovariance), silent = TRUE)
      if(any(class(tryChol) == "try-error")) {
        stepLength <- stepLength*.9
        steps <- steps + 1
        next
      }
      
      break
    }
    
    if(steps != 1) print(paste0("Used linesearch with ", steps, " steps. Final stepLength: ", stepLength))
    
    names(parameters_s) <- names(parameters)
    subsetParameters[[s]] <- parameters_s
    
    for(i in 1:length(subsets[[s]])){
      leaveOutFits[s] <- leaveOutFits[s] + individualMinus2LogLikelihood(par = subsetParameters[[s]], 
                                                                         SEM = SEM, 
                                                                         data = dataSet[subsets[[s]][i],], 
                                                                         raw = raw)
    }
    
  }
  
  return(list("leaveOutFits" = leaveOutFits,
              "subsets" = subsets,
              "scores" = scores,
              "subsetParameters" = subsetParameters
  ))
  
}


approximateCrossValidationLavaan <- function(lavaanModel, k, individualPenaltyFunction = NULL, raw = FALSE, ...){
  if(!is(lavaanModel, "lavaan")){
    stop("lavaanModel must be of class lavaan")
  }
  
  if(lavaanModel@Options$estimator != "ML") stop("lavaanModel must be fit with ml estimator.")
  
  data <- try(lavaan::lavInspect(lavaanModel, "data"))
  if(is(data, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  
  aCVSEM <- SEMFromLavaan(lavaanModel = lavaanModel, rawData = data, transformVariances = TRUE)
  
  return(approximateCrossValidationRcpp_SEMCpp(SEM = aCVSEM, 
                                               k = k,
                                               individualPenaltyFunction = individualPenaltyFunction, 
                                               raw = raw, 
                                               ... = ...))
}

approximateCrossValidationRegsem <- function(regsemModel, lavaanModel, k, individualPenaltyFunction = NULL, raw = FALSE, ...){
  if(!is(lavaanModel, "lavaan")){
    stop("lavaanModel must be of class lavaan")
  }
  if(is.null(individualPenaltyFunction)) stop("approximateCrossValidationRegsem requires the individualPenaltyFunction to be specified.")
  
  if(lavaanModel@Options$estimator != "ML") stop("lavaanModel must be fit with ml estimator.")
  
  data <- try(lavaan::lavInspect(lavaanModel, "data"))
  if(is(data, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  if(!is(regsemModel, "regsem")){stop("regsemModel must be of class regsem.")}
  
  aCVSEM <- SEMFromLavaan(lavaanModel = lavaanModel, rawData = data, transformVariances = TRUE)
  
  regsemParameters <- regsem2LavaanParameters(regsemModel = regsemModel, lavaanModel = lavaanModel)
  regularizedParameterLabels <- names(regsemParameters)[regsemModel$pars_pen]
  
  message("aCV4SEM assumes that the following parameters were regularized: ", 
          paste0(regularizedParameterLabels, sep = ", "), 
          ". Please make sure that this is correct.")
  
  stop("MISSING IMPLEMENTATION")
  return(approximateCrossValidationRcpp_SEMCpp(SEM = aCVSEM, 
                                               k = k,
                                               individualPenaltyFunction = individualPenaltyFunction, 
                                               raw = raw, 
                                               ... = ...))
}

approximateCrossValidationCvregsem <- function(cvregsemModel, lavaanModel, k, penalty, eps = 1e-4){
  if(!is(lavaanModel, "lavaan")){
    stop("lavaanModel must be of class lavaan")
  }
  if(!penalty %in% c("lasso", "ridge", "elasticNet")) stop("approximateCrossValidationCvregsem currently only supports lasso, ridge, or elastic net as penalty functions")
  
  if(lavaanModel@Options$estimator != "ML") stop("lavaanModel must be fit with ml estimator.")
  
  data <- try(lavaan::lavInspect(lavaanModel, "data"))
  if(is(data, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  if(!is(cvregsemModel, "cvregsem")) stop("cvregsemModel must be of class regsem.")
  
  aCVSEM <- SEMFromLavaan(lavaanModel = lavaanModel, rawData = data, transformVariances = TRUE)
  
  regsemParameters <- cvregsem2LavaanParameters(cvregsemModel = cvregsemModel, lavaanModel = lavaanModel)
  regularizedParameterLabels <- colnames(regsemParameters)[cvregsemModel$pars_pen]
  
  message("aCV4SEM assumes that the following parameters were regularized: ", 
          paste0(regularizedParameterLabels, collapse = ", "), 
          ". Please make sure that this is correct.")
  
  # extract tuning parameters
  lambdas <- cvregsemModel$fits[,"lambda"]
  if(penalty == "elasticNet") {
    if(is.null(cvregsemModel$call$alpha)) stop("Cannot find the tuning parameter alpha in cvregsemModel. Please make sure to call cv_regsem with alpha explicitly specified")
    alphas <- rep(cvregsemModel$call$alpha, length(lambdas)) # cv_regsem cannot use a grid of alphas
  }
  
  if(penalty == "elasticNet") {
    coln <- paste0("lambda=",lambdas, "; alpha=", alphas)
  }else{
    coln <- paste0("lambda=",lambdas)
  }
  aCVs <- matrix(NA,
                 nrow = k, 
                 ncol = nrow(regsemParameters),
                 dimnames = list(
                   paste0("sample", 1:k),
                   coln
                 ))
  
  for(p in 1:nrow(regsemParameters)){
    if(anyNA(regsemParameters[p,])) next
    aCVSEM <- setParameters(SEM = aCVSEM,
                            labels = colnames(regsemParameters),
                            values = regsemParameters[p,],
                            raw = FALSE)
    aCVSEM <- fit(aCVSEM)
    
    if(penalty == "elasticNet"){
      aCV <- approximateCrossValidationRcpp_SEMCpp(SEM = aCVSEM, 
                                                   k = k,
                                                   individualPenaltyFunction = smoothElasticNet, 
                                                   raw = FALSE, 
                                                   regularizedParameterLabels = regularizedParameterLabels,
                                                   lambda = lambdas[p],
                                                   alpha = alphas[p],
                                                   eps = eps)
      
    }
    if(penalty == "ridge"){
      aCV <- approximateCrossValidationRcpp_SEMCpp(SEM = aCVSEM, 
                                                   k = k,
                                                   individualPenaltyFunction = ridge, 
                                                   raw = FALSE, 
                                                   regularizedParameterLabels = regularizedParameterLabels,
                                                   lambda = lambdas[p])
      
    }
    if(penalty == "lasso"){
      aCV <- approximateCrossValidationRcpp_SEMCpp(SEM = aCVSEM, 
                                                   k = k,
                                                   individualPenaltyFunction = smoothLASSO, 
                                                   raw = FALSE, 
                                                   regularizedParameterLabels = regularizedParameterLabels,
                                                   lambda = lambdas[p],
                                                   eps = eps)
    }
    
    aCVs[paste0("sample", 1:k),p] <- aCV$leaveOutFits
  }
  
  return(list(
    "lambda" = lambdas,
    "regsemParameters" = regsemParameters,
    "approximateCV" = aCVs
  )
  )
}
