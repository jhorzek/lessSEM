smoothACVRcpp_SEMCpp <- function(SEM, k, 
                                 individualPenaltyFunction = NULL, 
                                 individualPenaltyFunctionGradient = NULL, 
                                 individualPenaltyFunctionHessian = NULL, 
                                 raw = FALSE, 
                                 penaltyFunctionArguments = NULL){
  if(!is.null(individualPenaltyFunction)){
    if(is.null(individualPenaltyFunctionGradient) || is.null(individualPenaltyFunctionHessian)){
      warning("You did not specify the individualPenaltyFunctionGradient and individualPenaltyFunctionHessian. We highly recommend that you do so for smooth approximations of non-differential penalties (lasso, adaptive lasso, etc.) as it may improve the results considerably")
    }
  }
  
  if(!is(SEM, "Rcpp_SEMCpp")){
    stop("SEM must be of class Rcpp_SEMCpp")
  }
  
  parameters <- getParameters(SEM = SEM, raw = raw)
  dataSet <- SEM$rawData
  N <- nrow(dataSet)
  
  # step 1: compute scores
  scores <- getScores(SEM = SEM, raw = raw)
  
  # compute penalty scores
  
  if(!is.null(individualPenaltyFunction)){
    for(i in 1:N){
      
      if(is.null(individualPenaltyFunctionGradient)){
        scores[i,] <- scores[i,] + do.call(what = numDeriv::grad, 
                                           args = c(list(
                                             "func" = individualPenaltyFunction,
                                             "x" = parameters[colnames(scores)],
                                             "method" = "simple",
                                             "method.args" = list(eps = 1e-7),
                                             "penaltyFunctionArguments" = penaltyFunctionArguments
                                           )
                                           ))
      }else{
        scores[i,] <- scores[i,] + individualPenaltyFunctionGradient(parameters = parameters[colnames(scores)],
                                                                     penaltyFunctionArguments = penaltyFunctionArguments)
        
      }
      
    }
  }
  
  hessian <- getHessian(SEM = SEM, raw = raw)
  
  if(!is.null(individualPenaltyFunction)){
    
    if(is.null(individualPenaltyFunctionHessian)){
      hessian + N*do.call(what = numDeriv::hessian, 
                          args = list(
                            "func" = individualPenaltyFunction,
                            "x" = parameters[rownames(hessian)],
                            "method.args" = list(eps = 1e-7),
                            "penaltyFunctionArguments" = penaltyFunctionArguments)
      )
    }else{
      for(i in 1:N){
        hessian <- hessian + individualPenaltyFunctionHessian(
          parameters = parameters[colnames(scores)],
          penaltyFunctionArguments = penaltyFunctionArguments
        )
      }
    }
    
    
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

GLMNETACVRcpp_SEMCpp <- function(SEM, 
                                 k, 
                                 penalty, 
                                 raw = FALSE, 
                                 penaltyFunctionArguments = NULL,
                                 maxIter = 100,
                                 epsBreak = 1e-10){
  
  if(!is(SEM, "Rcpp_SEMCpp")){
    stop("SEM must be of class Rcpp_SEMCpp")
  }
  
  if(!penalty %in% c("lasso", "adaptiveLasso")) stop("Currently only lasso and adaptive lasso supported!")
  
  parameters <- getParameters(SEM = SEM, raw = raw)
  dataSet <- SEM$rawData
  N <- nrow(dataSet)
  
  regularizedParameterLabels <- penaltyFunctionArguments$regularizedParameterLabels
  if(is.null(regularizedParameterLabels)) stop("penaltyFunctionArguments$regularizedParameterLabels missing")
  lambda <- penaltyFunctionArguments$lambda
  if(is.null(lambda)) stop("penaltyFunctionArguments$lambda missing")
  
  adaptiveLassoWeights <- penaltyFunctionArguments$adaptiveLassoWeights
  if(is.null(adaptiveLassoWeights) && penalty == "adaptiveLasso") stop("penaltyFunctionArguments$adaptiveLassoWeights missing")
  if(is.null(adaptiveLassoWeights)) adaptiveLassoWeights <- rep(1, length(parameters))
  
  # build subgroups
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
  
  # compute derivatives of -2log-Likelihood without penalty
  
  scores <- getScores(SEM = SEM, raw = raw)
  hessian <- getHessian(SEM = SEM, raw = raw)
  
  # Now do the GLMNET inner iteration for each sub-group
  stepdirections <- matrix(NA, nrow = k, ncol = length(parameters))
  colnames(stepdirections) <- names(parameters)
  
  for(s in 1:k){
    subGroupGradient <- apply(scores[-subsets[[s]],],2,sum) # gradients of all inidividuals but the ones in the subgroup
    subGroupHessian <- ((N-length(subsets[[s]]))/N)*hessian # approximated hessian for all but the subgroup
    subGroupLambda <- (N-length(subsets[[s]]))*lambda
    
    direction <- innerGLMNET(parameters = parameters, 
                             subGroupGradient = subGroupGradient, 
                             subGroupHessian = subGroupHessian, 
                             subGroupLambda = subGroupLambda, 
                             regularized = names(parameters)%in%regularizedParameterLabels, 
                             adaptiveLassoWeights = adaptiveLassoWeights, 
                             maxIter = maxIter, 
                             epsBreak = epsBreak)
    
    rownames(direction) = names(parameters)
    
    # return step direction
    stepdirections[s,names(parameters)] <- direction[names(parameters),]
  }
  
  subsetParameters <- vector("list",k)
  leaveOutFits <- rep(0, k)
  names(leaveOutFits) <- paste0("sample", 1:k)
  
  for(s in 1:k){
    # currently, we do not use line search as this would result in additional computational overhead.
    # therefore a fixed stepLength of 1 is used; this could be adapted to get better results.
    stepLength <- 1
    parameters_s <- parameters + 1*stepdirections[s,names(parameters)]
    
    SEM <- setParameters(SEM = SEM, labels = names(parameters), values = parameters_s, raw = raw)
    SEM <- try(fit(SEM), silent = TRUE)
    subsetParameters[[s]] <- getParameters(SEM = SEM, raw = FALSE)
    
    for(i in 1:length(subsets[[s]])){
      leaveOutFits[s] <- leaveOutFits[s] + individualMinus2LogLikelihood(par = subsetParameters[[s]], 
                                                                         SEM = SEM, 
                                                                         data = dataSet[subsets[[s]][i],], 
                                                                         raw = FALSE)
    }
    
  }
  
  return(
    list("leaveOutFits" = leaveOutFits,
         "subsets" = subsets,
         "scores" = scores,
         "subsetParameters" = subsetParameters)
  )
  
}