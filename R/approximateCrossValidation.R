#' GLMNETACVRcpp_SEMCpp
#' 
#' internal function for approximate cross-validation based on the internal model representation of aCV4SEM. The GLMNET part refers to the fact
#' that this function uses the GLMNET optimizer for non-differentiable penalty functions. 
#' 
#' @param SEM model of class Rcpp_SEMCpp. Models of this class
#' can be generated with the SEMFromLavaan-function.
#' @param subsets list with subsets created with createSubsets()
#' @param raw controls if the internal transformations of aCV4SEM should be used.
#' @param regularizedParameterLabels vector with labels of regularized parameters
#' @param lambda value of tuning parameter lambda
#' @param alpha value of tuning parameter alpha (for elastic net)
#' @param adaptiveLassoWeights vector with labeled adaptive lasso weights. Only required if penalty = "adaptiveLasso"
#' @param maxIter maximal number of iterations used in GLMENT
#' @param epsBreak breaking condition of GLMNET
GLMNETACVRcpp_SEMCpp <- function(SEM, 
                                 subsets, 
                                 raw = FALSE, 
                                 regularizedParameterLabels,
                                 lambda,
                                 alpha = NULL,
                                 adaptiveLassoWeights = NULL,
                                 hessianOfDifferentiablePart = NULL,
                                 maxIter = 100,
                                 epsBreak = 1e-10){
  
  if(!is(SEM, "Rcpp_SEMCpp")){
    stop("SEM must be of class Rcpp_SEMCpp")
  }
  
  parameters <- aCV4SEM:::getParameters(SEM = SEM, raw = raw)
  dataSet <- SEM$rawData
  N <- nrow(dataSet)
  k <- ncol(subsets)
  
  # compute derivatives of -2log-Likelihood without penalty
  
  scores <- aCV4SEM:::getScores(SEM = SEM, raw = raw)
  if(is.null(hessianOfDifferentiablePart)){
    hessian <- aCV4SEM:::getHessian(SEM = SEM, raw = raw)
  }else{
    hessian <- hessianOfDifferentiablePart
  }
  
  # Now do the GLMNET inner iteration for each sub-group
  stepdirections <- matrix(NA, nrow = k, ncol = length(parameters))
  colnames(stepdirections) <- names(parameters)
  
  #subsetParameters <- vector("list",k)
  subsetParameters <- matrix(NA, nrow = k, ncol = length(parameters))
  colnames(subsetParameters) <- names(parameters)
  rownames(subsetParameters) <- paste0("sample", 1:k)
  leaveOutFits <- rep(0, k)
  names(leaveOutFits) <- paste0("sample", 1:k)
  
  for(s in 1:k){

    subGroupGradient <- apply(scores[-c(which(subsets[,s])),],2,sum) # gradients of all inidividuals but the ones in the subgroup
    subGroupHessian <- ((N-sum(subsets[,s]))/N)*hessian # approximated hessian for all but the subgroup
    subGroupLambda <- alpha*(N-sum(subsets[,s]))*lambda
    
    # if elastic net is used, we approximate the ridge part as well and add this here:
    if(alpha != 1){
      # ridge or elastic net are used
      if(length(unique(adaptiveLassoWeights))!= 1) warning("Combining ridge and elastic net with adaptive lasso weights is unexplored territory.")
      # we multiply with the adaptive lasso weights.
      penaltyFunctionArguments <- list("regularizedParameterLabels" = regularizedParameterLabels,
                                       "lambda" = unique(adaptiveLassoWeights)*lambda*(1-alpha)*(N-length(subsets[[s]])))
      subGroupGradient <- subGroupGradient + aCV4SEM::ridgeGradient(parameters = parameters,
                                                                       penaltyFunctionArguments = penaltyFunctionArguments)
      if(is.null(hessianOfDifferentiablePart)){
        subGroupHessian <- subGroupHessian + aCV4SEM::ridgeHessian(parameters = parameters,
                                                                      penaltyFunctionArguments = penaltyFunctionArguments)
      }
    }
    
    startDirection <- Sys.time()
    direction <- aCV4SEM:::innerGLMNET(parameters = parameters, 
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
    
    # compute out of sample fit
    parameters_s <- parameters + stepdirections[s,names(parameters)]
    SEM <- aCV4SEM:::setParameters(SEM = SEM, labels = names(parameters), values = parameters_s, raw = raw)
    #SEM <- try(aCV4SEM:::fit(SEM), silent = TRUE)
    subsetParameters[s,names(parameters)] <- aCV4SEM:::getParameters(SEM = SEM, raw = FALSE)[names(parameters)]
   
    
    for(i in which(subsets[,s])){
      leaveOutFits[s] <- leaveOutFits[s] + aCV4SEM:::individualMinus2LogLikelihood(par = subsetParameters[s,], 
                                                                                   SEM = SEM, 
                                                                                   data = dataSet[i,], 
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


#' smoothACVRcpp_SEMCpp
#' 
#' internal function for approximate cross-validation based on the internal model representation of aCV4SEM. The smooth part refers to the fact
#' that this function expectes the penalty functions to be smooth. If the true penalty function is non-differentiable (e.g., lasso) a smooth
#' approximation of this function should be provided. See aCV4SEM::smoothLASSO as an example. Also, the gradients and the Hessian of this 
#' smooth approximation should be provided.
#' 
#' @param SEM model of class Rcpp_SEMCpp. Models of this class
#' can be generated with the SEMFromLavaan-function.
#' @param k the number of cross-validation folds. We recommend leave-one-out cross-validation; i.e. set k to the number of persons in the data set
#' @param individualPenaltyFunction penalty function which takes the current parameter values as first argument and the penaltyFunctionArguments as second argument and 
#' returns a single value - the value of the penalty function for a single person. If the true penalty function is non-differentiable (e.g., lasso) a smooth
#' approximation of this function should be provided.
#' @param individualPenaltyFunctionGradient gradients of the penalty function. Function should take the current parameter values as first argument and the penaltyFunctionArguments as second argument and 
#' return a vector of the same length as parameters. If the true penalty function is non-differentiable (e.g., lasso) a smooth
#' approximation of this function should be provided.
#' @param individualPenaltyFunctionHessian Hessian of the penalty function. Function should take the current parameter values as first argument and the penaltyFunctionArguments as second argument and 
#' return a matrix with (length as parameters)^2 number of elements. If the true penalty function is non-differentiable (e.g., lasso) a smooth
#' approximation of this function should be provided.
#' @param raw controls if the internal transformations of aCV4SEM should be used.
#' @param penaltyFunctionArguments can be anything that the functions individualPenaltyFunction, individualPenaltyFunctionGradient, or individualPenaltyFunctionHessian need. See aCV4SEM::smoothLASSO for an example.
smoothACVRcpp_SEMCpp <- function(SEM, 
                                 k, 
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
  
  parameters <- aCV4SEM:::getParameters(SEM = SEM, raw = raw)
  dataSet <- SEM$rawData
  N <- nrow(dataSet)
  
  # step 1: compute scores
  scores <- aCV4SEM:::getScores(SEM = SEM, raw = raw)
  
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
  
  hessian <- aCV4SEM:::getHessian(SEM = SEM, raw = raw)
  
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
      SEM <- aCV4SEM:::setParameters(SEM = SEM, labels = names(parameters), values = parameters_s, raw = raw)
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
      leaveOutFits[s] <- leaveOutFits[s] + aCV4SEM:::individualMinus2LogLikelihood(par = subsetParameters[[s]], 
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

#' createSubsets
#' 
#' create subsets for cross-validation
#' @param N number of samples in the data set
#' @param k number of subsets to create
#' @return matrix with subsets
createSubsets <- function(N,k){
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
  
  subsetMatrix <- matrix(NA, nrow = N, ncol = k,
                         dimnames = list(paste0("N",1:N), 
                                         paste0("subset", 1:k)))
  for(s in 1:length(subsets)){
    subsetMatrix[,s] <- 1:N %in% subsets[[s]]
  }
  
  if(any(apply(subsetMatrix,1,sum) != 1)) stop("Error while splitting data in subsets: Some persons are in multiple or none of the subsets")
  
  return(subsetMatrix)
}