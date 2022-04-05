approximateCrossValidation <- function(SEM, k, individualPenaltyFunction = NULL, raw = FALSE, ...){
  
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

