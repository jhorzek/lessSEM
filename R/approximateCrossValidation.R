approximateCrossValidation <- function(SEM, k, penaltyScoreFunction = NULL, penaltyFunction = NULL, penaltyFunctionDerivative = NULL, raw = FALSE, ...){

  if(sum(is.null(penaltyScoreFunction) + is.null(penaltyFunction)) == 1) stop("Both, penaltyScoreFunction and penaltyFunction must be provided.")
  additionalArguments <- list(...)
  
  parameters <- getParameters(SEM = SEM, raw = raw)
  dataSet <- SEM$data$rawData
  N <- nrow(dataSet)
  
  # step 1: compute scores
  scores <- computeAnalyticScores(par = parameters, SEM = SEM, raw = raw)
  # compute penalty scores
  if(!is.null(penaltyScoreFunction)) scores <- scores + do.call(what = penaltyScoreFunction, 
                                                                args = c(list("par" = parameters),
                                                                  additionalArguments
                                                                  ))[,colnames(scores)]
  
  # compute Hessian
  if(!is.null(penaltyScoreFunction)){
    
    hessian <- do.call(what = computeRegularizedHessian, 
                       args = c(list("par" = parameters, 
                                     "SEM" = SEM, 
                                     "raw" = raw,
                                     "penaltyFunction" = penaltyFunction, 
                                     "penaltyFunctionDerivative" = penaltyFunctionDerivative),
                                additionalArguments
                       ))
  }else{
    hessian <- computeHessianFromAnalytic(par = parameters, raw = raw, SEM = SEM)
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
      SEM_s <- setParameters(SEM = SEM, labels = names(parameters), values = parameters_s, raw = raw)
      SEM_s <- try(computeExpected(SEM_s), silent = TRUE)
      if(any(class(SEM_s) == "try-error")) {
        stepLength <- stepLength*.9
        steps <- steps + 1
        next
      }
      # check positive definiteness
      tryChol <- try(chol(SEM_s$model$expected$covariance), silent = TRUE)
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

