approximateCrossValidation.lavaan <- function(model, k = NULL, par = NULL){
  
  if(model@Data@data.type != "full") stop("requires model@Data@data.type == 'full'")
  if(model@Data@ngroups != 1) stop("model@Data@ngroups != 1 currently not supported")
  rawData <- model@Data@X[[1]]
  colnames(rawData) <- model@Data@ov.names[[1]]
  N <- nrow(rawData)
  
  parameters <- parameterEstimates(model, remove.nonfree = TRUE)
  parameterLabels <- paste0(parameters$lhs, parameters$op, parameters$rhs)
  parameterLabels[parameters$label != ""] <- parameters$label[parameters$label != ""]
  parameterValues <- parameters$est
  names(parameterValues) <- parameterLabels
  
  # the folling approach is taken from the ipcr package; see ipcr:::bread_ipcr.lavaan
  if(model@Model@eq.constraints) {
    constraintMatrix <- lavaan:::lav_constraints_R2K(model@Model)
  }else{
    constraintMatrix <- diag(length(parameterValues))
  }
  
  # remove duplicated values due to equality constraints
  parameterValuesUnique <- parameterValues[unique(parameterLabels)]
  
  # we want to change the parameter values to those in par, if they are provided.
  if(!is.null(par)){
    if(!all(names(par) %in% names(parameterValuesUnique))) stop("Not all labels of par were found in the model.")
    parameterValuesUnique <- par[names(parameterValuesUnique)]
  }
  
  # extract some model contents for convenience
  lavmodel = model@Model
  lavsamplestats = model@SampleStats
  lavdata = model@Data
  lavoptions <- model@Options
  
  # change parameters; lavaan expects a vector with duplicates for parameters with constraints
  parametersExtended <- c(constraintMatrix%*%matrix(parameterValuesUnique, ncol = 1))
  names(parametersExtended) <- parameterLabels
  lavmodel <- lav_model_set_parameters(lavmodel, parametersExtended)
  
  # compute hessian
  # Note: the hessian returned by lavaan is for -N^(-1) log likelihood,
  # we will multiply it with N
  hessian <- N*lavaan:::lav_model_hessian(lavmodel = lavmodel, 
                                          lavsamplestats = lavsamplestats, 
                                          lavdata = lavdata, 
                                          lavoptions = lavoptions)
  
  # reduce to actual parameters by summing duplicate elements
  hessian <- t(constraintMatrix)%*%hessian%*%(constraintMatrix)
  
  
  # also: update the model object
  model@Model <- lavmodel
  model@implied <- lavaan:::lav_model_implied(lavmodel = lavmodel)
  
  # now compute the scores
  scores <- lavaan:::lavScores(object = model)
  # just as a check: the column-wise sum of the scores should be identical to
  # -N times the gradients (N times because the scores are for the log Likelihood
  # and the gradients are for -N^(-1) log likelihood)
  gradients <- lavaan:::lav_model_gradient(lavmodel = lavmodel, 
                                           lavsamplestats = lavsamplestats, 
                                           lavdata = lavdata)
  if(any(abs(-N*gradients%*%constraintMatrix - apply(scores,2,sum)) > .00001)) warning("There might be something wrong with the gradients here...")
  
  # approximate the parameter estimates when a subgroup is removed:
  if(is.null(k)) k <- N
  
  if(k < N){
    
    randomCases <- sample(1:N,N)
    subsets <- split(randomCases, sort(randomCases%%k)) # https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks
    
  }else if(k == N){
    
    subsets <- vector("list",N)
    names(subsets) <- 1:N
    for(i in 1:N) subsets[[i]] <- i
    
  }
  
  subsetParameters <- vector("list",k)
  leaveOutFits <- rep(0, k)
  names(leaveOutFits) <- paste0("sample", 1:k)
  
  HessianInverse <- solve((N/(N-1))*hessian) # approximation for hessian when removing one person
  
  for(s in seq_len(k)){
    
    direction <- -apply(scores[-subsets[[s]],],2,sum)%*%HessianInverse 
    
    # Taking a step size of 1 can result in nonsensical parameters. Therefore, a rudimentary line search is used here:
    stepLength <- 1
    
    parameters_s <- parameterValuesUnique + stepLength*direction
    names(parameters_s) <- names(parameterValuesUnique)
    subsetParameters[[s]] <- parameters_s
    
    # set parameters to sub-group parameters
    parametersExtended_s <- c(constraintMatrix%*%matrix(parameters_s, ncol = 1))
    names(parametersExtended_s) <- parameterLabels
    lavmodel_s <- lav_model_set_parameters(lavmodel, parametersExtended_s)
    implied_s <- lavaan:::lav_model_implied(lavmodel = lavmodel_s)
    
    for(i in 1:length(subsets[[s]])){
      individualData <- rawData[subsets[[s]][i],]
      isMissing <- is.na(individualData)
      
      leaveOutFits[s] <- leaveOutFits[s] + computeIndividualM2LL(nObservedVariables = sum(!isMissing),  
                                                                 rawData = individualData[!isMissing], 
                                                                 expectedMeans = implied_s$mean[[1]][!isMissing], 
                                                                 expectedCovariance = implied_s$cov[[1]][!isMissing,!isMissing])
    }
    
  }
  
  return(list(
    "leaveOutFits" = leaveOutFits,
    "subsetParameters" = subsetParameters,
    "scores" = scores
  )
  )
}


approximateCrossValidationWithPenalty.lavaan <- function(model, k = NULL, par = NULL, penaltyScoreFunction = NULL, penaltyFunction = NULL, penaltyFunctionDerivative = NULL, ...){

  if(sum(is.null(penaltyScoreFunction) + is.null(penaltyFunction)) == 1) stop("Both, penaltyScoreFunction and penaltyFunction must be provided.")
  
  parameters <- getParameters(SEM = SEM)
  dataSet <- SEM$data$rawData
  N <- nrow(dataSet)
  
  # step 1: compute scores
  scores <- computeAnalyticScores(par = parameters, SEM = SEM)
  # compute penalty scores
  if(!is.null(penaltyScoreFunction)) scores <- scores + penaltyScoreFunction(parameters, ...)[,colnames(scores)]
  
  # compute Hessian
  if(!is.null(penaltyScoreFunction)){
    hessian <- computeRegularizedHessian(par = parameters, 
                                         SEM = SEM, 
                                         penaltyFunction = penaltyFunction, 
                                         penaltyFunctionDerivative = penaltyFunctionDerivative, 
                                         ... = ...)
  }else{
    hessian <- computeHessianFromAnalytic(par = parameters, SEM = SEM)
  }
  
  
  # step 2: subgroup-parameters
  if(k < N){
    
    randomCases <- sample(1:N,N)
    subsets <- split(randomCases, sort(randomCases%%k)) # https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks
    
  }else if(k == N){
    
    subsets <- vector("list",N)
    names(subsets) <- 1:N
    for(i in 1:N) subsets[[i]] <- i
    
  }
  
  subsetParameters <- vector("list",k)
  leaveOutFits <- rep(0, k)
  names(leaveOutFits) <- paste0("sample", 1:k)
  
  inverseHessian <- (N/(N-1))*solve(hessian)
  for(s in 1:k){
    
    direction <- -apply(scores[-subsets[[s]],],2,sum)%*%inverseHessian 
    
    # Taking a step size of 1 can result in nonsensical parameters. Therefore, a rudimentary line search is used here:
    stepLength <- 1
    steps <- 1
    while(steps<100){
      parameters_s <- parameters + stepLength*direction
      SEM_s <- setParameters(SEM = SEM, labels = names(parameters), values = parameters_s)
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
      leaveOutFits[s] <- leaveOutFits[s] + individualMinus2LogLikelihood(par = subsetParameters[[s]], SEM = SEM, data = dataSet[subsets[[s]][i],])
    }
    
  }
  
  return(list("leaveOutFits" = leaveOutFits,
              "subsets" = subsets,
              "scores" = scores,
              "subsetParameters" = subsetParameters
  ))
  
}