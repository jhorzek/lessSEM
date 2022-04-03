computeAnalyticScores <- function(par, SEM, raw = FALSE){
  
  SEM <- setParameters(SEM, labels = names(par), values = par, raw = raw)
  SEM <- computeExpected(SEM)
  
  parameters <- getParameters(SEM, raw = raw)
  parameterList <- SEM$model$parameters
  
  IminusAInverse <- solve(diag(nrow(SEM$model$matrices$Amatrix))-SEM$model$matrices$Amatrix)
  
  # for convenience:
  Amatrix <- SEM$model$matrices$Amatrix
  Smatrix <- SEM$model$matrices$Smatrix
  Fmatrix <- SEM$model$matrices$Fmatrix
  Mvector <- SEM$model$matrices$Mvector
  expectedCovariance <- SEM$model$expected$covariance
  expectedMeans <- SEM$model$expected$means
  
  # For each parameter, compute the derivative of the expected covariance
  derivativesOfCovariance <- computeExpectedCovarianceDerivatives(expectedCovariance = expectedCovariance, 
                                                                  parameters = parameters, 
                                                                  parameterList = parameterList, 
                                                                  Amatrix = Amatrix, 
                                                                  Smatrix = Smatrix, 
                                                                  Fmatrix = Fmatrix, 
                                                                  Mvector = Mvector, 
                                                                  IminusAInverse = IminusAInverse,
                                                                  raw = raw)
  
  # log-det-Sigma-derivative is independent of the person,
  # but depends on the missing structure! 
  logDetSigmas <- computeLogDetSigmas(uniqueMissingPatterns = SEM$data$internalData$uniqueMissingPatterns,
                                      parameters = parameters,
                                      parameterList = parameterList, 
                                      expectedCovariance = expectedCovariance,
                                      derivativesOfCovariance = derivativesOfCovariance)
  expectedCovInverses <- logDetSigmas$expectedCovInverses
  logDetSigmas <- logDetSigmas$logDetSigmas
  
  ## compute individual gradients
  
  N <- nrow(SEM$data$rawData)
  individualGradients <- matrix(0, 
                                nrow = N, 
                                ncol = length(parameters),
                                dimnames = list(paste0("person_", 1:N),
                                                names(parameters)))
  for(individual in 1:N){
    x <- SEM$data$rawData[individual,]
    isMissing <- SEM$data$internalData$uniqueMissingPatterns[SEM$data$internalData$individualMissingPatternID[individual],]
    individualGradients[individual, ] <- computeSingleSubjectGradient_Internal(x = x[!isMissing], 
                                                                               parameters = parameters, 
                                                                               parameterList = parameterList, 
                                                                               Amatrix = Amatrix, 
                                                                               Smatrix = Smatrix, 
                                                                               Fmatrix = Fmatrix, 
                                                                               Mvector = Mvector, 
                                                                               IminusAInverse = IminusAInverse,
                                                                               isMissing = isMissing,
                                                                               missingPattern = SEM$data$internalData$individualMissingPatternID[individual],
                                                                               expectedCovInverses = expectedCovInverses,
                                                                               expectedMeans = expectedMeans, 
                                                                               logDetSigmas = logDetSigmas, 
                                                                               derivativesOfCovariance = derivativesOfCovariance)
    
  }
  return(individualGradients)
}

computeAnalyticGradients <- function(par, SEM, raw = FALSE){
  return(apply(computeAnalyticGradientsByGroup(par = par, SEM = SEM, raw = raw),2,sum))
}

computeAnalyticGradientsByGroup <- function(par, SEM, raw){
  SEM <- setParameters(SEM, labels = names(par), values = par, raw = raw)
  SEM <- computeExpected(SEM)
  
  parameters <- getParameters(SEM)
  parameterList <- SEM$model$parameters
  
  IminusAInverse <- solve(diag(nrow(SEM$model$matrices$Amatrix))-SEM$model$matrices$Amatrix)
  
  # for convenience:
  Amatrix <- SEM$model$matrices$Amatrix
  Smatrix <- SEM$model$matrices$Smatrix
  Fmatrix <- SEM$model$matrices$Fmatrix
  Mvector <- SEM$model$matrices$Mvector
  expectedCovariance <- SEM$model$expected$covariance
  expectedMeans <- SEM$model$expected$means
  
  # For each parameter, compute the derivative of the expected covariance
  derivativesOfCovariance <- computeExpectedCovarianceDerivatives(expectedCovariance = expectedCovariance, 
                                                                  parameters = parameters, 
                                                                  parameterList = parameterList, 
                                                                  Amatrix = Amatrix, 
                                                                  Smatrix = Smatrix, 
                                                                  Fmatrix = Fmatrix, 
                                                                  Mvector = Mvector, 
                                                                  IminusAInverse = IminusAInverse, 
                                                                  raw = raw)
  
  # log-det-Sigma-derivative is independent of the person,
  # but depends on the missing structure! 
  logDetSigmas <- computeLogDetSigmas(uniqueMissingPatterns = SEM$data$internalData$uniqueMissingPatterns,
                                      parameters = parameters,
                                      parameterList = parameterList, 
                                      expectedCovariance = expectedCovariance,
                                      derivativesOfCovariance = derivativesOfCovariance)
  expectedCovInverses <- logDetSigmas$expectedCovInverses
  logDetSigmas <- logDetSigmas$logDetSigmas
  
  ## compute gradients for each subgroup
  NGroups <- nrow(SEM$data$internalData$uniqueMissingPatterns)
  groupGradients <- matrix(0, 
                           nrow = NGroups, 
                           ncol = length(parameters),
                           dimnames = list(paste0("group_", 1:NGroups),
                                           names(parameters)))
  for(gr in 1:NGroups){
    NInSubgroup <- SEM$data$internalData$missingSubsets[[gr]]$N
    
    if(NInSubgroup == 1){
      
      groupGradients[gr, ] <- computeSingleSubjectGradient_Internal(x = SEM$data$internalData$missingSubsets[[gr]]$raw, 
                                                                    parameters = parameters, 
                                                                    parameterList = parameterList, 
                                                                    Amatrix = Amatrix, 
                                                                    Smatrix = Smatrix, 
                                                                    Fmatrix = Fmatrix, 
                                                                    Mvector = Mvector, 
                                                                    IminusAInverse = IminusAInverse,
                                                                    isMissing = SEM$data$internalData$missingSubsets[[gr]]$missing,
                                                                    missingPattern = gr,
                                                                    expectedCovInverses = expectedCovInverses,
                                                                    expectedMeans = expectedMeans, 
                                                                    logDetSigmas = logDetSigmas, 
                                                                    derivativesOfCovariance = derivativesOfCovariance)
      next
    }
    
    if(NInSubgroup > 1){
      
      groupGradients[gr, ] <- computeSubgroupGradient_Internal(
        means = SEM$data$internalData$missingSubsets[[gr]]$means,
        covariance = SEM$data$internalData$missingSubsets[[gr]]$covariance,
        N = NInSubgroup,
        parameters = parameters, 
        parameterList = parameterList, 
        Amatrix = Amatrix, 
        Smatrix = Smatrix, 
        Fmatrix = Fmatrix, 
        Mvector = Mvector, 
        IminusAInverse = IminusAInverse,
        isMissing = SEM$data$internalData$missingSubsets[[gr]]$missing,
        missingPattern = gr,
        expectedCovInverses = expectedCovInverses,
        expectedMeans = expectedMeans, 
        logDetSigmas = logDetSigmas, 
        derivativesOfCovariance = derivativesOfCovariance)
      
      next
    }
    
  }
  
  
  
  return(groupGradients)
  
  
}

## Helper Functions:

computeExpectedCovarianceDerivatives <- function(expectedCovariance, 
                                                 parameters, 
                                                 parameterList, 
                                                 Amatrix, Smatrix, Fmatrix, Mvector, 
                                                 IminusAInverse = NULL,
                                                 raw){
  if(is.null(IminusAInverse)) IminusAInverse <- solve(diag(nrow(Amatrix)) - Amatrix)
  
  derivativesOfCovariance <- array(0, dim = c(nrow(expectedCovariance),
                                              ncol(expectedCovariance),
                                              length(parameters)))
  
  # pre-compute some repeatedly used elements
  FIminusAInverse <- Fmatrix%*%IminusAInverse
  tFIminusAInverse <- t(FIminusAInverse)
  Aelement1 <- Smatrix%*%
    t(IminusAInverse)%*%
    t(Fmatrix) 
  Aelement2 <-  Fmatrix%*%IminusAInverse%*%
    Smatrix
  
  for(p in 1:length(parameters)){
    # locate parameter
    relevantRows <- parameterList$label == names(parameters)[p]
    parIn <- unique(parameterList$location[relevantRows])
    
    if(parIn == "Mvector") next
    if(parIn == "Smatrix"){
      Sderiv <- Smatrix
      Sderiv[] <- 0
      for(ro in which(relevantRows)){
        Sderiv[parameterList$row[ro], parameterList$col[ro]] <- ifelse(parameterList$row[ro] == parameterList$col[ro] && raw, parameterList$value[ro], 1)
      }
      derivativesOfCovariance[,,p] <- FIminusAInverse%*%Sderiv%*%tFIminusAInverse
      next
    }
    if(parIn == "Amatrix"){
      Aderiv <- Amatrix
      Aderiv[] <- 0
      for(ro in which(relevantRows)){
        Aderiv[parameterList$row[ro], parameterList$col[ro]] <- 1
      }
      derivativesOfCovariance[,,p] <- Fmatrix%*%(
        IminusAInverse%*%
          Aderiv%*%
          IminusAInverse)%*%
        Aelement1 + 
        Aelement2%*%
        t(IminusAInverse%*%
            Aderiv%*%
            IminusAInverse
        )%*%
        t(Fmatrix)
      next
    }
  }
  return(derivativesOfCovariance)
}

computeLogDetSigmas <- function(uniqueMissingPatterns,parameters,parameterList, expectedCovariance, derivativesOfCovariance){
  logDetSigmas <- matrix(0, 
                         nrow = nrow(uniqueMissingPatterns),
                         ncol = length(parameters),
                         dimnames = list(paste0("missingPattern", 1:nrow(uniqueMissingPatterns)),
                                         names(parameters)))
  expectedCovInverses <- vector("list", nrow(logDetSigmas))
  for(missingPattern in 1:nrow(logDetSigmas)){
    isMissing <- uniqueMissingPatterns[missingPattern,]
    currentExpectedCovariance <- expectedCovariance[!isMissing, !isMissing]
    expectedCovInverses[[missingPattern]] <- solve(currentExpectedCovariance)
    
    for(p in 1:length(parameters)){
      # step 1: check location of parameter
      parIn <- unique(parameterList$location[parameterList$label == names(parameters)[p]])
      
      if(parIn == "Mvector") next
      if(parIn == "Smatrix" || parIn == "Amatrix"){
        logDetSigmas[missingPattern, p] <- sum(diag(expectedCovInverses[[missingPattern]]%*%derivativesOfCovariance[,,p][!isMissing, !isMissing]))
        next
      }
      stop("Encountered parameter in unknown location")
    }
  }
  return(list("logDetSigmas" = logDetSigmas,
              "expectedCovInverses" = expectedCovInverses)
  )
}


computeSingleSubjectGradient_Internal <- function(x, 
                                                  parameters, 
                                                  parameterList, 
                                                  Amatrix, 
                                                  Smatrix,
                                                  Fmatrix, 
                                                  Mvector, 
                                                  IminusAInverse = NULL,
                                                  isMissing,
                                                  missingPattern,
                                                  expectedCovInverses,
                                                  expectedMeans, 
                                                  logDetSigmas, 
                                                  derivativesOfCovariance){
  if(anyNA(x)) stop("x must be the data without NAs")
  if(length(x) != sum(!isMissing)){stop("x must be of length sum(!isMissing)")}
  if(is.null(IminusAInverse)) IminusAInverse <- solve(diag(nrow(Amatrix)) - Amatrix)
  
  individualGradient <- rep(NA, length(parameters))
  names(individualGradient) <- names(parameters)
  
  
  
  for(p in 1:length(parameters)){
    # step 1: check location of parameter
    parIn <- unique(parameterList$location[parameterList$label == names(parameters)[p]])
    
    if(parIn == "Mvector"){
      relevantRows <- parameterList$label == names(parameters)[p]
      Mderiv <- Mvector
      Mderiv[] <- 0
      for(ro in which(relevantRows)){
        Mderiv[parameterList$row[ro], parameterList$col[ro]] <- 1
      }
      individualGradient[names(parameters)[p]] <- 2*matrix(-(Fmatrix%*%IminusAInverse%*%Mderiv)[!isMissing], nrow = 1)%*%
        expectedCovInverses[[missingPattern]]%*%
        matrix(x - expectedMeans[!isMissing], ncol = 1)
      next
    }
    if(parIn == "Smatrix"){
      individualGradient[names(parameters)[p]] <- logDetSigmas[missingPattern, p] + matrix(x - expectedMeans[!isMissing], nrow = 1)%*%
        (-expectedCovInverses[[missingPattern]])%*%
        derivativesOfCovariance[,,p][!isMissing, !isMissing]%*%
        expectedCovInverses[[missingPattern]]%*%
        matrix(x - expectedMeans[!isMissing], ncol = 1)
      next
    }
    if(parIn == "Amatrix"){
      relevantRows <- parameterList$label == names(parameters)[p]
      Aderiv <- Amatrix
      Aderiv[] <- 0
      for(ro in which(relevantRows)){
        Aderiv[parameterList$row[ro], parameterList$col[ro]] <- 1
      }
      
      individualGradient[names(parameters)[p]] <- logDetSigmas[missingPattern, p] +
        2*t((-Fmatrix%*%IminusAInverse %*%Aderiv%*%IminusAInverse%*%Mvector)[!isMissing])%*%
        expectedCovInverses[[missingPattern]]%*%
        matrix(x - expectedMeans[!isMissing], ncol = 1) +
        matrix(x - expectedMeans[!isMissing], nrow = 1)%*%
        (-expectedCovInverses[[missingPattern]])%*%
        derivativesOfCovariance[,,p][!isMissing, !isMissing]%*%
        expectedCovInverses[[missingPattern]]%*%
        matrix(x - expectedMeans[!isMissing], ncol = 1)
      next
    }
    stop("Encountered parameter in unknown location")
  }
  return(individualGradient)
}

computeSubgroupGradient_Internal <- function(
  means,
  covariance,
  N,
  parameters, 
  parameterList, 
  Amatrix, 
  Smatrix, 
  Fmatrix, 
  Mvector, 
  IminusAInverse,
  isMissing,
  missingPattern,
  expectedCovInverses,
  expectedMeans, 
  logDetSigmas, 
  derivativesOfCovariance){
  
  if(is.null(IminusAInverse)) IminusAInverse <- solve(diag(nrow(Amatrix)) - Amatrix)
  
  groupGradient <- rep(NA, length(parameters))
  names(groupGradient) <- names(parameters)
  
  for(p in 1:length(parameters)){
    # step 1: check location of parameter
    parIn <- unique(parameterList$location[parameterList$label == names(parameters)[p]])
    
    if(parIn == "Mvector"){
      relevantRows <- parameterList$label == names(parameters)[p]
      Mderiv <- Mvector
      Mderiv[] <- 0
      for(ro in which(relevantRows)){
        Mderiv[parameterList$row[ro], parameterList$col[ro]] <- 1
      }
      groupGradient[names(parameters)[p]] <- N*2*matrix(-(Fmatrix%*%IminusAInverse%*%Mderiv)[!isMissing], nrow = 1)%*%
        expectedCovInverses[[missingPattern]]%*%
        matrix(means - expectedMeans[!isMissing], ncol = 1)
      
      next
    }
    if(parIn == "Smatrix"){
      B <- (-expectedCovInverses[[missingPattern]])%*%
        derivativesOfCovariance[,,p][!isMissing, !isMissing]%*%
        expectedCovInverses[[missingPattern]]
      
      groupGradient[names(parameters)[p]] <- N*logDetSigmas[missingPattern, p] + 
        N*sum(diag(covariance%*%B)) +
        N*
        matrix(means - expectedMeans[!isMissing], nrow = 1)%*%
        B%*%
        matrix(means - expectedMeans[!isMissing], ncol = 1)
      
      next
    }
    if(parIn == "Amatrix"){
      relevantRows <- parameterList$label == names(parameters)[p]
      Aderiv <- Amatrix
      Aderiv[] <- 0
      for(ro in which(relevantRows)){
        Aderiv[parameterList$row[ro], parameterList$col[ro]] <- 1
      }
      
      B <- (-expectedCovInverses[[missingPattern]])%*%
        derivativesOfCovariance[,,p][!isMissing, !isMissing]%*%
        expectedCovInverses[[missingPattern]]
      
      groupGradient[names(parameters)[p]] <- 
        N*logDetSigmas[missingPattern, p] +
        N*2*t((-Fmatrix%*%IminusAInverse %*%Aderiv%*%IminusAInverse%*%Mvector)[!isMissing])%*%
        expectedCovInverses[[missingPattern]]%*%
        matrix(means - expectedMeans[!isMissing], ncol = 1) +
        
        N*sum(diag(covariance%*%B)) + 
        
        N*matrix(means - expectedMeans[!isMissing], nrow = 1)%*%
        B%*%
        matrix(means - expectedMeans[!isMissing], ncol = 1)
      
      next
    }
    stop("Encountered parameter in unknown location")
  }
  return(groupGradient)
  
}

#### Hessian ####

computeHessianFromAnalytic <- function (par, SEM, raw = FALSE, eps = 1e-7){
  # THE FOLLOWING CODE IS ADAPTED FROM LAVAAN. 
  # SEE lavaan:::lav_model_hessian FOR THE IMPLEMENTATION
  # BY Yves Rosseel
  
  SEM <- setParameters(SEM = SEM, labels = names(par), values = par, raw = raw)
  SEM <- fit(SEM = SEM)
  parameters <- getParameters(SEM = SEM, raw = raw) # in case the par passed by the user is only a subset of the parameters of the model
  nParameters <- length(parameters)
  Hessian <- matrix(data = 0, 
                    nrow = nParameters, 
                    ncol = nParameters, 
                    dimnames = list(names(parameters),
                                    names(parameters)))
  
  for (p in 1:nParameters) {
    
    stepLeft <- twoStepLeft <- stepRight <- twoStepRight <- parameters
    stepLeft[p] <- stepLeft[p] - eps
    twoStepLeft[p] <- twoStepLeft[p] - 2 * eps
    stepRight[p] <- stepRight[p] + eps
    twoStepRight[p] <- twoStepRight[p] + 2 * eps
    
    gradientsStepLeft <- computeAnalyticGradients(par = stepLeft, SEM = SEM, raw = raw)
    gradientsTwoStepLeft <- computeAnalyticGradients(par = twoStepLeft, SEM = SEM, raw = raw)
    gradientsStepRight <- computeAnalyticGradients(par = stepRight, SEM = SEM, raw = raw)
    gradientsTwoStepRight <- computeAnalyticGradients(par = twoStepRight, SEM = SEM, raw = raw)
    
    Hessian[,p] <- (gradientsTwoStepLeft - 
                      8 * gradientsStepLeft + 
                      8 * gradientsStepRight - 
                      gradientsTwoStepRight)/(12 * eps)
  }
  # make symmetric
  Hessian <- (Hessian + t(Hessian))/2
  
  return(Hessian)
  
}

# computeAnalyticHessian <- function(par, SEM){
#   stop("Requires further work")
#   SEM <- setParameters(SEM, labels = names(par), values = par)
#   SEM <- computeExpected(SEM)
#   
#   parameters <- getParameters(SEM)
#   parameterList <- SEM$model$parameters
#   
#   IminusAInverse <- solve(diag(nrow(SEM$model$matrices$Amatrix))-SEM$model$matrices$Amatrix)
#   
#   # for convenience:
#   Amatrix <- SEM$model$matrices$Amatrix
#   Smatrix <- SEM$model$matrices$Smatrix
#   Fmatrix <- SEM$model$matrices$Fmatrix
#   Mvector <- SEM$model$matrices$Mvector
#   expectedCovariance <- SEM$model$expected$covariance
#   expectedMeans <- SEM$model$expected$means
#   
#   # For each parameter, compute the derivative of the expected covariance
#   derivativesOfCovariance <- computeExpectedCovarianceDerivatives(expectedCovariance = expectedCovariance, 
#                                                                   parameters = parameters, 
#                                                                   parameterList = parameterList, 
#                                                                   Amatrix = Amatrix, 
#                                                                   Smatrix = Smatrix, 
#                                                                   Fmatrix = Fmatrix, 
#                                                                   Mvector = Mvector, 
#                                                                   IminusAInverse = IminusAInverse)
#   
#   # log-det-Sigma-derivative is independent of the person,
#   # but depends on the missing structure! 
#   logDetSigmas <- computeLogDetSigmas(uniqueMissingPatterns = SEM$data$internalData$uniqueMissingPatterns,
#                                       parameters = parameters,
#                                       parameterList = parameterList, 
#                                       expectedCovariance = expectedCovariance,
#                                       derivativesOfCovariance = derivativesOfCovariance)
#   expectedCovInverses <- logDetSigmas$expectedCovInverses
#   logDetSigmas <- logDetSigmas$logDetSigmas
#   
#   ## compute individual gradients
#   
#   N <- nrow(SEM$data$rawData)
#   individualGradients <- matrix(0, 
#                                 nrow = N, 
#                                 ncol = length(parameters),
#                                 dimnames = list(paste0("person_", 1:N),
#                                                 names(parameters)))
#   for(individual in 1:N){
#     x <- SEM$data$rawData[individual,]
#     individualGradients[individual, ] <- computeSingleSubjectGradient_Internal(x = x[!isMissing], 
#                                                                                parameters = parameters, 
#                                                                                parameterList = parameterList, 
#                                                                                Amatrix = Amatrix, 
#                                                                                Smatrix = Smatrix, 
#                                                                                Fmatrix = Fmatrix, 
#                                                                                Mvector = Mvector, 
#                                                                                IminusAInverse = IminusAInverse,
#                                                                                missingPattern = SEM$data$internalData$individualMissingPatternID[individual], 
#                                                                                uniqueMissingPatterns = SEM$data$internalData$uniqueMissingPatterns, 
#                                                                                expectedCovInverses = expectedCovInverses,
#                                                                                expectedMeans = expectedMeans, 
#                                                                                logDetSigmas = logDetSigmas, 
#                                                                                derivativesOfCovariance = derivativesOfCovariance)
#     
#   }
#   return(individualGradients)
#   
#   
# }
# 
# computeLogDetSigmasHessian <- function(uniqueMissingPatterns,
#                                        parameters,
#                                        parameterList, 
#                                        expectedCovariance, 
#                                        derivativesOfCovariance){
#   
#   logDetSigmasHessian <- array(0, dim = c(length(parameters), # first derivative
#                                           length(parameters), # second derivative
#                                           nrow(uniqueMissingPatterns))) # missingness patterns
#   
#   expectedCovInverses <- vector("list", nrow(logDetSigmasHessian))
#   for(missingPattern in 1:nrow(logDetSigmasHessian)){
#     isMissing <- uniqueMissingPatterns[missingPattern,]
#     currentExpectedCovariance <- expectedCovariance[!isMissing, !isMissing]
#     expectedCovInverses[[missingPattern]] <- solve(currentExpectedCovariance)
#     
#     for(p1 in 1:length(parameters)){
#       # step 1: check location of parameter
#       par1In <- unique(parameterList$location[parameterList$label == names(parameters)[p1]])
#       if(par1In == "Mvector") next
#       if(par1In == "Smatrix"){
#         relevantRows <- parameterList$label == names(parameters)[p1]
#         Sderiv <- Smatrix
#         Sderiv[] <- 0
#         for(ro in which(relevantRows)){
#           Sderiv[parameterList$row[ro], parameterList$col[ro]] <- 1
#         }
#       }
#       
#       for(p2 in 1:length(parameters)){
#         # step 1: check location of parameter
#         par2In <- unique(parameterList$location[parameterList$label == names(parameters)[p2]])
#         
#         if(par2In == "Mvector") next
#         if(par1In == "Smatrix" && par2In == "Smatrix"){
#           logDetSigmasHessian[p1, p2, missingPattern] <- sum(diag(
#             -expectedCovInverses[[missingPattern]]%*%
#               derivativesOfCovariance[,,p2][!isMissing, !isMissing]%*%
#               expectedCovInverses[[missingPattern]]%*%
#               Fmatrix%*%
#               IminusAInverse%*%
#               Sderiv%*%
#               t(IminusAInverse)%*%
#               t(Fmatrix)
#           ))
#           next
#         }
#         if(par1In == "Smatrix" && par2In == "Amatrix"){
#           
#           relevantRows <- parameterList$label == names(parameters)[p2]
#           Aderiv <- Amatrix
#           Aderiv[] <- 0
#           for(ro in which(relevantRows)){
#             Aderiv[parameterList$row[ro], parameterList$col[ro]] <- 1
#           }
#           
#           logDetSigmasHessian[p1, p2, missingPattern] <- sum(diag(
#             -expectedCovInverses[[missingPattern]]%*%
#               derivativesOfCovariance[,,p2][!isMissing, !isMissing]%*%
#               expectedCovInverses[[missingPattern]]%*%
#               Fmatrix%*%
#               IminusAInverse%*%
#               Sderiv%*%
#               t(IminusAInverse)%*%
#               t(Fmatrix) +
#               
#               expectedCovInverses[[missingPattern]]%*%
#               Fmatrix%*%
#               IminusAInverse%*%
#               Aderiv%*%
#               IminusAInverse%*%
#               Sderiv%*%
#               t(IminusAInverse)%*%
#               t(Fmatrix) +
#               
#               expectedCovInverses[[missingPattern]]%*%
#               Fmatrix%*%
#               IminusAInverse%*%
#               Sderiv%*%
#               t(IminusAInverse%*%
#                   Aderiv%*%
#                   IminusAInverse)%*%
#               Fmatrix
#           ))
#           next
#         }
#         
#       }
#       
#       #if(parIn == "Mvector") next
#       #if(parIn == "Smatrix" || parIn == "Amatrix"){
#       
#       next
#       #}
#       #stop("Encountered parameter in unknown location")
#     }
#   }
#   return(list("logDetSigmasHessian" = logDetSigmasHessian,
#               "expectedCovInverses" = expectedCovInverses)
#   )
# }