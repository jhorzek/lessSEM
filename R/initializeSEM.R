SEMFromLavaan <- function(model, rawData, transformVariances = TRUE){
  
  # extract basic features
  meanstructure <- model@Options$meanstructure
  
  latentNames <- lavNames(model, type = "lv")
  manifestNames <- lavNames(model, type = "ov")
  
  nLatent <- length(latentNames)
  nManifest <- length(manifestNames)
  
  # data
  ## subset data
  rawData <- subset(rawData, select = manifestNames)
  if(!is.matrix(rawData)) rawData <- as.matrix(rawData)
  if(!is.numeric(rawData)) stop("rawData must be numeric")
  
  ## missing data
  FIML <- model@Options$missing == "ml"
  if(!FIML && anyNA(rawData)){
    
    warning("Deleting cases with missing values. Refit the lavaan object with missing = 'ML' to use full information maximum likelihood")
    rawData <- rawData[!apply(rawData,1,anyNA),]
    
  }
  
  internalData <- SEMdata(rawData)
  
  # Extract Model
  
  lavaanParameterTable <- lavInspect(model, what = "parTable")
  
  # translate to RAM notation
  
  ## directed paths
  AmatrixElements <- setAMatrix(model = model, 
                                lavaanParameterTable = lavaanParameterTable, 
                                nLatent = nLatent, 
                                nManifest = nManifest, 
                                latentNames = latentNames, 
                                manifestNames = manifestNames)
  
  ## undirected paths
  SmatrixElements <- setSMatrix(model = model, 
                                lavaanParameterTable = lavaanParameterTable, 
                                nLatent = nLatent, 
                                nManifest = nManifest, 
                                latentNames = latentNames, 
                                manifestNames = manifestNames)
  
  
  ## Mean structure
  MvectorElements <- setMVector(model = model, 
                                lavaanParameterTable = lavaanParameterTable, 
                                nLatent = nLatent, 
                                nManifest = nManifest, 
                                latentNames = latentNames, 
                                manifestNames = manifestNames,
                                rawData = rawData)
  
  ## Filter matrix
  Fmatrix <- matrix(0, 
                    nrow = nManifest, 
                    ncol = nLatent + nManifest, 
                    dimnames = list(manifestNames, 
                                    c(latentNames, manifestNames))
  )
  diag(Fmatrix[manifestNames, manifestNames]) <- 1
  
  ## Combine in list:
  modelMatrices <- list("Amatrix" = AmatrixElements$Amatrix,
                        "Smatrix" = SmatrixElements$Smatrix,
                        "Mvector" = MvectorElements$Mvector,
                        "Fmatrix" = Fmatrix)
  if(meanstructure){
    modelParameters <- list("Amatrix" = AmatrixElements$AmatrixParameters,
                            "Smatrix" = SmatrixElements$SmatrixParameters,
                            "Mvector" = MvectorElements$MvectorParameters)
  }else{
    modelParameters <- list("Amatrix" = AmatrixElements$AmatrixParameters,
                            "Smatrix" = SmatrixElements$SmatrixParameters)
  }
  
  # extract parameters
  parameterIDs <- model@ParTable$id[model@ParTable$free != 0]
  ops <- model@ParTable$op
  parameterLabels <- paste0(model@ParTable$lhs, ops, model@ParTable$rhs)
  parameterLabels[model@ParTable$label!=""] <- model@ParTable$label[model@ParTable$label!=""]
  parameterLabels <- parameterLabels[model@ParTable$free != 0]
  parameterValues <- model@ParTable$est[model@ParTable$free != 0]
  
  opsAlternative <- model@ParTable$op
  opsAlternative[opsAlternative == "=~"] <- "->"
  alternativeLabels <- paste(model@ParTable$lhs, opsAlternative, model@ParTable$rhs)
  alternativeLabels[opsAlternative == "~1"] <- paste("1 ->", model@ParTable$lhs[opsAlternative == "~1"])
  alternativeLabels[model@ParTable$label!=""] <- model@ParTable$label[model@ParTable$label!=""]
  alternativeLabels <- alternativeLabels[model@ParTable$free != 0]
  
  # construct internal representation of parameters
  nParameters <- length(parameterLabels)
  parameterTable <- data.frame("label" = c(),
                               "location" = c(), 
                               "row" = c(), 
                               "col" = c(),
                               "value" = c())
  
  matrixNames <- names(modelParameters)
  
  for(parameter in 1:nParameters){
    
    for(matrixName in matrixNames){
      
      parameterAt <- which(modelParameters[[matrixName]] == parameterIDs[parameter])
      
      if(parameterIDs[parameter] %in% modelParameters[[matrixName]]){
        
        for(ro in 1:nrow(modelParameters[[matrixName]])){
          for(co in 1:ncol(modelParameters[[matrixName]])){
            if(modelParameters[[matrixName]][ro,co] != parameterIDs[parameter]) next
            rawParameterValue <- parameterValues[parameter]
            # keep variances positive:
            if(transformVariances && ro==co && matrixName=="Smatrix") rawParameterValue <- ifelse(rawParameterValue < 0, log(.01), log(rawParameterValue))
            parameterTable <- rbind(parameterTable,
                                    data.frame(
                                      "label" = parameterLabels[parameter],
                                      "alternativeLabel" = alternativeLabels[parameter],
                                      "location" = matrixName, 
                                      "row" = ro, 
                                      "col" = co,
                                      "value" = parameterValues[parameter],
                                      "rawValue" = rawParameterValue
                                    )
            )
          }
        }
      }
    }
  }
  
  # if no mean structure: add 
  if(!meanstructure) {
    for(manifestName in manifestNames){
      parameterTable <- rbind(parameterTable,
                              data.frame(
                                "label" = paste0("nu_", manifestName),
                                "alternativeLabel" = paste0("nu_", manifestName),
                                "location" = "Mvector", 
                                "row" = which(rownames(MvectorElements$Mvector) == manifestName), 
                                "col" = 1,
                                "value" = MvectorElements$Mvector[manifestName,1],
                                "rawValue" = MvectorElements$Mvector[manifestName,1]
                              )
      )
    }
  }
  
  # saturated model fit
  N <- nrow(rawData)
  saturatedFit <- rep(NA,N)
  obsMeans <- NULL
  obsCov <- NULL
  
  if(!anyNA(rawData)) {
    
    nParameters <- ncol(rawData)
    
    obsMeans <- apply(rawData,2,mean)
    obsCov <- ((N - 1)/N)*cov(rawData)
    
    for(i in 1:N){
      saturatedFit[i] <- computeIndividualM2LL(nObservedVariables = nParameters, 
                                               rawData = rawData[i,],
                                               impliedMeans = obsMeans,
                                               impliedCovariance = obsCov)
    }
  }
  
  # build model
  
  SEMCpp <- new(SEMCpp)
  
  # set matrices and vector
  SEMCpp$setMatrix("A", modelMatrices$Amatrix)
  SEMCpp$setMatrix("S", modelMatrices$Smatrix)
  SEMCpp$setMatrix("F", modelMatrices$Fmatrix)
  SEMCpp$setVector('m', modelMatrices$Mvector)
  
  # set parameters
  SEMCpp$initializeParameters(parameterTable$label,
                              parameterTable$location,
                              parameterTable$row-1, # c++ starts at 0 
                              parameterTable$col-1, # c++ starts at 0 
                              parameterTable$value,
                              parameterTable$rawValue)
  
  # set derivative elements
  for(p in unique(parameterTable$label)){
    
    select <- parameterTable$label == p
    uniqueLocation <- unique(parameterTable$location[select])
    if(length(uniqueLocation) != 1) stop("Any parameter can only be in either A, S, or m")
    rows <- parameterTable$row[select]
    cols <- parameterTable$col[select]

    if(uniqueLocation == "Amatrix"){
      positionMatrix <- modelMatrices$Amatrix
      positionMatrix[] <- 0
      positionMatrix[cbind(rows,cols)] <- 1
      isVariance <- FALSE
    }else if(uniqueLocation == "Smatrix"){
      positionMatrix <- modelMatrices$Smatrix
      positionMatrix[] <- 0
      positionMatrix[cbind(rows,cols)] <- 1
      if(all(rows == cols)) isVariance <- TRUE
    }else if(uniqueLocation == "Mvector"){
      positionMatrix <- modelMatrices$Mvector
      positionMatrix[] <- 0
      positionMatrix[cbind(rows,cols)] <- 1
      isVariance <- FALSE
    }
    
    SEMCpp$addDerivativeElement(p, 
                                uniqueLocation, 
                                isVariance, 
                                positionMatrix)
  }
  
  # set data
  SEMCpp$addRawData(rawData, manifestNames)
  
  for(s in 1:length(internalData$missingSubsets)){
    SEMCpp$addSubset(internalData$missingSubsets[[s]]$N,
                     internalData$missingSubsets[[s]]$observed,
                     internalData$missingSubsets[[s]]$notMissing,
                     internalData$missingSubsets[[s]]$covariance,
                     internalData$missingSubsets[[s]]$means,
                     internalData$missingSubsets[[s]]$rawNoNA)
  }
  
  # check model
  if(round(SEMCpp$fit() - (-2*logLik(model)), 4) !=0) stop("Error translating lavaan to internal model representation: Different fit in SEMCpp and lavaan")
  
  return(SEMCpp)
}

setAMatrix <- function(model, lavaanParameterTable, nLatent, nManifest, latentNames, manifestNames){
  
  Amatrix <- matrix(0, 
                    nrow = nLatent + nManifest, 
                    ncol = nLatent + nManifest, 
                    dimnames = list(c(latentNames, manifestNames), 
                                    c(latentNames, manifestNames))
  )
  
  AmatrixParameters <- Amatrix
  
  if(!is.null(model@Model@GLIST$beta))  {
    Amatrix[latentNames, latentNames] <- model@Model@GLIST$beta # latent directed effects
    AmatrixParameters[latentNames, latentNames] <- lavaanParameterTable$beta
  }
  if(!is.null(model@Model@GLIST$lambda)) {
    AmatrixParameters[manifestNames, latentNames] <- lavaanParameterTable$lambda
    Amatrix[manifestNames, latentNames] <- model@Model@GLIST$lambda # loadings
  }
  
  return(list("Amatrix" = Amatrix, "AmatrixParameters" = AmatrixParameters))
  
}

setSMatrix <- function(model, lavaanParameterTable, nLatent, nManifest, latentNames, manifestNames){
  
  Smatrix <- matrix(0, 
                    nrow = nLatent + nManifest, 
                    ncol = nLatent + nManifest, 
                    dimnames = list(c(latentNames, manifestNames), 
                                    c(latentNames, manifestNames))
  )
  SmatrixParameters <- Smatrix
  
  if(!is.null(model@Model@GLIST$psi)){
    Smatrix[latentNames, latentNames] <- model@Model@GLIST$psi
    SmatrixParameters[latentNames, latentNames] <- lavaanParameterTable$psi
  }
  if(!is.null(model@Model@GLIST$theta)){
    Smatrix[manifestNames, manifestNames] <- model@Model@GLIST$theta
    SmatrixParameters[manifestNames, manifestNames] <- lavaanParameterTable$theta
  }
  
  return(list("Smatrix" = Smatrix, "SmatrixParameters" = SmatrixParameters))
  
}

setMVector <- function(model, lavaanParameterTable, nLatent, nManifest, latentNames, manifestNames, rawData){
  
  Mvector <- matrix(0, 
                    nrow = nLatent + nManifest, 
                    ncol = 1, 
                    dimnames = list(c(latentNames, manifestNames), 
                                    "1")
  )
  
  MvectorParameters <- Mvector
  
  if(!is.null(model@Model@GLIST$alpha)) {
    Mvector[latentNames, ] <- model@Model@GLIST$alpha # latent means
    MvectorParameters[latentNames, ] <- lavaanParameterTable$alpha
  }else{
    Mvector[latentNames, ] <- rep(0, length(latentNames))
  }
  
  if(!is.null(model@Model@GLIST$nu)) {
    Mvector[manifestNames, ] <- model@Model@GLIST$nu # manifest means
    MvectorParameters[manifestNames, ] <- lavaanParameterTable$nu
  }else{
    Mvector[manifestNames, ] <- apply(rawData, 2, mean, na.rm = TRUE)
  }
  
  return(list("Mvector" = Mvector, "MvectorParameters" = MvectorParameters))
}

