#' .SEMFromLavaan
#' 
#' internal function. Translates an object of class lavaan to the internal model representation.
#' 
#' @param lavaanModel model of class lavaan
#' @param whichPars which parameters should be used to initialize the model. If set to "est", the parameters will be set to the
#' estimated parameters of the lavaan model. If set to "start", the starting values of lavaan will be used. The latter can be useful if parameters are to
#' be optimized afterwards as setting the parameters to "est" may result in the model getting stuck in a local minimum.
#' @param fit should the model be fitted and compared to the lavaanModel?
#' @param activeSet Option to only use a subset of the individuals in the data set. Logical vector of length N indicating which subjects should remain in the sample.
#' @param addMeans If lavaanModel has meanstructure = FALSE, addMeans = TRUE will add a mean structure. FALSE will set the means of the observed variables to the average
#' @param dataSet optional: Pass an alternative data set to lessSEM:::.SEMFromLavaan which will replace the original data set in lavaanModel.
#' @returns Object of class Rcpp_SEMCpp
#' @keywords internal
.SEMFromLavaan <- function(lavaanModel, 
                           whichPars = "est",
                           fit = TRUE,
                           addMeans = TRUE,
                           activeSet = NULL,
                           dataSet = NULL,
                           transformations = NULL,
                           transformationList = list()){
  if(!is(lavaanModel, "lavaan")) stop("lavaanModel must be of class lavaan.")
  
  SEMList <- list(
    matrices = list("A" = NULL,
                    "S" = NULL,
                    "F" = NULL,
                    "M" = NULL),
    parameters = list("label" = NULL,
                      "location" = NULL,
                      "row" = NULL,
                      "col" = NULL,
                      "value" = NULL,
                      "rawValue" = NULL,
                      "isTransformation" = NULL),
    DerivativeElements = list(),
    rawData = list("rawData" = NULL,
                   "personInSubset" = NULL,
                   "manifestNames" = NULL),
    subsets = list()
  )
  
  if(is.null(dataSet)){
    rawData <- try(lavaan::lavInspect(lavaanModel, "data"))
    if(is(rawData, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
    checkFit <- TRUE
  }else{
    lavaanData <- try(lavaan::lavInspect(lavaanModel, "data"))
    # make sure that the sorting of rawData is correct:
    rawData <- dataSet[,colnames(lavaanData),drop = FALSE]
    if(any(!colnames(lavaanData) %in% colnames(rawData))) 
      stop("Not all variables are in rawData.")
    checkFit <- FALSE
  }
  
  if(!is.null(activeSet)){
    if(length(activeSet) != nrow(rawData)) stop("length of activeSet must be identical to the rows in the data set.")
    if(!is.logical(activeSet)) stop("activeSet must be logical.")
    rawData <- rawData[activeSet,,drop = FALSE]
    checkFit <- FALSE
  }
  
  # extract basic features
  meanstructure <- lavaanModel@Options$meanstructure
  fixedX <- lavaan::lavInspect(lavaanModel, "fixed.x")
  if(fixedX) {
    warning("lavaanModel has option fixed.x set to TRUE. This is currently not fully supported by lessSEM. Be sure to check the results.")
    checkFit <- FALSE
  }
  
  lavaanParameterTable <- lavaan::lavInspect(lavaanModel, what = "parTable")
  
  latentNames <- colnames(lavaanParameterTable$lambda)  #lavNames(lavaanModel, type = "lv")
  manifestNames <- rownames(lavaanParameterTable$lambda) # lavNames(lavaanModel, type = "ov")
  
  # if regressions manifest -> latent are specified, lavaan will duplicate some latent and
  # manifest names; this may be confusing in RAM notation, so we take care of this here (see testing-mediation.R):
  latentNames[latentNames%in%manifestNames] <- paste0(latentNames[latentNames%in%manifestNames], "_lv")
  nLatent <- length(latentNames)
  nManifest <- length(manifestNames)
  
  # data
  ## subset data
  rawData <- subset(rawData, select = manifestNames)
  if(!is.matrix(rawData)) rawData <- as.matrix(rawData)
  if(!is.numeric(rawData)) stop("rawData must be numeric")
  
  ## missing data
  FIML <- lavaanModel@Options$missing == "ml"
  if(!FIML && anyNA(rawData)){
    
    warning("Deleting cases with missing values. Refit the lavaan object with missing = 'ML' to use full information maximum likelihood")
    rawData <- rawData[!apply(rawData,1,anyNA),]
    
  }
  
  internalData <- .SEMdata(rawData)
  
  # translate to RAM notation
  
  ## directed paths
  AmatrixElements <- .setAMatrix(model = lavaanModel, 
                                 lavaanParameterTable = lavaanParameterTable, 
                                 nLatent = nLatent, 
                                 nManifest = nManifest, 
                                 latentNames = latentNames, 
                                 manifestNames = manifestNames)
  
  ## undirected paths
  SmatrixElements <- .setSMatrix(model = lavaanModel, 
                                 lavaanParameterTable = lavaanParameterTable, 
                                 nLatent = nLatent, 
                                 nManifest = nManifest, 
                                 latentNames = latentNames, 
                                 manifestNames = manifestNames)
  
  
  ## Mean structure
  MvectorElements <- .setMVector(model = lavaanModel, 
                                 lavaanParameterTable = lavaanParameterTable, 
                                 nLatent = nLatent, 
                                 nManifest = nManifest, 
                                 latentNames = latentNames, 
                                 manifestNames = manifestNames,
                                 rawData = rawData)
  
  ## Filter matrix
  Fmatrix <- diag(nManifest)
  rownames(Fmatrix) <- manifestNames
  colnames(Fmatrix) <- manifestNames
  Fmatrix <- cbind(
    matrix(0, nrow = nManifest, ncol = nLatent, dimnames = list(manifestNames, 
                                                                latentNames)),
    Fmatrix)
  
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
  parameterIDs <- lavaanModel@ParTable$id[lavaanModel@ParTable$free != 0]
  ops <- lavaanModel@ParTable$op
  parameterLabels <- paste0(lavaanModel@ParTable$lhs, ops, lavaanModel@ParTable$rhs)
  parameterLabels[lavaanModel@ParTable$label!=""] <- lavaanModel@ParTable$label[lavaanModel@ParTable$label!=""]
  parameterLabels <- parameterLabels[lavaanModel@ParTable$free != 0]
  if(whichPars == "est"){
    parameterValues <- lavaanModel@ParTable$est[lavaanModel@ParTable$free != 0]
  }else if(whichPars == "start"){
    parameterValues <- lavaanModel@ParTable$start[lavaanModel@ParTable$free != 0]
  }else{
    stop(paste0("Could not set the parameters of the model. Set whichPars to one of: 'est', 'start'. See ?.SEMFromLavaan for more details."))
  }
  
  # construct internal representation of parameters
  nParameters <- length(parameterLabels)
  parameterTable <- data.frame("label" = c(),
                               "location" = c(), 
                               "row" = c(), 
                               "col" = c(),
                               "value" = c(),
                               "isTransformation" = c())
  
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
            if(ro==co && matrixName=="Smatrix") {
              if(rawParameterValue < 0){
                rawParameterValue <- log(.01)
                warning("lavaanModel has negative variances. Cannot compare fit to lavaanModel")
                fit <- FALSE
              }else{
                rawParameterValue <- log(rawParameterValue)
              }
              
            }
            parameterTable <- rbind(parameterTable,
                                    data.frame(
                                      "label" = parameterLabels[parameter],
                                      "location" = matrixName, 
                                      "row" = ro, 
                                      "col" = co,
                                      "value" = parameterValues[parameter],
                                      "rawValue" = rawParameterValue,
                                      "isTransformation" = FALSE
                                    )
            )
          }
        }
      }
    }
  }
  
  # if no mean structure: add 
  if(!meanstructure && addMeans) {
    for(manifestName in manifestNames){
      parameterTable <- rbind(parameterTable,
                              data.frame(
                                "label" = paste0(manifestName, "~1"),
                                "location" = "Mvector", 
                                "row" = which(rownames(MvectorElements$Mvector) == manifestName), 
                                "col" = 1,
                                "value" = MvectorElements$Mvector[manifestName,1],
                                "rawValue" = MvectorElements$Mvector[manifestName,1],
                                "isTransformation" = FALSE
                              )
      )
    }
  }
  
  # check for transformations
  if(!is.null(transformations)){
    
    hasTransformations <- TRUE
    fit <- FALSE
    
    if(!is(transformations, "character")) stop("transformations must be a string.")
    
    transformationFunctions <- .compileTransformations(syntax = transformations,
                                                       parameterLabels = parameterLabels)
    
    transformationFunctionPointer <- transformationFunctions$getPtr()
    
    # let's check for new parameters which have to be added to the model:
    addParameters <- transformationFunctions$parameters[!transformationFunctions$parameters %in% parameterTable$label]
    if(length(addParameters) != 0)
      parameterTable <- rbind(
        parameterTable,
        data.frame(
          "label" = addParameters,
          "location" = "transformation", 
          "row" = 1, 
          "col" = 1,
          "value" = 0,
          "rawValue" = 0,
          "isTransformation" = FALSE
        )
      )
    
    parameterTable$isTransformation[parameterTable$label %in% transformationFunctions$isTransformation] <- TRUE
    
    if(all(!is.na(transformationFunctions$startingValues))){
      for(i in 1:length(transformationFunctions$startingValues)){
        whichPar <- parameterTable$label == names(transformationFunctions$startingValues)[i]
        parameterTable$rawValue[whichPar] <- transformationFunctions$startingValues[i]
        parameterTable$value[whichPar] <- transformationFunctions$startingValues[i]
      }
    }
    
  }else{
    hasTransformations <- FALSE
  }
  
  # build model
  
  SEMList$matrices$A <- modelMatrices$Amatrix
  SEMList$matrices$S <- modelMatrices$Smatrix
  SEMList$matrices$M <- modelMatrices$Mvector
  SEMList$matrices$F <- modelMatrices$Fmatrix
  
  SEMList$parameters$label <- parameterTable$label
  SEMList$parameters$location <- parameterTable$location
  SEMList$parameters$row <- parameterTable$row - 1 # c++ starts at 0 
  SEMList$parameters$col <- parameterTable$col - 1 # c++ starts at 0 
  SEMList$parameters$value <- parameterTable$value
  SEMList$parameters$rawValue <- parameterTable$rawValue
  SEMList$parameters$isTransformation <- parameterTable$isTransformation
  
  ## set derivative elements
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
      isVariance <- all(rows == cols)
    }else if(uniqueLocation == "Mvector"){
      positionMatrix <- modelMatrices$Mvector
      positionMatrix[] <- 0
      positionMatrix[cbind(rows,cols)] <- 1
      isVariance <- FALSE
    }else if(uniqueLocation == "transformation"){
      next
    }else{
      stop("unknown location")
    }
    SEMList$DerivativeElements[[p]] <- list(label = p, 
                                         location = uniqueLocation, 
                                         isVariance = isVariance, 
                                         positionMatrix = positionMatrix)
    rm(positionMatrix, isVariance) # avoid any carry over
  }
  
  ## set data
  SEMList$rawData$rawData <- rawData
  SEMList$rawData$manifestNames <- manifestNames
  SEMList$rawData$personInSubset <- internalData$individualMissingPatternID - 1
  
  ## add subsets
  for(s in 1:length(internalData$missingSubsets)){
    SEMList$subsets[[s]] <- list(
                           N = internalData$missingSubsets[[s]]$N,
                           persons = which(internalData$individualMissingPatternID == s)-1,
                           observed = internalData$missingSubsets[[s]]$observed,
                           notMissing = internalData$missingSubsets[[s]]$notMissing,
                           covariance = internalData$missingSubsets[[s]]$covariance,
                           means = internalData$missingSubsets[[s]]$means,
                           rawData = internalData$missingSubsets[[s]]$rawData)
  }
  
  ## initialize C++ model
  
  mySEM <- new(SEMCpp)
  
  mySEM$fill(SEMList)
  
  if(hasTransformations){
    mySEM$addTransformation(transformationFunctionPointer, transformationList)
  }
  
  # the following step is necessary if the parameters of the lavaanModel do not correspond to those in
  # the matrices of the lavaan object. This is, for instance, the case if the starting values
  # instead of the estimates are used.
  parameters <- .getParameters(SEM = mySEM, raw = TRUE)
  mySEM <- .setParameters(SEM = mySEM, 
                          labels = names(parameters), 
                          values = parameters, 
                          raw = TRUE)
  
  if(fit){
    mySEM <- .fit(SEM = mySEM)
    if(whichPars == "est" && checkFit){
      # check model fit
      if(round(mySEM$m2LL - (-2*logLik(lavaanModel)), 4) !=0) 
        stop("Error translating lavaan to internal model representation: Different fit in SEMCpp and lavaan")
    }
  }
  
  return(mySEM)
}

#' .setAMatrix
#' 
#' internal function. Populates the matrix with directed effects in RAM notation
#' 
#' @param model model of class lavaan
#' @param lavaanParameterTable parameter table from lavaan
#' @param nLatent number of latent variables
#' @param nManifest number of manifest variables 
#' @param latentNames names of latent variables
#' @param manifestNames names of manifest variables
#' @keywords internal
.setAMatrix <- function(model, lavaanParameterTable, nLatent, nManifest, latentNames, manifestNames){
  
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

#' .setSMatrix
#' 
#' internal function. Populates the matrix with undirected paths in RAM notation
#' 
#' @param model model of class lavaan
#' @param lavaanParameterTable parameter table from lavaan
#' @param nLatent number of latent variables
#' @param nManifest number of manifest variables 
#' @param latentNames names of latent variables
#' @param manifestNames names of manifest variables
#' @keywords internal
.setSMatrix <- function(model, lavaanParameterTable, nLatent, nManifest, latentNames, manifestNames){
  
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

#' .setMVector
#' 
#' internal function. Populates the vector with means in RAM notation
#' 
#' @param model model of class lavaan
#' @param lavaanParameterTable parameter table from lavaan
#' @param nLatent number of latent variables
#' @param nManifest number of manifest variables 
#' @param latentNames names of latent variables
#' @param manifestNames names of manifest variables
#' @param rawData matrix with raw data
#' @keywords internal
.setMVector <- function(model, lavaanParameterTable, nLatent, nManifest, latentNames, manifestNames, rawData){
  
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

