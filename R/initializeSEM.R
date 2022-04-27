#' SEMFromLavaan
#' 
#' internal function. Translates an object of class lavaan to the internal model representation.
#' 
#' @param lavaanModel model of class lavaan
#' @param whichPars which parameters should be used to initialize the model. If set to "est", the parameters will be set to the
#' estimated parameters of the lavaan model. If set to "start", the starting values of lavaan will be used. The latter can be useful if parameters are to
#' be optimized afterwards as setting the parameters to "est" may result in the model getting stuck in a local minimum.
#' @param transformVariances set to TRUE to use the internal transformation of variances. This will make sure that estimates for variances can never be negative
#' @param fit should the model be fitted and compared to the lavaanModel?
SEMFromLavaan <- function(lavaanModel, 
                          whichPars = "est",
                          transformVariances = TRUE, 
                          fit = TRUE){
  if(!is(lavaanModel, "lavaan")) stop("lavaanModel must be of class lavaan.")
  
  rawData <- try(lavaan::lavInspect(lavaanModel, "data"))
  if(is(rawData, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  
  # extract basic features
  meanstructure <- lavaanModel@Options$meanstructure
  
  latentNames <- lavNames(lavaanModel, type = "lv")
  manifestNames <- lavNames(lavaanModel, type = "ov")
  
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
  
  internalData <- aCV4SEM:::SEMdata(rawData)
  
  # Extract Model
  
  lavaanParameterTable <- lavInspect(lavaanModel, what = "parTable")
  
  # translate to RAM notation
  
  ## directed paths
  AmatrixElements <- aCV4SEM:::setAMatrix(model = lavaanModel, 
                                          lavaanParameterTable = lavaanParameterTable, 
                                          nLatent = nLatent, 
                                          nManifest = nManifest, 
                                          latentNames = latentNames, 
                                          manifestNames = manifestNames)
  
  ## undirected paths
  SmatrixElements <- aCV4SEM:::setSMatrix(model = lavaanModel, 
                                          lavaanParameterTable = lavaanParameterTable, 
                                          nLatent = nLatent, 
                                          nManifest = nManifest, 
                                          latentNames = latentNames, 
                                          manifestNames = manifestNames)
  
  
  ## Mean structure
  MvectorElements <- aCV4SEM:::setMVector(model = lavaanModel, 
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
    stop(paste0("Could not set the parameters of the model. Set whichPars to one of: 'est', 'start'. See ?aCV4SEM:::SEMFromLavaan for more details."))
  }
  
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
            if(transformVariances && ro==co && matrixName=="Smatrix") {
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
                                "label" = paste0(manifestName, "~1"),
                                "location" = "Mvector", 
                                "row" = which(rownames(MvectorElements$Mvector) == manifestName), 
                                "col" = 1,
                                "value" = MvectorElements$Mvector[manifestName,1],
                                "rawValue" = MvectorElements$Mvector[manifestName,1]
                              )
      )
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
      isVariance <- all(rows == cols)
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
    rm(positionMatrix, isVariance) # avoid any carry over
  }
  
  # set data
  SEMCpp$addRawData(rawData, manifestNames, internalData$individualMissingPatternID-1)
  for(s in 1:length(internalData$missingSubsets)){
    SEMCpp$addSubset(internalData$missingSubsets[[s]]$N,
                     which(internalData$individualMissingPatternID == s)-1,
                     internalData$missingSubsets[[s]]$observed,
                     internalData$missingSubsets[[s]]$notMissing,
                     internalData$missingSubsets[[s]]$covariance,
                     internalData$missingSubsets[[s]]$means,
                     internalData$missingSubsets[[s]]$rawData)
  }
  
  # the following step is necessary if the parameters of the lavaanModel do not correspond to those in
  # the matrices of the lavaan object. This is, for instance, the case if the starting values
  # instead of the estimates are used.
  parameters <- aCV4SEM:::getParameters(SEM = SEMCpp, raw = TRUE)
  SEMCpp <- aCV4SEM:::setParameters(SEM = SEMCpp, labels = names(parameters), values = parameters, raw = TRUE)
  
  if(fit){
    SEMCpp <- aCV4SEM:::fit(SEM = SEMCpp)
    if(whichPars == "est"){
      # check model fit
      if(round(SEMCpp$m2LL - (-2*logLik(lavaanModel)), 4) !=0) stop("Error translating lavaan to internal model representation: Different fit in SEMCpp and lavaan")
    }
  }
  
  return(SEMCpp)
}

#' setAMatrix
#' 
#' internal function. Populates the matrix with directed effects in RAM notation
#' 
#' @param model model of class lavaan
#' @param lavaanParameterTable parameter table from lavaan
#' @param nLatent number of latent variables
#' @param nManifest number of manifest variables 
#' @param latentNames names of latent variables
#' @param manifestNames names of manifest variables
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

#' setSMatrix
#' 
#' internal function. Populates the matrix with undirected paths in RAM notation
#' 
#' @param model model of class lavaan
#' @param lavaanParameterTable parameter table from lavaan
#' @param nLatent number of latent variables
#' @param nManifest number of manifest variables 
#' @param latentNames names of latent variables
#' @param manifestNames names of manifest variables
#' @export
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

#' setMVector
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
#' @export
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

