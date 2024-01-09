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
#' @param transformationGradientStepSize step size used to compute the gradients of the
#' transformations
#' @returns Object of class Rcpp_SEMCpp
#' @keywords internal
.SEMFromLavaan <- function(lavaanModel, 
                           whichPars = "est",
                           fit = TRUE,
                           addMeans = TRUE,
                           activeSet = NULL,
                           dataSet = NULL,
                           transformations = NULL,
                           transformationList = list(),
                           transformationGradientStepSize = 1e-6){
  
  SEMObj <- .extractSEMFromLavaan(lavaanModel = lavaanModel,
                                  whichPars = whichPars,
                                  fit = fit,
                                  addMeans = addMeans,
                                  activeSet = activeSet,
                                  dataSet = dataSet,
                                  transformations = transformations)
  
  ## initialize C++ model
  
  mySEM <- new(SEMCpp)
  
  mySEM$fill(SEMObj$SEMList)
  
  if(!is.null(transformations)){
    mySEM$addTransformation(SEMObj$transformationFunctionPointer, transformationList)
    mySEM$setTransformationGradientStepSize(transformationGradientStepSize)
  }
  
  # the following step is necessary if the parameters of the lavaanModel do not correspond to those in
  # the matrices of the lavaan object. This is, for instance, the case if the starting values
  # instead of the estimates are used.
  parameters <- .getParameters(SEM = mySEM, raw = TRUE)
  mySEM <- .setParameters(SEM = mySEM, 
                          labels = names(parameters), 
                          values = parameters, 
                          raw = TRUE)
  
  if(SEMObj$fit){
    mySEM <- .fit(SEM = mySEM)
    
    if(whichPars == "est" && SEMObj$checkFit && lavaanModel@Options$do.fit){
      # check model fit
      if(mySEM$getEstimator() == "fiml"){
        if(round(mySEM$objectiveValue - (-2*lavaan::logLik(lavaanModel)), 4) !=0) 
          stop("Error translating lavaan to internal model representation: Different fit in SEMCpp and lavaan")
      }else if(mySEM$getEstimator() == "wls"){
        # to stay consistent with lavaan, we have to use the following
        # scaling:
        lessSEMObjective <- 0.5 * (mySEM$sampleSize-1)/(mySEM$sampleSize^2)* mySEM$objectiveValue
        if(round(lessSEMObjective - (lavaan::fitMeasures(lavaanModel, "fmin")), 4) !=0)
          stop("Error translating lavaan to internal model representation: Different fit in SEMCpp and lavaan")
        
      }
    }
  }
  
  return(mySEM)
}

#' .multiGroupSEMFromLavaan
#' 
#' internal function. Translates a vector of objects of class lavaan to the 
#' internal model representation.
#' 
#' @param lavaanModels vector with lavaan models
#' @param whichPars which parameters should be used to initialize the model. If set to "est", the parameters will be set to the
#' estimated parameters of the lavaan model. If set to "start", the starting values of lavaan will be used. The latter can be useful if parameters are to
#' be optimized afterwards as setting the parameters to "est" may result in the model getting stuck in a local minimum.
#' @param fit should the model be fitted
#' @param addMeans If lavaanModel has meanstructure = FALSE, addMeans = TRUE will add a mean structure. FALSE will set the means of the observed variables to the average
#' @param transformations string with transformations
#' @param transformationList list for transformations
#' @param transformationGradientStepSize step size used to compute the gradients of the
#' transformations
#' @returns Object of class Rcpp_mgSEMCpp
#' @keywords internal
.multiGroupSEMFromLavaan <- function(lavaanModels, 
                                     whichPars = "est",
                                     fit = TRUE,
                                     addMeans = TRUE,
                                     transformations = NULL,
                                     transformationList = list(),
                                     transformationGradientStepSize = 1e-6){
  
  SEMs <- vector("list", length(lavaanModels))
  
  for(m in 1:length(lavaanModels)){
    
    if(!lavaanModels[[m]]@Options$meanstructure && addMeans){
      warning("Model ", m, " does not have an explicit mean structure. If the items in the ",
              "submodels have the same names, lessSEM will assume that the means for ",
              "the different groups are also the same!")
    }
    # Note: transformation will be defined for the full model, not the submodels.
    SEMs[[m]] <- .extractSEMFromLavaan(lavaanModel = lavaanModels[[m]],
                                       whichPars = whichPars,
                                       fit = FALSE,
                                       addMeans = addMeans,
                                       activeSet = NULL,
                                       dataSet = NULL,
                                       transformations = NULL)$SEMList
    
  }
  
  ## initialize C++ model
  
  mySEM <- new(mgSEM, length(lavaanModels))
  mySEM$addModels(SEMs)
  # if(length(lavaanModels) > 10){
  #   cat("Thats a lot of models... Setup may take some time as all R objects are copied to C++. Will print progress:\n")
  #   it <- 0
  #   pb <- utils::txtProgressBar(min = it, max = length(lavaanModels), style = 3)
  # }
  # 
  # for(m in SEMs){
  #   if(length(lavaanModels) > 10){
  #     it <- it + 1
  #     utils::setTxtProgressBar(pb, it)
  #   }
  #   mySEM$addModel(m)
  # }
  
  # extract parameters
  parameters <- .getParameters(SEM = mySEM, raw = TRUE)
  
  if(!is.null(transformations)){
    transforms <- .createMultiGroupTransformations(transformations = transformations, 
                                                   parameterValues = parameters)
    
    mySEM$addTransformation(transforms$parameterValues, 
                            transforms$isTransformation,
                            transforms$transformationFunctionPointer,
                            transformationList)
    
    mySEM$setTransformationGradientStepSize(transformationGradientStepSize)
  }
  
  # the following step is necessary if the parameters of the lavaanModel do not correspond to those in
  # the matrices of the lavaan object. This is, for instance, the case if the starting values
  # instead of the estimates are used.
  mySEM <- .setParameters(SEM = mySEM, 
                          labels = names(parameters), 
                          values = parameters, 
                          raw = TRUE)
  
  return(mySEM)
}

#' .extractSEMFromLavaan
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
#' @param transformations optional: transform parameter values.
#' @returns list with SEMList (model in RAM representation) and fit (boolean indicating if the model should
#' be fit and compared to lavaan)
#' @keywords internal
.extractSEMFromLavaan <- function(lavaanModel,
                                  whichPars = "est",
                                  fit = TRUE,
                                  addMeans = TRUE,
                                  activeSet = NULL,
                                  dataSet = NULL,
                                  transformations = NULL
){
  if(!is(lavaanModel, "lavaan")) stop("lavaanModel must be of class lavaan.")
  
  SEMList <- list(
    estimator = NULL,
    WLSWeightsInverse = NULL,
    WLSMeanstructure = NULL,
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
  
  if(tolower(lavaanModel@Options$estimator) %in% c("ml", "fiml","mlm", "mlmv", "mlmvs", "mlf", "mlr")){
    SEMList$estimator <- "fiml"
  }else if(tolower(lavaanModel@Options$estimator) %in% c("uls","wls", "dwls", "gls")){
    SEMList$estimator <- "wls"
    SEMList$WLSMeanstructure <- lavaanModel@Options$meanstructure
  }else{
    stop("Currenlty only models estimated with maximum likelihood (estimator = 'ml') or WLS (estimator = 'WLS') are supported.")
  }
  
  rawDataList <- .getRawData(lavaanModel = lavaanModel, dataSet = dataSet, estimator = SEMList$estimator)
  rawData <- rawDataList$rawData
  checkFit <- rawDataList$checkFit
  
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
  
  if(SEMList$estimator == "fiml"){
    internalData <- .SEMdata(rawData)
  }else if(SEMList$estimator == "wls"){
    internalData <- .SEMdataWLS(rawData, lavaanModel)
    SEMList$WLSWeightsInverse <- lavInspect(object = lavaanModel, what = "wls.v")
  }else{
    stop("Unknown estimator")
  }
  
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
  Fmatrix <- .setFmatrix(nManifest = nManifest, 
                         manifestNames = manifestNames, 
                         nLatent = nLatent, 
                         latentNames = latentNames)
  
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
  parameterTable <- .createParameterTable(parameterValues = parameterValues,
                                          parameterLabels = parameterLabels, 
                                          modelParameters = modelParameters, 
                                          parameterIDs = parameterIDs)
  
  # if no mean structure: add 
  if(!meanstructure && addMeans) {
    
    if((!lavaanModel@Options$meanstructure) & (SEMList$estimator == "wls"))
      stop("Cannot post-hoc add a meanstructure to a lavaan model fitted with WLS. Please add the meanstructure to the lavaan model using meanstructure = TRUE.")
    
    parameterTable <- .addMeanStructure(parameterTable = parameterTable, 
                                        manifestNames = manifestNames, 
                                        MvectorElements = MvectorElements)
  }
  
  # add transformations
  
  if(!is.null(transformations)){ 
    hasTransformations <- TRUE
    fit <- FALSE
    if(!is(transformations, "character")) stop("transformations must be a string.")
    
    transformsList <- .createTransformations(transformations = transformations,
                                             parameterLabels = parameterLabels,
                                             parameterTable = parameterTable
    )
    
    transformationFunctionPointer <- transformsList$transformationFunctionPointer
    parameterTable <- transformsList$parameterTable
    
  }else{
    transformationFunctionPointer <- NULL
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
  SEMList <- .defineDerivatives(SEMList = SEMList,
                                parameterTable = parameterTable, 
                                modelMatrices = modelMatrices)
  
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
  
  return(list(SEMList = SEMList,
              transformationFunctionPointer = transformationFunctionPointer,
              fit = fit, 
              checkFit = checkFit))
}

#' .getRawData
#' 
#' Extracts the raw data from lavaan or adapts a user supplied data set to
#' the structure of the lavaan data
#' @param lavaanModel model fitted with lavaan
#' @param dataSet user supplied data set
#' @param estimator which estimator is used?
#' @return raw data
#' @keywords internal
.getRawData <- function(lavaanModel, dataSet, estimator){
  
  if(!is.null(dataSet) && (estimator != "fiml")){
    stop("Providing a different data set is currently only supported for full information maximum likelihood estimation.")
  }
  
  if(is.null(dataSet)){
    
    rawData <- try(lavaan::lavInspect(lavaanModel, "data"))
    if(is(rawData, "try-error")) 
      stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")

    checkFit <- TRUE
    
  }else{
    
    lavaanData <- try(lavaan::lavInspect(lavaanModel, "data"))
    # make sure that the sorting of rawData is correct:
    rawData <- dataSet[,colnames(lavaanData),drop = FALSE]
    if(any(!colnames(lavaanData) %in% colnames(rawData))) 
      stop("Not all variables are in rawData.")
    checkFit <- FALSE
    
  }

  # remove empty rows:
  if(any(apply(rawData, 1, function(x) all(is.na(x))))){
      warning("Your data set has rows where all observations are missing. lessSEM will ",
              "remove those rows, but it is recommended to do so before fitting the models.")
      rawData <- rawData[!apply(rawData, 1, function(x) all(is.na(x))),,drop = FALSE]
  }

  return(list(rawData = rawData, checkFit = checkFit))
  
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

#' .setFmatrix
#' 
#' returns the filter matrix of a RAM
#' 
#' @param nManifest number of manifest variables
#' @param manifestNames names of manifest variables
#' @param nLatent number of latent variables
#' @param latentNames names of latent variables
#' @return matrix
#' @keywords internal
.setFmatrix <- function(nManifest, manifestNames, nLatent, latentNames){
  Fmatrix <- diag(nManifest)
  rownames(Fmatrix) <- manifestNames
  colnames(Fmatrix) <- manifestNames
  Fmatrix <- cbind(
    matrix(0, nrow = nManifest, ncol = nLatent, dimnames = list(manifestNames, 
                                                                latentNames)),
    Fmatrix)
  return(Fmatrix)
}

#' .addMeanStructure
#' 
#' adds a mean strucuture to the parameter table
#' @param parameterTable table with parameters
#' @param manifestNames names of manifest variables
#' @param MvectorElements elements of the means vector
#' @return parameterTable
#' @keywords internal
.addMeanStructure <- function(parameterTable, manifestNames, MvectorElements){
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
  return(parameterTable)
}


#' .createParameterTable
#' 
#' create a parameter table using the elements extracted from lavaan
#' 
#' @param parameterValues values of parameters
#' @param parameterLabels names of the parameters
#' @param modelParameters model parameters from lavaan
#' @param parameterIDs unique parameter ids from lavaan -> identify each parameter 
#' with a unique number
#' @return parameter table for lessSEM
#' @keywords internal 
.createParameterTable <- function(parameterValues, parameterLabels, modelParameters, parameterIDs){
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
  
  return(parameterTable)
}

#' .defineDerivatives
#' 
#' adds all elements required to compute the derivatives of the fitting function
#' with respect to the parameters to the SEMList
#' 
#' @param SEMList list representing SEM
#' @param parameterTable table with parameters
#' @param modelMatrices matrices of the RAM model
#' @return SEMList
#' @keywords internal
.defineDerivatives <- function(SEMList, parameterTable, modelMatrices){
  
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
  
  return(SEMList)
}
