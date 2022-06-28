#' ai4elasticNet
#' 
#' Provides an approximate influence for models of class regularizedSEM. These models can be fit with ai4elasticNet() (see ?linr::ai4elasticNet)
#' in this package.
#' 
#' @param regularizedSEM model of class regularizedSEM
#' @param k the number of subset fold. We recommend leave-one-out influence functions; i.e. set k to the number of persons in the data set. Alternatively, 
#' a matrix with pre-defined subsets can be passed to the function. See ?linr::aCV4regularizedSEM for an example
#' @param recomputeHessian if set to false, the Hessians from the quasi newton optimization with GLMNET will be used. Otherwise the Hessian will be recomputed.
#' @param control parameters passed to the GLMNET optimizer. Note that only arguments of the inner iteration are used. See ?controlGLMNET for more details
#' @examples
#' library(linr)
#' 
#' # Let's first set up a regularized model. The following steps are
#' # explained in detail in ?linr::regularizeSEM
#' dataset <- simulateExampleData()
#' 
#' lavaanSyntax <- "
#' f =~ l1*y1 + l2*y2 + l3*y3 + l4*y4 + l5*y5 + 
#'      l6*y6 + l7*y7 + l8*y8 + l9*y9 + l10*y10 + 
#'      l11*y11 + l12*y12 + l13*y13 + l14*y14 + l15*y15
#' f ~~ 1*f
#' "
#' 
#' lavaanModel <- lavaan::sem(lavaanSyntax,
#'                            data = dataset,
#'                            meanstructure = TRUE,
#'                            std.lv = TRUE)
#' 
#' # Optional: Plot the model
#' # semPlot::semPaths(lavaanModel, 
#' #                   what = "est",
#' #                   fade = FALSE)
#' 
#' regsem <- regularizeSEM(
#'   lavaanModel = lavaanModel,
#'   regularizedParameterLabels = paste0("l", 6:15),
#'   penalty = "lasso",
#'   nLambdas = 5)
#' plot(regsem)
#' 
#' ## The approximate influence function can be computed with:
#' aI <- aI4regularizedSEM(regularizedSEM = regsem,
#'                           # we compute the influence of each individual
#'                           k = nrow(dataset)
#'                           )
#' 
#' # let's plot the influence
#' plot(aI) # The points are the fits for the model
#' # when one individual is removed. The lines are the
#' # fits with all individuals.
#' 
#' # To get more information, install the plotly - package and use::
#' # plot(aI, interactive  = TRUE)
#' # This creates an interactive plot which can be explored with your mouse
#' @export
ai4elasticNet <- function(regularizedSEM, 
                              k,
                              recomputeHessian = TRUE, 
                              control = controlGlmnet()){
  if(!is(regularizedSEM, "regularizedSEM")){
    stop("regularizedSEM must be of class regularizedSEM")
  }
  if(!is(control, "controlGlmnet")){
    stop("control must be of class controlGlmnet These objects can be generated with the controlGlmnet function. See ?linr::controlGlmnet")
  }
  
  data <- try(lavaan::lavInspect(regularizedSEM@inputArguments$lavaanModel, "data"))
  if(is(data, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  
  N <- nrow(data)
  
  aCVSEM <- linr:::SEMFromLavaan(lavaanModel = regularizedSEM@inputArguments$lavaanModel, 
                                 transformVariances = TRUE, 
                                 fit = FALSE)
  
  # create subsets 
  if(is.matrix(k)){
    subsets <- k
    if(!is.logical(subsets)) stop("k must be a matrix with booleans (TRUE, FALSE)")
    if(nrow(subsets) != N) stop(paste0("k must have as many rows as there are subjects in your data set (", N, ")."))
    k <- ncol(subsets)
  }else{
    subsets <- linr:::createSubsets(N = N, k = k)
  }
  
  # extract elements for easier access
  fits <- regularizedSEM@fits
  parameters <- regularizedSEM@parameters
  parameterLabels <- regularizedSEM@parameterLabels
  
  weights <- regularizedSEM@inputArguments$weights
  
  tuningParameters <- data.frame(lambda = fits$lambda, alpha = fits$alpha)
  
  subsetElements<- expand.grid(removedSubset = 1:k, 
                               lambda = unique(tuningParameters$lambda), 
                               alpha = unique(tuningParameters$alpha)
  )
  subsetFits <- cbind(subsetElements, 
                      matrix(NA, nrow(subsetElements), 
                             ncol = 1, 
                             dimnames = list(NULL, "fit")))
  
  subsetParameters <- cbind(subsetElements, 
                            matrix(NA, nrow(subsetElements), 
                                   ncol = length(parameterLabels), 
                                   dimnames = list(NULL, parameterLabels)))
  
  
  progressbar = txtProgressBar(min = 0,
                               max = nrow(tuningParameters), 
                               initial = 0, 
                               style = 3)
  
  for(ro in 1:nrow(tuningParameters)){
    
    lambda <- tuningParameters$lambda[ro]
    alpha <- tuningParameters$alpha[ro]
    pars <- coef(regularizedSEM, lambda = lambda, alpha = alpha)
    
    aCVSEM <- linr:::setParameters(SEM = aCVSEM,
                                   labels = names(pars),
                                   values = pars,
                                   raw = FALSE)
    aCVSEM <- linr:::fit(aCVSEM)
    
    if(recomputeHessian){
      hessianOfDifferentiablePart <- NULL
    }else{
      if(!regularizedSEM@inputArguments$control$saveHessian) 
        stop("Hessians were not saved in the regularizedSEM object. This is the default as saving the Hessians can take a lot of disk space. To save the Hessians, use the controlGLMNET argument (see ?controlGLMNET).")
      select <- regularizedSEM@internalOptimization$HessiansOfDifferentiablePart$lambda == lambda &
        regularizedSEM@internalOptimization$HessiansOfDifferentiablePart$alpha == alpha
      hessianOfDifferentiablePart <- regularizedSEM@internalOptimization$HessiansOfDifferentiablePart$Hessian[[which(select)]]
    }
    
    aInfluence <- try(linr:::acv4enet_GLMNET_SEMCpp(SEM = aCVSEM, 
                                                    subsets = subsets,
                                                    raw = TRUE, 
                                                    weights = weights, 
                                                    lambda = lambda,
                                                    alpha = alpha,
                                                    hessianOfDifferentiablePart = hessianOfDifferentiablePart,
                                                    control = control))
    subsetFits$fit[subsetFits$lambda == lambda & subsetFits$alpha == alpha] <- 
      aInfluence$leaveOutFits
    subsetParameters[subsetFits$lambda == lambda & subsetFits$alpha == alpha, parameterLabels] <- 
      aInfluence$subsetParameters[,parameterLabels]
    
    setTxtProgressBar(progressbar,ro)
    
  }
  
  return(
    new("approximateInfluence",
        subsets=subsets,
        subsetParameters = subsetParameters,
        tuningParameters = tuningParameters,
        subsetFits = subsetFits,
        parameters = regularizedSEM@parameters,
        parameterLabels = parameterLabels,
        regularized = regularizedParameterLabels)
  )
}