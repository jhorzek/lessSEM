#' regularizeSEMWithCustomPenaltyRsolnp
#' 
#' Optimize a SEM with custom penalty function using the Rsolnp optimizer (see ?Rsolnp::solnp). This optimizer is the default in regsem (see ?regsem::cv_regsem).
#' 
#' @param lavaanModel model of class lavaan 
#' @param individualPenaltyFunction penalty function which takes the current parameter values as first argument, the tuning parameters as second, and the penaltyFunctionArguments as third argument and 
#' returns a single value - the value of the penalty function for a single person. If the true penalty function is non-differentiable (e.g., lasso) a smooth
#' approximation of this function should be provided.
#' @param tuningParameters data.frame with tuning parameter values. Important: The function will iterate over the rows of these tuning parameters and pass them to your penalty function
#' @param penaltyFunctionArguments arguments passed to individualPenaltyFunction, individualPenaltyFunctionGradient, and individualPenaltyFunctionHessian
#' @param startingValues option to provide initial starting values. Only used for the first lambda. Three options are supported. Setting to "est" will use the estimates
#' from the lavaan model object. Setting to "start" will use the starting values of the lavaan model. Finally, a labeled vector with parameter
#' values can be passed to the function which will then be used as starting values.
#' @param carryOverParameters should parameters from the previous iteration be used as starting values of
#' the next iteration?
#' @param control option to set parameters of the optimizer; see ?Rsolnp::solnp
#' @examples 
#' library(lessSEM)
#' 
#' # Identical to regsem, lessSEM builds on the lavaan
#' # package for model specification. The first step
#' # therefore is to implement the model in lavaan.
#' 
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
#' ## Defining a custom penalty function is a bit more complicated than
#' # using the default ones provided in regularizeSEM (see ?lessSEM::regularizeSEM).
#' # We start with the definition of the penalty function. Make sure that the derivatives
#' # of this function are well defined (i.e., the function is smooth)!
#' # We will use the lasso penalty as an example.
#' 
#' # The penalty function MUST accept three arguments; how you use them is up to you.
#' 
#' # The first argument are the parameters. lessSEM will pass a vector with the current
#' # parameter values and their names to your function. The parameter labels will be
#' # exactly the same as those used by lavaan! So you can check them before with:
#' 
#' lessSEM::getLavaanParameters(lavaanModel)
#' 
#' # The vector passed to your function will look like the vector above.
#' 
#' # The second argument is called tuningParameters MUST be a data.frame with the
#' # tuning-parameters. lessSEM will automatically iterate over the rows of this object.
#' # In case of LASSO regularization there is only one tuning parameter: lambda.
#' # Therefore, we specify the tuningParameters object as:
#' 
#' tuningParameters <- data.frame(lambda = seq(0,1,.1)) # we will test 11 lambdas here
#' print(tuningParameters)
#' 
#' # The third argument is called penaltyFunctionArguments and we can pass anything
#' # we want here. For the lasso penalty, we need two additional things:
#' # 1) We need the smoothing-parameter epsilon, which makes sure that our
#' # penalty is differentiable
#' # 2) We need the labels of the regularized parameters, so that we can only
#' # penalize those while all others remain unpenalized.
#' 
#' penaltyFunctionArguments <- list(
#'   eps = 1e-10,
#'   regularizedParameterLabels = paste0("l", 6:15)
#' )
#' 
#' # Now, it is time to specify our custom penalty function:
#' 
#' smoothLASSO <- function(
#'   # here are our three arguments:
#'   parameters,
#'   tuningParameters,
#'   penaltyFunctionArguments
#' ){
#'   # to make it easier to see what is going on:
#'   lambda <- tuningParameters$lambda # tuningParameters will be ONE ROW OF the
#'   # tuningParameters object we created before -> it will contain just
#'   # one lambda in this case!
#'   eps <- penaltyFunctionArguments$eps
#'   regularizedParameterLabels <- penaltyFunctionArguments$regularizedParameterLabels
#' 
#'   regularizedParameters <- parameters[regularizedParameterLabels]
#' 
#'   # now, let's define our penalty function:
#'   penaltyLasso <- lambda*sum(sqrt(regularizedParameters^2 + eps))
#' 
#'   return(penaltyLasso)
#' }
#' # Important: This penalty function is assumed to be for a single individual only.
#' # lessSEM will multiply it with sample size N to get the penalty value of the
#' # full sample!
#' 
#' #### Now we are ready to optimize! ####
#' regsemApprox <- regularizeSEMWithCustomPenaltyRsolnp(lavaanModel = lavaanModel,
#'                                                individualPenaltyFunction = smoothLASSO,
#'                                                tuningParameters = tuningParameters,
#'                                                penaltyFunctionArguments = penaltyFunctionArguments)
#' 
#' # let's compare the results to an exact optimization:
#' 
#' regsemExact <- lasso(
#'   lavaanModel = lavaanModel,
#'   regularized = paste0("l", 6:15),
#'   lambdas = tuningParameters$lambda)
#' 
#' head(regsemExact@parameters[,regsemExact@parameterLabels] -
#'        regsemApprox@parameters[,regsemExact@parameterLabels])
#' # Note that the parameter estimates are basically identical.
#' @export
regularizeSEMWithCustomPenaltyRsolnp <- function(lavaanModel, 
                                                 individualPenaltyFunction,
                                                 tuningParameters,
                                                 penaltyFunctionArguments,
                                                 startingValues = "est",
                                                 carryOverParameters = TRUE,
                                                 control = list("trace" = 0)){
  
  inputArguments <- as.list(environment())
  
  if(!is.data.frame(tuningParameters)) stop("tuningParameters must be a data frame (e.g., data.frame(lambda = seq(0,1,.1))).")
  
  if(!is(lavaanModel, "lavaan")){
    stop("lavaanModel must be of class lavaan")
  }
  
  if(lavaanModel@Options$estimator != "ML") stop("lavaanModel must be fit with ml estimator.")
  
  rawData <- try(lavaan::lavInspect(lavaanModel, "data"))
  if(is(rawData, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  
  sampleSize <- nrow(rawData)
  
  ### initialize model ####
  if(any(startingValues == "est")){
    SEM <- lessSEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                   transformVariances = TRUE,
                                   whichPars = "est",
                                   addMeans = control$addMeans, 
                                   activeSet = control$activeSet)
  }else if(any(startingValues == "start")){
    SEM <- lessSEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                   transformVariances = TRUE,
                                   whichPars = "start",
                                   addMeans = control$addMeans, 
                                   activeSet = control$activeSet)
  }else if(is.numeric(startingValues)){
    
    if(!all(names(startingValues) %in% names(lessSEM::getLavaanParameters(lavaanModel)))) stop("Parameter names of startingValues do not match those of the lavaan object. See lessSEM::getLavaanParameters(lavaanModel).")
    SEM <- lessSEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                   transformVariances = TRUE,
                                   whichPars = "start", 
                                   fit = FALSE,
                                   addMeans = control$addMeans, 
                                   activeSet = control$activeSet)
    SEM <- lessSEM:::setParameters(SEM = SEM, labels = names(startingValues), value = startingValues, raw = FALSE)
    SEM <- try(lessSEM:::fit(SEM))
    if(is(SEM, "try-error") || !is.finite(SEM$m2LL)) stop("Infeasible starting values.")
    
  }else{
    stop("Invalid startingValues passed to regularizeSEMWithCustomPenaltyRsolnp.")
  }
  
  # get parameters
  parameters <- lessSEM:::getParameters(SEM, raw = TRUE)
  
  # define fitting function
  fitfun <- function(parameters, 
                     SEM, 
                     sampleSize,
                     individualPenaltyFunction, 
                     tuningParameters,
                     penaltyFunctionArguments){
    SEM <- try(lessSEM:::setParameters(SEM = SEM, 
                                       labels = names(parameters), 
                                       values = parameters, 
                                       raw = TRUE), 
               silent = TRUE)
    if(is(SEM, "try-error")){
      return(99999999999999999)
    }
    SEM <- try(lessSEM:::fit(SEM), silent = TRUE)
    if(is(SEM, "try-error") || !is.finite(SEM$m2LL)){
      return(99999999999999999)
    }
    
    # check if covariance matrix is positive definite
    if(any(eigen(SEM$S, only.values = TRUE)$values < 0)){
      return(99999999999999999)
    }
    if(any(eigen(SEM$impliedCovariance, only.values = TRUE)$values < 0)){
      return(99999999999999999)
    }
    
    penalty <- sampleSize*individualPenaltyFunction(parameters, 
                                                    tuningParameters, 
                                                    penaltyFunctionArguments)
    if(is(penalty, "try-error") || !is.finite(penalty)){
      return(99999999999999999)
    }
    
    # division by N to be closer to the implementation in regsem
    return((SEM$m2LL + penalty)/sampleSize)
  }
  
  fits <- data.frame(
    "m2LL" = rep(NA,nrow(tuningParameters)),
    "regM2LL"= rep(NA,nrow(tuningParameters)),
    "convergence" = rep(NA,nrow(tuningParameters))
  )
  fits <- cbind(tuningParameters,
                fits)
  parameterEstimates <- as.data.frame(matrix(NA,nrow = nrow(tuningParameters), ncol = length(parameters)))
  colnames(parameterEstimates) <- names(parameters)
  parameterEstimates <- cbind(
    tuningParameters,
    parameterEstimates
  )
  
  Hessians <- list(NULL)
  
  progressbar = utils::txtProgressBar(min = 0, 
                               max = nrow(tuningParameters), 
                               initial = 0, 
                               style = 3)
  
  parametersInit <- parameters
  
  for(i in 1:nrow(tuningParameters)){
    
    utils::setTxtProgressBar(progressbar,i)
    
    currentTuningParameters <- tuningParameters[i,,drop = FALSE]
    
    result <- try(Rsolnp::solnp(pars = parametersInit, fun = fitfun, SEM = SEM, 
                                sampleSize = sampleSize,
                                individualPenaltyFunction = individualPenaltyFunction, 
                                tuningParameters = currentTuningParameters,
                                penaltyFunctionArguments = penaltyFunctionArguments,
                                control = control))
    if(is(result, "try-error")) next
    # regsem does not re-use the parameters for the next iteration;
    # this is considerably slower, but it does help the optimizer.
    # The optimizer builds an approximation of the Hessian which
    # will be better if more iterations have to be done until convergence
    if(carryOverParameters) parametersInit <- result$pars
    
    SEM <- try(lessSEM:::setParameters(SEM = SEM, labels = names(result$pars), 
                                       values = result$pars, raw = TRUE), 
               silent = TRUE)
    SEM <- try(lessSEM:::fit(SEM), silent = TRUE)
    
    parameterEstimates[i, names(parameters)] <- lessSEM:::getParameters(SEM, raw = FALSE)[names(parameters)]
    fits$m2LL[i] <- SEM$m2LL
    fits$regM2LL[i] <- SEM$m2LL + sampleSize*individualPenaltyFunction(result$pars, 
                                                                       currentTuningParameters, 
                                                                       penaltyFunctionArguments)
    #result$value[length(result$value)]
    fits$convergence[i] <- result$convergence == 0
    
    if(result$convergence != 0) warning(paste0("Optimizer did not converge for ", paste0(names(currentTuningParameters), " = ",currentTuningParameters) , collapse = ", "))
  }
  
  internalOptimization <- list(
    "HessiansOfDifferentiablePart" = Hessians
  )
  
  results <- new("regularizedSEMWithCustomPenalty",
                 parameters = parameterEstimates,
                 fits = fits,
                 parameterLabels = names(parameters),
                 internalOptimization = internalOptimization,
                 inputArguments = inputArguments)
  
  return(results)
  
}
