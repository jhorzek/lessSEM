#' stabilitySelection
#' 
#' Provides rudimentary stability selection for regularized SEM. Stability
#' selection has been proposed by Meinshausen & Bühlmann (2010) and was
#' extended to SEM by Li & Jacobucci (2021). The problem that stabiltiy selection
#' tries to solve is the instability of regularization procedures: Small changes in
#' the data set may result in different parameters being selected. To address
#' this issue, stability selection uses random subsamples from the initial data
#' set and fits models in these subsamples. For each parameter, we can now check
#' how often it is included in the model for a given set of tuning parameters. 
#' Plotting these probabilities can provide an overview over which of the parameters
#' are often removed and which remain in the model most of the time. To get
#' a final selection, a threshold t can be defined: If a parameter is in the model
#' t% of the time, it is retained.
#' 
#' # References
#' 
#' - Li, X., & Jacobucci, R. (2021). Regularized structural equation modeling with 
#' stability selection. Psychological Methods, 27(4), 497–518. https://doi.org/10.1037/met0000389
#' 
#' - Meinshausen, N., & Bühlmann, P. (2010). Stability selection. Journal of the 
#' Royal Statistical Society: Series B (Statistical Methodology), 72(4), 417–473. 
#' https://doi.org/10.1111/j.1467-9868.2010.00740.x
#' 
#' @param  modelSpecification a call to one of the penalty functions in lessSEM. See
#' examples for details
#' @param subsampleSize number of subjects in each subsample. Must be smaller than
#' the number of subjects in the original data set
#' @param numberOfSubsamples number of times the procedure should subsample and
#' recompute the model. According to Meinshausen & Bühlmann (2010), 100 seems to
#' work quite well and is also the default in regsem
#' @param threshold percentage of models, where the parameter should be contained in order
#' to be in the final model
#' @param maxTries fitting models in a subset may fail. maxTries sets the maximal
#' number of subsets to try.
#' @return estimates for each subsample and aggregated percentages for each parameter
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
#' # Stability selection
#' stabSel <- stabilitySelection(
#'   # IMPORTANT: Wrap your call to the penalty function in an rlang::expr-Block:
#'   modelSpecification = 
#'     rlang::expr(
#'       lasso(
#'         # pass the fitted lavaan model
#'         lavaanModel = lavaanModel,
#'         # names of the regularized parameters:
#'         regularized = paste0("l", 6:15),
#'         # in case of lasso and adaptive lasso, we can specify the number of lambda
#'         # values to use. lessSEM will automatically find lambda_max and fit
#'         # models for nLambda values between 0 and lambda_max. For the other
#'         # penalty functions, lambdas must be specified explicitly
#'         nLambdas = 50)
#'     ),
#'   subsampleSize = 80,
#'   numberOfSubsamples = 5, # should be set to a much higher number (e.g., 100)
#'   threshold = 70
#' )
#' stabSel
#' plot(stabSel)
#' @export
stabilitySelection <- function(modelSpecification,
                               subsampleSize,
                               numberOfSubsamples = 100,
                               threshold = 70,
                               maxTries = 10*numberOfSubsamples){
  
  if(!is(modelSpecification, "call"))
    stop("modelSpecification must be a call object. Wrap your call to the penalty",
         "in an expr() - Block, See ?lessSEM::stabilitySelection for an example.")
  
  if((threshold < 0) | (threshold > 100))
    stop("threshold is a percentage value and should be between 0 and 100.")
  
  if(!is.null(modelSpecification$nLambdas)){
    rlang::inform("Running model once to compute lambdas")
    initialFit <- eval(modelSpecification)
    modelSpecification$lambdas <- initialFit@fits$lambda
    modelSpecification$nLambdas <- NULL
  }
  
  dataset <- lavaan::lavInspect(
    object = eval(modelSpecification$lavaanModel),
    what = "data")
  
  N <- nrow(dataset)
  if(subsampleSize >= N)
    stop("subsampleSize must be smaller than the sample size in the data set.")
  
  parLabels <- eval(modelSpecification$lavaanModel) |>
    getLavaanParameters() |>
    names()
  
  parameterEstimates <- matrix(NA, 
                               nrow = 0,
                               ncol = length(parLabels) + 1, 
                               dimnames = list(NULL, 
                                               c("subsample", parLabels)))
  it <- 0
  sucessful <- 0
  pb <- utils::txtProgressBar(min = 0, max = numberOfSubsamples, style = 3)
  
  
  while(it < maxTries){
    it <- it + 1
    if(sucessful >= numberOfSubsamples)
      break
    
    # sample
    subset <- dataset[sample(1:nrow(dataset), subsampleSize),,drop = FALSE]
    
    # replace
    modelSpecificationIt <- modelSpecification
    modelSpecificationIt$modifyModel$dataSet <- subset
    
    # fit
    invisible(
      utils::capture.output(
        fitIt <- try(eval(modelSpecificationIt), 
                     silent = TRUE)
      )
    )
    
    if(is(fitIt, "try-error"))
      next
    
    sucessful <- sucessful + 1
    utils::setTxtProgressBar(pb = pb, 
                             value = sucessful)
    
    parameterEstimates <- rbind(
      parameterEstimates,
      cbind(
        sucessful,
        estimates(fitIt)
      )
    )
    
  }
  
  stabilityPaths <- matrix(0, 
                           nrow = nrow(estimates(fitIt)), 
                           ncol = length(parLabels),
                           dimnames = list(NULL, parLabels))
  for(i in unique(parameterEstimates[,"subsample"])){
    stabilityPaths <- stabilityPaths +
      (parameterEstimates[parameterEstimates[,"subsample"] == i, -1] != 0)
  }
  stabilityPaths <- stabilityPaths/sucessful
  
  results <- new("stabSel",
                 regularized = eval(modelSpecification$regularized),
                 tuningParameters = coef(fitIt)@tuningParameters,
                 stabilityPaths = stabilityPaths,
                 percentSelected = apply(parameterEstimates[,-1],2,function(x) mean(x!=0)),
                 selectedParameters = apply(parameterEstimates[,-1],2,function(x) mean(x!=0) >= threshold/100),
                 settings = list(subsampleSize = subsampleSize,
                                 numberOfSubsamples = numberOfSubsamples,
                                 threshold = threshold,
                                 maxTries = maxTries)
  )
  
  return(results)
  
}