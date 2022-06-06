#' SEM
#' 
#' Fit an unregularized SEM. The function builds on the brilliant lavaan package which 
#' simplifies setting up models considerably. The main objective of SEM is to
#' test the ensmallen optimizers.
#' 
#' @param lavaanSyntax syntax used to generate a model for lavaan
#' @param rawData matrix with data used for the model
#' @param addMeans estimate the means?
#' @param optimizer optimizer to use. SEM relies on optimizers from ensmallen
#' @param fchange convergence criterion value: change in fit
#' @param verbose boolean: Print every tenth iteration?
#' @return object of class SEMCpp
#' @examples 
#' library(lavaan)
#' library(aCV4SEM)
#' set.seed(123)
#' modelSyntax <- ' 
#'   # latent variable definitions
#'      ind60 =~ x1 + x2 + x3
#'      dem60 =~ y1 + a*y2 + b*y3 + c*y4
#'      dem65 =~ y5 + a*y6 + b*y7 + c*y8
#' 
#'   # regressions
#'     dem60 ~ ind60
#'     dem65 ~ ind60 + dem60
#' 
#'   # residual correlations
#'     y1 ~~ y5
#'     y2 ~~ y4 
#'     y3 ~~ y7
#'     y4 ~~ y8
#'     y6 ~~ y8
#' '
#' 
#' SEM <- aCV4SEM:::SEM(lavaanSyntax = modelSyntax,
#'                      rawData = PoliticalDemocracy)

SEM <- function(lavaanSyntax,
                rawData,
                optimizer = "lbfgs",
                fchange = 1e-8,
                verbose = FALSE){
  if(verbose) cat("Setting up model...\n")
  lavaanModel <- sem(model = modelSyntax, 
                     data = rawData, 
                     meanstructure = TRUE,
                     missing = "ml",
                     do.fit = FALSE)
  SEMCpp <- aCV4SEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                    whichPars = "start", 
                                    fit = FALSE
                                    )
  parameterValues <- aCV4SEM:::getParameters(SEMCpp, raw = TRUE)
  parameterLabels <- names(parameterValues)
  
  if(verbose) cat("Starting optimization...\n")
  optimizedParameters <- c(SEMCpp$optimize(matrix(parameterValues, ncol = 1), 
                                           parameterLabels, 
                                           optimizer,
                                           fchange,
                                           verbose))
  
  SEMCpp <- aCV4SEM:::setParameters(SEM = SEMCpp, 
                                    values = optimizedParameters,
                                    labels = parameterLabels, 
                                    raw = TRUE)
  SEMCpp <- aCV4SEM:::fit(SEMCpp)
  if(verbose) cat("Done.\n")
  return(SEMCpp)
  
}