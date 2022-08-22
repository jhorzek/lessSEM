#' .compileTransformations
#' 
#' compile user defined parameter transformations to a 
#' pass to a SEM
#' @param syntax string with user defined transformations
#' @param parameterLabels names of parameters in the model
#' @returns list with parameter names and two Rcpp functions: (1) the transformation function and 
#' (2) a function to create a pointer to the transformation function. 
#' If starting values were defined, these are returned as well.
#' @keywords internal
.compileTransformations <- function(syntax,
                                    parameterLabels){
  
  syntax <- .reduceSyntax(syntax = syntax)
  
  .checkSyntax(syntax = syntax)
  
  parameters <- .extractParametersFromSyntax(syntax = syntax,
                                             parameterLabels = parameterLabels)
  
  armaFunction <- .createRcppTransformationFunction(syntax = syntax, 
                                                    parameters = parameters$parameters)
  
  cat("Compiling the transformation function ... ")
  
  Rcpp::sourceCpp(code = armaFunction)
  cat("done.\n")
  
  return(
    list("parameters" = parameters$parameters,
         "isTransformation" = parameters$parameters[parameters$isTransformation],
         "startingValues" = parameters$startingValues,
         "getPtr" = getPtr,
         "transformationFunction" = transformationFunction     
    )
  )
}

#' .reduceSyntax
#' 
#' reduce user defined parameter transformation syntax to basic elements
#' @param syntax string with user defined transformations
#' @returns a cut and simplified version of the syntax
#' @keywords internal
.reduceSyntax <- function(syntax){
  
  # first, split rows and remove everything we don't need
  # split rows
  syntax <- stringr::str_split(string = syntax, 
                               pattern = "\\n")[[1]]
  # remove white space
  syntax <- stringr::str_replace_all(string = syntax, 
                                     pattern = "\\s",
                                     replacement = "")
  # remove comments
  hasComment <- stringr::str_locate(syntax,
                                    "#|!")
  for(i in 1:length(syntax)){
    if(!is.na(hasComment[i,1])){
      if(hasComment[i,1] == 1){
        syntax[i] <- ""
        next
      }
      syntax[i] <- stringr::str_trunc(
        syntax[i], width = hasComment[i,1]-1, ellipsis = ""
      )
    }
  }
  
  # remove empty
  syntax <- syntax[syntax != ""]
  
  return(syntax)
}

#' .checkSyntax
#' 
#' check the syntax for parameter transformations
#' @param syntax syntax for parameter transformations
#' @return nothing
#' @keywords internal
.checkSyntax <- function(syntax){
  if(!grepl(pattern = "parameters:", syntax[1])) 
    stop("The syntax for parameter transformations must start with a statement similar to paramters: a, b, c, ...")
}

#' .extractParametersFromSyntax
#' 
#' extract the names of the parameters in a syntax
#' @param syntax syntax for parameter transformations
#' @param parameterLabels names of parameters in the model
#' @return vector with names of parameters used in the syntax and vector with
#' boolean indicating if parameter is transformation result
#' @keywords internal
.extractParametersFromSyntax <- function(syntax,
                                         parameterLabels){
  # check which line starts with "parameters:" or "start:"
  parametersAt <- NA
  startingValuesAt <- NA
  for(i in 1:length(syntax)){
    if(grepl("parameters:", syntax[i])){
      parametersAt <- i
    }
    if(grepl("start:", syntax[i])){
      startingValuesAt <- i
    }
  }
  
  if(is.na(parametersAt)) stop("Could not find a statement with 'parameters:' in your transformations")
  
  parameters <- syntax[parametersAt]
  parameters <- gsub(x = parameters, 
                     pattern = "parameters:",
                     replacement = "", 
                     fixed = TRUE)
  parameters <- stringr::str_split(string = parameters, 
                                   pattern = ",")[[1]]
  isTransformation <- rep(FALSE, length(parameters))
  names(isTransformation) <- parameters
  
  if(!is.na(startingValuesAt)){
    startingValues <- syntax[startingValuesAt]
    startingValues <- gsub(x = startingValues, 
                           pattern = "start:",
                           replacement = "", 
                           fixed = TRUE)
    startingValues <- stringr::str_split(string = startingValues, 
                                         pattern = ",")[[1]]
    startingValuesNames <- rep(NA, length(startingValues))
    startingValuesVector <- rep(NA, length(startingValues))
    for(i in 1:length(startingValues)){
      startingValuesNames[i] <- stringr::str_split(string = startingValues[i], 
                                                   pattern = "=")[[1]][1]
      startingValuesVector[i] <- as.numeric(stringr::str_split(string = startingValues[i], 
                                                               pattern = "=")[[1]][2])
    }
    names(startingValuesVector) <- startingValuesNames
  }else{
    startingValuesVector <- NA
  }
  
  for(i in 1:length(syntax)){
    
    if(i == parametersAt) next
    if((!is.na(startingValuesAt)) && i == startingValuesAt) next
    
    isEquation <- grepl(pattern = "=", syntax[i])
    if(isEquation){
      # check left hand side
      lhs <- stringr::str_split(string = syntax[i], 
                                pattern = "=")[[1]][1]
      if(! lhs %in% parameters){
        stop(paste0("Could not find ", lhs, " in parameter: specification of transformations."))
      }
      if(! lhs %in% parameterLabels){
        stop(paste0("Could not find ", lhs, " in model parameters. Every left hand side of your equations must be a model parameter."))
      }
      isTransformation[lhs] <- TRUE
    }
  }
  
  return(list("parameters" = parameters,
              "isTransformation" = isTransformation,
              "startingValuesVector" = startingValuesVector))
  
}

#' .createRcppTransformationFunction
#' 
#' create an Rcpp function which uses the user-defined parameter transformation
#' @param syntax syntax with user defined transformations
#' @param parameters labels of parameters used in these transformations
#' @return string with functions for compilations with Rcpp
#' @keywords internal
.createRcppTransformationFunction <- function(syntax, parameters){
  
  # check which line starts with "parameters:" or "start:"
  parametersAt <- NA
  startingValuesAt <- NA
  for(i in 1:length(syntax)){
    if(grepl("parameters:", syntax[i])){
      parametersAt <- i
    }
    if(grepl("start:", syntax[i])){
      startingValuesAt <- i
    }
  }
  
  # remove parameter and starting values statement from syntax
  
  if(!is.na(startingValuesAt)) 
  {
    syntax <- syntax[-c(parametersAt, startingValuesAt)]
  }else{
    syntax <- syntax[-c(parametersAt)] 
  }
  
  functionHead <- "
  // [[Rcpp::depends(RcppArmadillo)]]
  #include <RcppArmadillo.h>
  // [[Rcpp::export]]
  Rcpp::NumericVector transformationFunction(Rcpp::NumericVector& parameterValues)
  {
  using namespace Rcpp;
  using namespace arma;
  
  // extract required parameters from parameterValues
  "
  functionBody <- c()
  for(p in parameters){
    functionBody <- c(functionBody,
                      paste0('double ', p ,' = parameterValues["', p, '"];') 
    )
  }
  functionBody <- c(functionBody, 
                    "\n\n// add user defined functions",
                    paste0(syntax, ";"), 
                    "\n\n// update parameters"
  )
  for(p in parameters){
    functionBody <- c(functionBody,
                      paste0('parameterValues["', p, '"] = ', p, ';') 
    )
  }
  
  functionEnd <- "
  
  return(parameterValues);
  }"
  
  # we also have to create a function which defines pointers:
  
  ptrFunction <- "
  
  
  // https://gallery.rcpp.org/articles/passing-cpp-function-pointers/
typedef Rcpp::NumericVector (*transformationFunctionPtr)(Rcpp::NumericVector&); //parameters

typedef Rcpp::XPtr<transformationFunctionPtr> transformationFunctionPtr_t;

// [[Rcpp::export]]
transformationFunctionPtr_t getPtr() {
        return(transformationFunctionPtr_t(new transformationFunctionPtr(&transformationFunction)));
}
"

return(
  paste0(c(functionHead, functionBody, functionEnd,
           ptrFunction), collapse = "\n")
)
}
