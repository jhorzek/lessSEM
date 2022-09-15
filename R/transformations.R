#' .compileTransformations
#' 
#' compile user defined parameter transformations to a 
#' pass to a SEM
#' @param syntax string with user defined transformations
#' @param parameterLabels names of parameters in the model
#' @param compile if set to FALSE, the function will not be compiled -> for visual inspection
#' @returns list with parameter names and two Rcpp functions: (1) the transformation function and 
#' (2) a function to create a pointer to the transformation function. 
#' If starting values were defined, these are returned as well.
#' @keywords internal
.compileTransformations <- function(syntax,
                                    parameterLabels,
                                    compile = TRUE){
  
  syntax <- .reduceSyntax(syntax = syntax)
  
  syntax <- .makeSingleLine(syntax = syntax, what = "parameters")
  if(all(is.na(syntax))) stop("Could not find a parameter: statement in your transformations")
  if(!any(is.na(.makeSingleLine(syntax = syntax, what = "start")))){
    syntax <- .makeSingleLine(syntax = syntax, what = "start")
  }
  
  parameters <- .extractParametersFromSyntax(syntax = syntax,
                                             parameterLabels = parameterLabels)
  
  armaFunction <- .createRcppTransformationFunction(syntax = syntax, 
                                                    parameters = parameters$parameters)
  
  if(!compile){
    return(
      list("parameters" = parameters$parameters,
           "isTransformation" = parameters$parameters[parameters$isTransformation],
           "startingValues" = parameters$startingValues,
           "armaFunction" = armaFunction
      )
    )
  }
  
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
  
  # remove empty elements
  syntax <- syntax[!grepl(pattern = "^\\s*$", x = syntax)]
  
  # # check if left hand side of an equation has white space -> will be 
  # # a data type and a variable name
  # isDefinition <- grepl(pattern = "[a-zA-Z:]+\\s+[a-zA-Z:]+\\s*=",
  #                          x = syntax) &
  #   !grepl(pattern = "start\\s*:\\s*",
  #         x = syntax)
  
  return(syntax)
}

#' .makeSingleLine
#' 
#' checks if a parameter: or a start: statement spans multiple lines 
#' and reduces it to one line.
#' 
#' @param syntax reduced syntax
#' @param what which statement to look for (parameters or start)
#' @return a syntax where multi-line statements are condensed to one line
#' @keywords internal
.makeSingleLine <- function(syntax, what){
  
  parameterStatement <- which(grepl(pattern = paste0(what,"\\s*:"), x = syntax))
  if(length(parameterStatement) != 1) return(NA)
  if(grepl(pattern = ",\\s*$", x = syntax[parameterStatement])){
    # is multi-line
    endsAt <- NA
    for(i in parameterStatement:length(syntax)){
      if(!grepl(pattern = ",\\s*$", x = syntax[i])){
        endsAt <- i
        break
      }
    }
    if(is.na(endsAt)) stop(paste0(what, " statement starts but does not end (you may have a trailing ',')."))
  }else{
    endsAt <- parameterStatement
  }
  
  # combine lines
  syntaxNew <- syntax[-c(parameterStatement:(endsAt))]
  syntaxNew <- c(
    paste0(syntax[c(parameterStatement:(endsAt))], collapse = ""),
    syntaxNew
  )
  
  return(syntaxNew)
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
    if(grepl("parameters\\s*:", syntax[i])){
      parametersAt <- i
    }
    if(grepl("start\\s*:", syntax[i])){
      startingValuesAt <- i
    }
  }
  
  if(is.na(parametersAt)) stop("Could not find a statement with 'parameters:' in your transformations")
  
  parameters <- syntax[parametersAt]
  #remove white space
  parameters <- gsub(pattern = "\\s",
                     x = parameters,
                     replacement = "")
  parameters <- gsub(x = parameters, 
                     pattern = "parameters:",
                     replacement = "")
  parameters <- stringr::str_split(string = parameters, 
                                   pattern = ",")[[1]]
  isTransformation <- rep(FALSE, length(parameters))
  names(isTransformation) <- parameters
  
  if(!is.na(startingValuesAt)){
    startingValues <- syntax[startingValuesAt]
    startingValues <- gsub(pattern = "\\s",
                           x = startingValues,
                           replacement = "")
    startingValues <- gsub(x = startingValues, 
                           pattern = "start:",
                           replacement = "")
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
                                pattern = "\\s*=")[[1]][1]
      # remove leading spaces
      lhs <- gsub(pattern = "^\\s*", replacement = "", x = lhs)
      if(lhs %in% parameterLabels){
        isTransformation[lhs] <- TRUE
      }
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
  
  # check for dangling braces
  isDanglingBrace <- grepl(pattern = "^\\s*)\\s*$", x = syntax)
  
  if(any(isDanglingBrace)){
    for(i in rev(which(isDanglingBrace))){
      # going in reverse because there may be layers of braces
      if(i == 0) stop("Your transformations seems to start with a closing brace.")
      syntax[i-1] <- paste0(syntax[i-1], syntax[i])
    }
    syntax <- syntax[!isDanglingBrace]
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
  endsWithOperator <- grepl(pattern = "[+-//*,(=]\\s*$", x = syntax)
  lineEndings <- ifelse(endsWithOperator, "", ";")
  functionBody <- c(functionBody, 
                    "\n\n// add user defined functions",
                    paste0(syntax, lineEndings), 
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
