#' .compileTransformations
#' 
#' compile user defined parameter transformations to a 
#' pass to a SEM
#' @param syntax string with user defined transformations
#' @param parameterLabels names of parameters in the model
#' @param compile if set to FALSE, the function will not be compiled -> for visual inspection
#' @param notes option to pass a notes to function. All notes of the current
#' function will be added
#' @returns list with parameter names and two Rcpp functions: (1) the transformation function and 
#' (2) a function to create a pointer to the transformation function. 
#' If starting values were defined, these are returned as well.
#' @importFrom Rcpp sourceCpp
#' @keywords internal
.compileTransformations <- function(syntax,
                                    parameterLabels,
                                    compile = TRUE,
                                    notes = NULL){
  
  if(is.null(notes)){
    printNotes <- TRUE
    notes <- c("Notes:")
  }else{
    # notes already exists and we only append the new ones
    printNotes <- FALSE
  }
  
  
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
  
  # We always print the following note:
  cat("\n")
  rlang::inform(
    c("Compilation Note:",
      paste0("Compiling the transformation function with armadillo version ",
             paste0(  RcppArmadillo::armadillo_version(single = FALSE), collapse = "."),
             ". This may take a few seconds.")
    )
  )
  
  Rcpp::sourceCpp(code = armaFunction)
  
  # The following two definitions are there for R cmd check. The functions
  # getPtr and transformationFunction are created by sourceCpp
  if(!exists("getPtr")){
    getPtr <- function(){
      stop("Failed to compile.")
    }
  }
  
  if(!exists("transformationFunction")){
    transformationFunction <- function(){
      stop("Failed to compile.")
    }
  }
  
  notes <- unique(notes)
  if(printNotes & (length(notes) > 1)){
    cat("\n")
    rlang::inform(notes)
  }
  
  return(
    list("parameters" = parameters$parameters,
         "isTransformation" = parameters$parameters[parameters$isTransformation],
         "startingValues" = parameters$startingValues,
         "getPtr" = getPtr, # this function is created when compiling the C++ code.
         "transformationFunction" = transformationFunction, # this function is created when compiling the C++ code.
         "notes" = notes
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
  hasComment <- grepl(pattern = "#",
                      x = syntax)
  if(any(hasComment)){
    cat("\n")
    rlang::inform(c("Note","Found a # in your transformations. Did you want to write a comment? Please use the C++ comment syntax (e.g., \\\\ my comment)"))
  }
  
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
  
  parameterStatement <- which(grepl(pattern = paste0("^\\s*",what,"\\s*:"), x = syntax))
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
    if(grepl("^\\s*parameters\\s*:", syntax[i])){
      parametersAt <- i
    }
    if(grepl("^\\s*start\\s*:", syntax[i])){
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
    if(grepl("^\\s*parameters:", syntax[i])){
      parametersAt <- i
    }
    if(grepl("^\\s*start:", syntax[i])){
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
  
  # check for missing semicolons
  missingSemicolon <- grepl(pattern = "[\\)a-zA-Z0-9_]$", x = syntax) &
    !grepl(pattern = "^\\s*\\/\\/", x = syntax)
  for(ms in which(missingSemicolon)){
    cat("\n")
    rlang::inform(c("Note",paste0("Found the following statement:\n  > ", syntax[ms], "\nDid you forget a semicolon?")))
  }
  
  functionHead <- "
  // [[Rcpp::depends(RcppArmadillo)]]
  #include <RcppArmadillo.h>
  // [[Rcpp::export]]
  Rcpp::NumericVector transformationFunction(Rcpp::NumericVector& parameterValues, Rcpp::List transformationList)
  {
  using namespace Rcpp;
  using namespace arma;
  
  // extract required parameters from parameterValues
  "
  functionBody <- c()
  for(p in parameters){
    if(any(grepl(pattern = p, x = syntax)))
      functionBody <- c(functionBody,
                        paste0('double ', p ,' = parameterValues["', p, '"];') 
      )
  }
  functionBody <- c(functionBody, 
                    "\n\n// add user defined functions",
                    syntax, 
                    "\n\n// update parameters"
  )
  for(p in parameters){
    if(any(grepl(pattern = p, x = syntax)))
      functionBody <- c(functionBody,
                        paste0('parameterValues["', p, '"] = ', p, ';') 
      )
  }
  
  functionEnd <- "
  
  return(parameterValues);
  }"
  
  # we also have to create a function which defines pointers:
  
  ptrFunction <- "
  
  
  // Dirk Eddelbuettel at
  // https://gallery.rcpp.org/articles/passing-cpp-function-pointers/
typedef Rcpp::NumericVector (*transformationFunctionPtr)(Rcpp::NumericVector&, //parameters
Rcpp::List // transformationList
);

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
