#' Class for the coefficients estimated by lessSEM.
#' @slot tuningParameters tuning parameters
#' @slot estimates parameter estimates
#' @slot transformations transformations of parameters
#' @export
setClass("lessSEMCoef",
         representation = representation(
           tuningParameters = "data.frame",
           estimates = "matrix",
           transformations = "matrix"
         ))

#' show
#' 
#' @param object object of class lessSEMCoef
#' @return No return value, just prints estimates
setMethod("show", "lessSEMCoef", function (object) {
  # This function is adapted from the Matrix package. See
  # https://github.com/cran/Matrix/blob/ff84c47d2cb6bca859a841b8d1358ab33ac2a71a/R/sparseMatrix.R#L797.
  
  # Extract estimates and replace zeros:
  res_estimates <- res_estimates_out <- object@estimates
  res_estimates_out[res_estimates == 0] <- "."
  # round to four digits:
  res_estimates_out[res_estimates != 0] <- format(round(res_estimates[res_estimates != 0], 4), 
                                                  nsmall = 4)
  
  if(ncol(object@transformations) != 0){
    # Add transformations and replace zeros:
    res_transformations <- res_transformations_out <- object@transformations
    res_transformations_out[res_transformations == 0] <- "."
    # round to four digits:
    res_transformations_out[res_transformations != 0] <- format(round(res_transformations[res_transformations != 0], 4), 
                                                                nsmall = 4)
  }
  
  # round tuning parameters to 4 digits as well:
  res_tuning_out <- format(round(as.matrix(object@tuningParameters),4), nsmall = 4)
  
  # now, we combine estimates and tuning parameters and separate them by ||--||
  res_out <- cbind(res_tuning_out, 
                   matrix("||--||", 
                          nrow = nrow(object@tuningParameters),
                          ncol = 1, 
                          dimnames = list(NULL, "||--||")),
                   res_estimates_out)
  
  if(ncol(object@transformations) != 0){
    res_out <- cbind(res_out,
                     matrix("||--||", 
                            nrow = nrow(object@tuningParameters),
                            ncol = 1, 
                            dimnames = list(NULL, "||--||")),
                     res_transformations_out
    )
  }
  
  # we extract the column names; they will become their own row
  cnames <- matrix(colnames(res_out), nrow = 1, dimnames = list("", NULL))
  
  # we determine the maximum number of chars in each of the columns to 
  # know how many separators (e.g., =, -), we have to print:
  nchar_tuning <- apply(nchar(
    rbind(colnames(res_tuning_out),
          res_tuning_out)),2,max)
  nchar_tuning[nchar_tuning <= nchar("Tuning")] <- nchar("Tuning") + 1 
  
  nchar_estimates <- apply(nchar(
    rbind(colnames(res_estimates_out),
          res_estimates_out)),2,max)
  nchar_estimates[nchar_estimates <= nchar("Estimates")] <- nchar("Estimates") + 1 
  
  if(ncol(object@transformations) != 0){
    nchar_tranformations <- apply(nchar(
      rbind(colnames(res_transformations_out),
            res_transformations_out)),2,max)
    nchar_tranformations[nchar_tranformations <= nchar("Transform")] <- nchar("Transform") + 1 
  }
  
  # here, we add a == under each parameter name to separate parameter name from parameter value
  if(ncol(object@transformations) != 0){
    res_out <- rbind(
      c(unlist(sapply(nchar_tuning, function(times) paste0(rep(x = "=", times), collapse = ""))), 
        "||--||",
        unlist(sapply(nchar_estimates, function(times) paste0(rep(x = "=", times), collapse = ""))),
        "||--||",
        unlist(sapply(nchar_tranformations, function(times) paste0(rep(x = "=", times), collapse = "")))
      ),
      res_out)
  }else{
    res_out <- rbind(
      c(unlist(sapply(nchar_tuning, function(times) paste0(rep(x = "=", times), collapse = ""))), 
        "||--||",
        unlist(sapply(nchar_estimates, function(times) paste0(rep(x = "=", times), collapse = "")))),
      res_out)
  }
  
  # remove all rownames
  rownames(res_out) <- rep("", nrow(res_out))
  
  # add column names as first row
  if(ncol(object@transformations) != 0){
    res_out <- rbind(
      c("Tuning", rep("", ncol(object@tuningParameters)-1), 
        "||--||", "Estimates", rep("", ncol(object@estimates)-1),
        "||--||", "Transform", rep("", ncol(object@transformations)-1)
      ),
      c(unlist(sapply(nchar_tuning, function(times) paste0(rep(x = "-", times), collapse = ""))), 
        "||--||",
        unlist(sapply(nchar_estimates, function(times) paste0(rep(x = "-", times), collapse = ""))),
        "||--||",
        unlist(sapply(nchar_tranformations, function(times) paste0(rep(x = "-", times), collapse = "")))),
      cnames,
      res_out
    )
  }else{
  res_out <- rbind(
    c("Tuning", rep("", ncol(object@tuningParameters)-1), 
      "||--||", "Estimates", rep("", ncol(object@estimates)-1)
    ),
    c(unlist(sapply(nchar_tuning, function(times) paste0(rep(x = "-", times), collapse = ""))), 
      "||--||",
      unlist(sapply(nchar_estimates, function(times) paste0(rep(x = "-", times), collapse = "")))),
    cnames,
    res_out
  )
  }
  
  # add -- as column names
  colnames(res_out) <- rep("", ncol(res_out))
  print(res_out, 
        quote = FALSE, 
        right = TRUE)
})