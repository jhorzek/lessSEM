#' Class for the coefficients estimated by lessSEM.
#' @slot tuningParameters tuning parameters
#' @slot estimates parameter estimates
setClass("lessSEMCoef",
         representation = representation(
           tuningParameters = "data.frame",
           estimates = "matrix"
         ))

#' show
#' 
#' @param object object of class lessSEMCoef
setMethod("show", "lessSEMCoef", function (object) {
  # This function is adapted from the Matrix package. See
  # https://github.com/cran/Matrix/blob/ff84c47d2cb6bca859a841b8d1358ab33ac2a71a/R/sparseMatrix.R#L797.
  
  # Extract estimates and replace zeros:
  res_estimates <- res_estimates_out <- object@estimates
  res_estimates_out[res_estimates == 0] <- "."
  # round to four digits:
  res_estimates_out[res_estimates != 0] <- format(round(res_estimates[res_estimates != 0], 4), 
                                                  nsmall = 4)
  
  # round tuning parameters to 4 digits as well:
  res_tuning_out <- format(as.matrix(object@tuningParameters), nsmall = 4)
  
  # we determine the maximum number of chars in each of the columns to 
  # know how many separators (e.g., =, -), we have to print:
  nchar_estimates <- apply(nchar(res_estimates_out),2,max)
  nchar_tuning <- apply(nchar(res_tuning_out),2,max)
  
  # now, we combine estimates and tuning parameters and separate them by ||--||
  res_out <-cbind(format(as.matrix(object@tuningParameters), nsmall = 4), 
                  matrix("||--||", nrow = nrow(object@tuningParameters), ncol = 1, dimnames = list(NULL, "||--||")),
                  res_estimates_out)
  
  # we extract the column names; they will become their own row
  cnames <- matrix(colnames(res_out), nrow = 1, dimnames = list("", NULL))
  
  # the following will be used to print an additional header:
  header <- paste0(
    paste0("\n  Tuning", paste0(rep(" ", sum(nchar_tuning) - nchar("Tuning")), collapse = "")),
    " ||--|| ", " Estimates\n")
  
  # here, we add a == under each parameter name to separate parameter name from parameter value
  res_out <- rbind(
    c(unlist(sapply(nchar_tuning, function(times) paste0(rep(x = "=", times), collapse = ""))), 
      "||--||",
      unlist(sapply(nchar_estimates, function(times) paste0(rep(x = "=", times), collapse = "")))),
    res_out)
  
  # remove all rownames
  rownames(res_out) <- rep("", nrow(res_out))
  
  # add colum names as first row
  res_out <- rbind(cnames,
                   res_out
                   )
  
  # add -- as column names
  colnames(res_out) <- c(
    unlist(sapply(nchar_tuning, function(times) paste0(rep(x = "-", times), collapse = ""))), 
                         "||--||",
                         unlist(sapply(nchar_estimates, function(times) paste0(rep(x = "-", times), collapse = ""))))
  cat(header)
  print(res_out, 
        quote = FALSE, 
        right = TRUE)
})