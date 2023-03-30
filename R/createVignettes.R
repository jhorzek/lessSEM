#' .knitVignettes
#' 
#' Takes vignettes of format .lessmd and knits them to .Rmd files to be used
#' as vignettes. The Reason for this two-step approach is to reduce the runtime
#' on CRAN. The function is adapted from Stefan Kloppenborg at 
#' https://www.kloppenborg.ca/2021/06/long-running-vignettes/
#' 
#' @param dir directory, where the vignettes are located.
#' @return creates Rmd vignettes
.knitVignettes <- function(dir = "vignettes"){
  
  if(!("knitr" %in% rownames(utils::installed.packages())))
    stop("Requires knitr.")
  
  
  pkgs <- c('lavaan',
            'Rcpp',
            'RcppArmadillo',
            'RcppParallel',
            'ggplot2',
            'tidyr',
            'stringr',
            'methods',
            'numDeriv',
            'utils',
            'stats',
            'graphics',
            'knitr',
            'plotly',
            'rmarkdown',
            'Rsolnp',
            'testthat',
            'glmnet',
            'ncvreg',
            'regsem',
            'lslx',
            'mvtnorm',
            'Matrix',
            'OpenMx',
            'ctsemOMX')
  
  for(p in pkgs){
    if(!requireNamespace(package = p))
      stop("Package ", p, " required to build vignettes.")
  }
  
  currentDir <- getwd()
  on.exit(setwd(currentDir))
  
  setwd(dir)
  
  files <- list.files()
  files <- files[grepl(pattern = ".lessmd$", 
                       x = files)]
  
  for(f in files){
    cat("Knitting", f, "\n")
    
    # knit file to markdown
    outFile <- gsub(pattern = ".lessmd$", 
                    replacement = "", 
                    x = f)
    knitr::knit(input = paste0(f), 
                output = paste0(outFile, ".Rmd"))
  }
}