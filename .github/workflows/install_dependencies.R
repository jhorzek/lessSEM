# installing all dependencies one by one while also checking if they could be loaded

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

for(pkg in pkgs)
{
  
  print(paste0("Installing package ", pkg))
  
  # check if the package has already been installed
  installed_packages <- installed.packages()
  if(pkg %in% installed_packages)
    next
  
  install.packages(pkg, 
                   dependencies = c("Depends", 
                                    "Imports", 
                                    "LinkingTo"),
                   repos = 'http://cran.rstudio.com/')
  
  # if the package could not be installed: exit
  if(!require(pkg, character.only	= TRUE))
  {
    stop(paste0("Could not install package ", pkg, "."))
    quit(status = 1, save = "no")
  }
  
}
