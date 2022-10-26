FROM rocker/tidyverse

RUN R -e "install.packages(c('Rcpp', 'RcppArmadillo', 'numDeriv', 'knitr', 'plotly', 'rmarkdown', 'Rsolnp', 'testthat', 'glmnet', 'ncvreg', 'regsem', 'lslx', 'mvtnorm', 'Matrix', 'OpenMx', 'ctsemOMX'), dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN apt-get update && apt-get install -y qpdf
