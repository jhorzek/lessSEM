# linr

linr (linr is not [regsem](https://github.com/Rjacobucci/regsem)) is an R package which provides regularized structural equation modeling (regularized SEM) building on lavaan; that is, it makes SEM *leaner*. As the name suggests, linr is heavily inspired by the [regsem](https://github.com/Rjacobucci/regsem) package.

The objectives of linr are:

1. to compare exact and approximate optimization of regularized SEM
2. to provide optimizers for other SEM packages which can be used with an interface similar to optim

linr also also provides experimental functions for approximate cross-validation and approximate influence functions. 

**Warning**: The package is relatively new and you may find more stable implementations of regularized SEM in the R packages [regsem](https://github.com/Rjacobucci/regsem) and [lslx](https://github.com/psyphh/lslx). 

The following features are implemented in linr:

- regularization of SEM with ridge, lasso, adaptive lasso, and elastic net (see ?linr::regularizeSEM)
- automatic selection of $\lambda$ values for lasso and adaptive lasso (see ?linr::regularizeSEM)
- approximate optimization of SEM with custom penalty functions using a BFGS optimizer or Rsolnp (see ?linr::regularizeSEMWithCustomPenalty and ?linr::regularizeSEMWithCustomPenaltyRsolnp)
- approximate cross-validation (see ?linr::aCV4lavaan, ?linr::aCV4RegularizedSEM, and ?linr::aCV4regularizedSEMWithCustomPenalty)
- exact cross-validation for regularized models (see ?linr::CV4regularizedSEM)
- approximate influence functions for regularized SEM (see ?linr::aI4RegularizedSEM)

# Installation

If you want to install linr from GitHub, use the following commands in R:

    if(!require(devtools))install.packages("devtools")

    devtools::install_github("jhorzek/linr")
    

# Example

    library(lavaan) # linr builds on lavaan models
    library(linr)
    
    ## Approximate leave one out cross-validation for a lavaan model
    ### set up model in lavaan
    HS.model <- ' visual  =~ x1 + x2 + x3
                        textual =~ x4 + x5 + x6
                        speed   =~ x7 + x8 + x9 '
    HS <- HolzingerSwineford1939[,paste0("x",1:9)]
    
    fit <- cfa(HS.model, data = HS, meanstructure = TRUE)
    
    ### approximate cross-validation
    aLOOCV <- aCV4lavaan(lavaanModel = fit, k = nrow(HS))
    
    ### we want to compare this to a true leave one out cross-validation
    exactLOOCV <- rep(NA, nrow(HS))
    for(i in 1:nrow(HS)){
      fit = sem(HS.model, HS[-i,], meanstructure = TRUE)
      exactLOOCV[i] <- linr:::computeIndividualM2LL(
        nObservedVariables = ncol(HS), 
        rawData = as.numeric(HS[i,]),
        impliedMeans = fit@implied$mean[[1]], 
        impliedCovariance = fit@implied$cov[[1]])
    }
    
    # The plot shows the relation between exact and approximate cross-validation.
    # If the points are on the line, the approximate and exact cross-validation
    # produce relatively similar results
    plot(exactLOOCV, exactLOOCV, type = "l",
         xlab = "exact loocv", ylab = "approximated loocv")
    points(exactLOOCV, aLOOCV$leaveOutFits, col = "red")
    
    ## Example for regularized SEM
    
    ### Simulate a data set
    set.seed(123)
    dataset <- simulateExampleData()
    
    ### Fit model with lavaan
    modelSyntax <- paste0('f =~ 1*', 
                          colnames(dataset)[1], ' + ', 
                          paste0(colnames(dataset)[2:ncol(dataset)], 
                                 collapse = " + "))
    modelFit = cfa(modelSyntax, dataset, meanstructure = TRUE)
    
    ### optional: plot model
    # semPlot::semPaths(modelFit)
    
    ### Use lasso regularization for a subset of the loadings.
    regularize <- paste0("f=~y", 6:15)
    
    regularizedModel <- regularizeSEM(lavaanModel = modelFit,
                                      penalt = "lasso",
                                      regularizedParameterLabels = regularize,
                                      lambda = seq(0,1,.1)
    )
    
    plot(regularizedModel)
    
    ### approximate leave one out cross-validation:
    aCV <- aCV4regularizedSEM(regularizedModel, k = nrow(dataset))
    plot(aCV)
    
    ### Missing data
    isMissing <- matrix(
      sample(x = c(TRUE,FALSE),
           size = prod(dim(dataset)),
           prob = c(.2,.8), # 20 % missing data
           replace = TRUE),
      nrow = nrow(dataset),
      ncol = ncol(dataset)
      )
    dataset[isMissing] <- NA
    head(dataset)
    
    # fit lavaan model with missing = "ml"
    modelFit = cfa(modelSyntax, dataset, meanstructure = TRUE, missing = "ml")
    regularizedModel <- regularizeSEM(lavaanModel = modelFit,
                                      penalt = "lasso",
                                      regularizedParameterLabels = regularize,
                                      lambda = seq(0,1,.1)
    )
    
    plot(regularizedModel)

# References

* Jacobucci, R., Grimm, K. J., & McArdle, J. J. (2016). Regularized Structural Equation Modeling. Structural Equation Modeling: A Multidisciplinary Journal, 23(4), 555–566. https://doi.org/10.1080/10705511.2016.1154793
* Jacobucci, R., Grimm, K. J., Brandmaier, A. M., Serang, S., Kievit, R. A., & Scharf, F. (2019). regsem: Regularized Structural Equation Modeling. https://CRAN.R-project.org/package=regsem
* Huang, P.-H., Chen, H., & Weng, L.-J. (2017). A Penalized Likelihood Method for Structural Equation Modeling. Psychometrika, 82(2), 329–354. https://doi.org/10.1007/s11336-017-9566-9
* Huang, P.-H. (2020). lslx: Semi-Confirmatory Structural Equation Modeling via Penalized Likelihood. Journal of Statistical Software, 93(7). https://doi.org/10.18637/jss.v093.i07

# Important Notes

THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 