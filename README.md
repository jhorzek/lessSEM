# lessSEM

lessSEM (**l**essSEM **es**timates **s**parse **SEM**) is an R package which provides regularized structural equation modeling (regularized SEM) with non-smooth penalty functions (e.g., lasso) building on lavaan. lessSEM is heavily inspired by the [regsem](https://github.com/Rjacobucci/regsem) package which provides similar functionality.

The objectives of lessSEM are:

1. to compare exact and approximate optimization of regularized SEM
2. to provide optimizers for other SEM packages which can be used with an interface similar to optim

lessSEM also provides experimental functions for approximate cross-validation and approximate influence functions. 

**Warning**: The package is relatively new and you may find more stable implementations of regularized SEM in the R packages [regsem](https://github.com/Rjacobucci/regsem) and [lslx](https://github.com/psyphh/lslx). 

The following features are implemented in lessSEM:

- regularization of SEM with ridge, lasso, adaptive lasso, and elastic net (see ?lessSEM::ridge, ?lessSEM::lasso, ?lessSEM::adaptiveLasso, ?lessSEM::elasticNet)
- automatic selection of $\lambda$ values for lasso and adaptive lasso (see ?lessSEM::lasso, ?lessSEM::adaptiveLasso)
- approximate optimization of SEM with custom penalty functions using a BFGS optimizer or Rsolnp (see ?lessSEM::smoothLasso, ?lessSEM::smoothAdaptiveLasso, ?lessSEM::smoothElasticNet, ?lessSEM::regularizeSEMWithCustomPenaltyRsolnp)
- exact cross-validation for regularized models (see ?lessSEM::cv4ridge, ?lessSEM::cv4lasso, ?lessSEM::cv4adaptiveLasso, ?lessSEM::cv4elasticNet)
- UNDER CONSTRUCTION: approximate cross-validation (see ?lessSEM::acv4elasticNet and ?lessSEM::aCV4regularizedSEMWithCustomPenalty)
- UNDER CONSTRUCTION: approximate influence functions for regularized SEM (see ?lessSEM::ai4elasticNet)

Currently, lessSEM has the following optimizers:

- (variants of) iterative shrinkage and thresholding 
- glmnet

These are also available for other packages. There are two ways to implement them:

1. using the R interface: (e.g., ?lessSEM::gpLasso, ?lessSEM::gpAdaptiveLasso, ?lessSEM::gpElasticNet). This interface is similar to the optim optimizers in R
2. All optimizers are implemented as C++ header-only files in lessSEM. Thus, they can be accessed from other packages using C++. The documentation for this approach will follow soon.

# Installation

If you want to install lessSEM from GitHub, use the following commands in R:

    if(!require(devtools))install.packages("devtools")
  
    devtools::install_github("jhorzek/lessSEM")
    

# Example

    library(lessSEM)
    
    # Identical to regsem, lessSEM builds on the lavaan
    # package for model specification. The first step
    # therefore is to implement the model in lavaan.
    
    dataset <- simulateExampleData()
    
    lavaanSyntax <- "
    f =~ l1*y1 + l2*y2 + l3*y3 + l4*y4 + l5*y5 + 
         l6*y6 + l7*y7 + l8*y8 + l9*y9 + l10*y10 + 
         l11*y11 + l12*y12 + l13*y13 + l14*y14 + l15*y15
    f ~~ 1*f
    "
    
    lavaanModel <- lavaan::sem(lavaanSyntax,
                               data = dataset,
                               meanstructure = TRUE,
                               std.lv = TRUE)
    
    # Optional: Plot the model
    # semPlot::semPaths(lavaanModel, 
    #                   what = "est",
    #                   fade = FALSE)
    
    regsem <- lasso(
      # pass the fitted lavaan model
      lavaanModel = lavaanModel,
      # names of the regularized parameters:
      regularized = paste0("l", 6:15),
      # in case of lasso and adaptive lasso, we can specify the number of lambda
      # values to use. lessSEM will automatically find lambda_max and fit
      # models for nLambda values between 0 and lambda_max. For the other
      # penalty functions, lambdas must be specified explicitly
      nLambdas = 50)
    
    # use the plot-function to plot the regularized parameters:
    plot(regsem)
    
    # elements of regsem can be accessed with the @ operator:
    regsem@parameters[1,]
    
    # AIC and BIC:
    AIC(regsem)
    BIC(regsem)
    
    # The best parameters can also be extracted with:
    coef(regsem, criterion = "AIC")
    coef(regsem, criterion = "BIC")
    
    #### Advanced ###
    # Switching the optimizer # 
    # Use the "method" argument to switch the optimizer. The control argument
    # must also be changed to the corresponding function:
    regsemGlmnet <- lasso(
      lavaanModel = lavaanModel,
      regularized = paste0("l", 6:15),
      nLambdas = 50,
      method = "glmnet",
      control = controlGlmnet())
    
    # Note: The results are basically identical:
    regsemGlmnet@parameters - regsem@parameters

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