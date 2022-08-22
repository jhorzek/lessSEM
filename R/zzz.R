#' lessSEM
#'
#' lessSEM (**l**essSEM **e**stimates **s**parse **SEMr**) 
#' is an R package which provides regularized 
#' structural equation modeling (regularized SEM) 
#' building on lavaan.
#' 
#' 
#' lessSEM is heavily inspired by the [regsem](https://github.com/Rjacobucci/regsem) 
#' and [lslx](https://github.com/psyphh/lslx) packages. The objectives of lessSEM are (1) to compare exact and approximate optimization of regularized SEM
#' and (2) to provide optimizers for other SEM packages which can be used with an 
#' interface similar to optim.
#'  
#' There are multiple vignettes for specific topics that may be interesting to you:
#'  
#' - If you are new to **lessSEM**, check out the vignette 'lessSEM' 
#' (`vignette("lessSEM", package = "lessSEM")`).
#' - Some advanced regularization procedures can be achieved by using model transformations 
#' (e.g., regularizing differences between parameters). See 
#' `vignette("Parameter-transformations", package = "lessSEM")` for more details.
#' - `vignette("General-Purpose-Optimization", package = "lessSEM")`  explains 
#' the use of the optimizers implemented in **lessSEM** for models other than SEMs
#' - `vignette("The-Structural-Equation-Model", package = "lessSEM")`  explains 
#' the how SEMs are implemented in **lessSEM**
#' - `vignette("The-optimizer-interface", package = "lessSEM")`  explains 
#' the details of the underlying optimizer framework
#'   
#' **Warning**: The package is relatively new and you may 
#' find more stable implementations of regularized SEM in the 
#' R packages [regsem](https://github.com/Rjacobucci/regsem) 
#' and [lslx](https://github.com/psyphh/lslx). 
#'  
#' ## References
#' 
#' More details can be found in the following publications:
#' 
#' ### R - Packages / Software
#' 
#' - [lavaan](https://github.com/yrosseel/lavaan) Rosseel, Y. (2012). lavaan: An R Package for Structural Equation Modeling. Journal of Statistical Software, 48(2), 1–36. https://doi.org/10.18637/jss.v048.i02
#' - [regsem](https://github.com/Rjacobucci/regsem): Jacobucci, R. (2017). regsem: 
#'   Regularized Structural Equation Modeling. ArXiv:1703.08489 [Stat]. http://arxiv.org/abs/1703.08489
#' - [lslx](https://github.com/psyphh/lslx): Huang, P.-H. (2020). lslx: 
#'   Semi-confirmatory structural equation modeling via penalized likelihood. Journal 
#' of Statistical Software, 93(7). https://doi.org/10.18637/jss.v093.i07
#' - [fasta](https://cran.r-project.org/web/packages/fasta/index.html): 
#'   Another implementation of the fista algorithm (Beck & Teboulle, 2009).
#' - [ensmallen](https://ensmallen.org/): Curtin, R. R., Edel, M., Prabhu, R. G., 
#' Basak, S., Lou, Z., & Sanderson, C. (2021). The ensmallen library for ﬂexible 
#' numerical optimization. Journal of Machine Learning Research, 22, 1–6.
#' 
#' ### Regularized Structural Equation Modeling
#' 
#' - Huang, P.-H., Chen, H., & Weng, L.-J. (2017). A Penalized Likelihood Method 
#' for Structural Equation Modeling. Psychometrika, 82(2), 329–354. https://doi.org/10.1007/s11336-017-9566-9
#' - Jacobucci, R., Grimm, K. J., & McArdle, J. J. (2016). Regularized Structural 
#' Equation Modeling. Structural Equation Modeling: A Multidisciplinary Journal, 23(4), 
#' 555–566. https://doi.org/10.1080/10705511.2016.1154793
#' 
#' ### Penalty Functions
#' 
#' - Candès, E. J., Wakin, M. B., & Boyd, S. P. (2008). Enhancing Sparsity by 
#' Reweighted l1 Minimization. Journal of Fourier Analysis and Applications, 14(5–6), 
#' 877–905. https://doi.org/10.1007/s00041-008-9045-x
#' - Fan, J., & Li, R. (2001). Variable selection via nonconcave penalized 
#' likelihood and its oracle properties. Journal of the American Statistical 
#' Association, 96(456), 1348–1360. https://doi.org/10.1198/016214501753382273
#' - Hoerl, A. E., & Kennard, R. W. (1970). Ridge Regression: Biased Estimation 
#' for Nonorthogonal Problems. Technometrics, 12(1), 55–67. https://doi.org/10.1080/00401706.1970.10488634
#' - Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. 
#' Journal of the Royal Statistical Society. Series B (Methodological), 58(1), 267–288.
#' - Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax concave penalty. 
#' The Annals of Statistics, 38(2), 894–942. https://doi.org/10.1214/09-AOS729
#' - Zhang, T. (2010). Analysis of Multi-stage Convex Relaxation for Sparse Regularization. 
#' Journal of Machine Learning Research, 11, 1081–1107.
#' - Zou, H. (2006). The adaptive lasso and its oracle properties. Journal of the 
#' American Statistical Association, 101(476), 1418–1429. https://doi.org/10.1198/016214506000000735
#' - Zou, H., & Hastie, T. (2005). Regularization and variable selection via the 
#' elastic net. Journal of the Royal Statistical Society: Series B, 67(2), 301–320. 
#' https://doi.org/10.1111/j.1467-9868.2005.00503.x
#' 
#' ### Optimizer
#' 
#' #### GLMNET 
#' 
#' - Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization paths for 
#' generalized linear models via coordinate descent. Journal of Statistical 
#' Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' - Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for 
#' l1-regularized logistic regression. The Journal of Machine Learning Research, 
#' 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
#' 
#' #### Variants of ISTA
#' 
#' - Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1), 
#' 183–202. https://doi.org/10.1137/080716542
#' - Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). A general iterative 
#' shrinkage and thresholding algorithm for non-convex regularized optimization problems. 
#' Proceedings of the 30th International Conference on Machine Learning, 28(2)(2), 37–45.
#' - Parikh, N., & Boyd, S. (2013). Proximal Algorithms. Foundations and 
#' Trends in Optimization, 1(3), 123–231.
#'  
#' ## Licence
#' 
#'  THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#'  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
#'  FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
#'  COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
#'  AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
#'  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
#'  
#' @examples 
#' library(lessSEM)
#' 
#' # Identical to regsem, lessSEM builds on the lavaan
#' # package for model specification. The first step
#' # therefore is to implement the model in lavaan.
#' 
#' dataset <- simulateExampleData()
#' 
#' lavaanSyntax <- "
#'       f =~ l1*y1 + l2*y2 + l3*y3 + l4*y4 + l5*y5 + 
#'            l6*y6 + l7*y7 + l8*y8 + l9*y9 + l10*y10 + 
#'            l11*y11 + l12*y12 + l13*y13 + l14*y14 + l15*y15
#'       f ~~ 1*f
#'       "
#' 
#' lavaanModel <- lavaan::sem(lavaanSyntax,
#'                            data = dataset,
#'                            meanstructure = TRUE,
#'                            std.lv = TRUE)
#' 
#' # Optional: Plot the model
#' # semPlot::semPaths(lavaanModel, 
#' #                   what = "est",
#' #                   fade = FALSE)
#' 
#' regsem <- lasso(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   # in case of lasso and adaptive lasso, we can specify the number of lambda
#'   # values to use. lessSEM will automatically find lambda_max and fit
#'   # models for nLambda values between 0 and lambda_max. For the other
#'   # penalty functions, lambdas must be specified explicitly
#'   nLambdas = 50)
#' 
#' # use the plot-function to plot the regularized parameters:
#' plot(regsem)
#' 
#' # elements of regsem can be accessed with the @ operator:
#' regsem@parameters[1,]
#' 
#' # AIC and BIC:
#' AIC(regsem)
#' BIC(regsem)
#' 
#' # The best parameters can also be extracted with:
#' coef(regsem, criterion = "AIC")
#' coef(regsem, criterion = "BIC")
#' 
#' # cross-validation
#' cv <- cvLasso(lavaanModel = lavaanModel,
#'               regularized = paste0("l", 6:15),
#'               lambdas = seq(0,1,.1),
#'               standardize = TRUE)
#' 
#' # get best model according to cross-validation:
#' coef(cv)
#' 
#' #### Advanced ###
#' # Switching the optimizer # 
#' # Use the "method" argument to switch the optimizer. The control argument
#' # must also be changed to the corresponding function:
#' regsemGlmnet <- lasso(
#'   lavaanModel = lavaanModel,
#'   regularized = paste0("l", 6:15),
#'   nLambdas = 50,
#'   method = "glmnet",
#'   control = controlGlmnet())
#' 
#' # Note: The results are basically identical:
#' regsemGlmnet@parameters - regsem@parameters
#'
#' @docType package
#' @author Jannik Orzek <orzek@mpib-berlin.mpg.de>
#' @import Rcpp
#' @importFrom Rcpp sourceCpp
#' @useDynLib lessSEM
#' @exportPattern("^[[:alpha:]]+")
#' @name lessSEM
#' @keywords internal
"_PACKAGE"

# SEM Module
Rcpp::loadModule("SEM_cpp", TRUE)

# SEM optimization
Rcpp::loadModule("istaEnet_cpp", TRUE)
Rcpp::loadModule("istaEnetGeneralPurpose_cpp", TRUE)
Rcpp::loadModule("istaCappedL1_cpp", TRUE)
Rcpp::loadModule("istaLSP_cpp", TRUE)
Rcpp::loadModule("istaScad_cpp", TRUE)
Rcpp::loadModule("istaMcp_cpp", TRUE)
Rcpp::loadModule("glmnetEnet_cpp", TRUE)
Rcpp::loadModule("bfgsEnet_cpp", TRUE)

# General Purpose
Rcpp::loadModule("glmnetEnetGeneralPurpose_cpp", TRUE)
Rcpp::loadModule("istaEnetGeneralPurpose_cpp", TRUE)
Rcpp::loadModule("istaLspGeneralPurpose_cpp", TRUE)
Rcpp::loadModule("istaMcpGeneralPurpose_cpp", TRUE)
Rcpp::loadModule("istaScadGeneralPurpose_cpp", TRUE)
Rcpp::loadModule("istaCappedL1GeneralPurpose_cpp", TRUE)

# general Purpose with Cpp
Rcpp::loadModule("glmnetEnetGeneralPurposeCpp_cpp", TRUE)
Rcpp::loadModule("istaEnetGeneralPurposeCpp_cpp", TRUE)
Rcpp::loadModule("istaLspGeneralPurposeCpp_cpp", TRUE)
Rcpp::loadModule("istaMcpGeneralPurposeCpp_cpp", TRUE)
Rcpp::loadModule("istaScadGeneralPurposeCpp_cpp", TRUE)
Rcpp::loadModule("istaCappedL1GeneralPurposeCpp_cpp", TRUE)