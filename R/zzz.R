#' lessSEM
#'
#' lessSEM (**l**inr **i**mplements **n**on-smooth **r**egsem) 
#' is an R package which provides regularized 
#' structural equation modeling (regularized SEM) 
#' building on lavaan; that is, it makes SEM *leaner*.
#'  As the name suggests, lessSEM is heavily inspired by
#'  the [regsem](https://github.com/Rjacobucci/regsem) package.
#'  
#'  The objectives of lessSEM are (1) to compare exact and approximate optimization of regularized SEM
#'  and (2) to provide optimizers for other SEM packages which can be used with an 
#'  interface similar to optim
#'  
#'  lessSEM also provides experimental functions for 
#'  approximate cross-validation and approximate influence functions.
#'  
#'  **Warning**: The package is relatively new and you may 
#'  find more stable implementations of regularized SEM in the 
#'  R packages [regsem](https://github.com/Rjacobucci/regsem) 
#'  and [lslx](https://github.com/psyphh/lslx). 
#'  
#'  THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#'  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
#'  FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
#'  COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
#'  AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
#'  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
#'
#' @docType package
#' @author Jannik Orzek <orzek@mpib-berlin.mpg.de>
#' @import Rcpp
#' @importFrom Rcpp sourceCpp
#' @useDynLib lessSEM
#' @exportPattern("^[[:alpha:]]+")
#' @name lessSEM
NULL

Rcpp::loadModule("SEM_cpp", TRUE)
Rcpp::loadModule("istaEnet_cpp", TRUE)
Rcpp::loadModule("istaEnetGeneralPurpose_cpp", TRUE)
Rcpp::loadModule("istaCappedL1_cpp", TRUE)
Rcpp::loadModule("istaLSP_cpp", TRUE)
Rcpp::loadModule("istaScad_cpp", TRUE)
Rcpp::loadModule("istaMcp_cpp", TRUE)
Rcpp::loadModule("glmnetEnet_cpp", TRUE)
Rcpp::loadModule("glmnetEnetGeneralPurpose_cpp", TRUE)
Rcpp::loadModule("bfgsEnet_cpp", TRUE)