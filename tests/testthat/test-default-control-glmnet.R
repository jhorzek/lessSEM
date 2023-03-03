test_that("testing default control for glmnet", {
  testthat::skip_on_cran()
tmp <-   try({
  library(Rcpp)
  library(RcppArmadillo)
  library(lessSEM)
  
  funCpp <- '
// [[Rcpp::depends(RcppArmadillo,lessSEM)]]

#include <RcppArmadillo.h>
#include "lessSEM.h"

double sumSquaredError(
    arma::colvec b, // the parameter vector
    arma::colvec y, // the dependent variable
    arma::mat X // the design matrix
){
  // compute the sum of squared errors:
  arma::mat sse = arma::trans(y-X*b)*(y-X*b);
  
  // other packages, such as glmnet, scale the sse with 
  // 1/(2*N), where N is the sample size. We will do that here as well
  
  sse *= 1.0/(2.0 * y.n_elem);
  
  // note: We must return a double, but the sse is a matrix
  // To get a double, just return the single value that is in 
  // this matrix:
  return(sse(0,0));
}

arma::rowvec sumSquaredErrorGradients(
    arma::colvec b, // the parameter vector
    arma::colvec y, // the dependent variable
    arma::mat X // the design matrix
){
  // note: we want to return our gradients as row-vector; therefore,
  // we have to transpose the resulting column-vector:
  arma::rowvec gradients = arma::trans(-2.0*X.t() * y + 2.0*X.t()*X*b);
  
  // other packages, such as glmnet, scale the sse with 
  // 1/(2*N), where N is the sample size. We will do that here as well
  
  gradients *= (.5/y.n_rows);
  
  return(gradients);
}

// THE FOLLOWING CODE IS ADAPTED FROM LAVAAN. 
// SEE lavaan:::lav_model_hessian FOR THE IMPLEMENTATION
// BY Yves Rosseel. The code is under GPL (>= 2)

// [[Rcpp::export]]
arma::mat approximateHessian(arma::colvec b, // the parameter vector
                             arma::colvec y, // the dependent variable
                             arma::mat X, // the design matrix
                             double eps // controls the exactness of the approximation
){
  int nPar = b.n_elem;
  arma::mat hessian(nPar, nPar, arma::fill::zeros);
  
  arma::colvec stepLeft = b, 
    twoStepLeft = b, 
    stepRight = b, 
    twoStepRight = b;
  
  arma::rowvec gradientsStepLeft(nPar);
  arma::rowvec gradientsTwoStepLeft(nPar);
  arma::rowvec gradientsStepRight(nPar);
  arma::rowvec gradientsTwoStepRight(nPar);
  
  for(int p = 0; p < nPar; p++) {
    
    stepLeft.at(p) -= eps;
    twoStepLeft.at(p) -= 2*eps;
    stepRight.at(p) += eps;
    twoStepRight.at(p) += 2*eps;
    
    // step left
    gradientsStepLeft = sumSquaredErrorGradients(stepLeft, y, X);
    
    // two step left
    gradientsTwoStepLeft = sumSquaredErrorGradients(twoStepLeft, y, X);
    
    // step right
    gradientsStepRight = sumSquaredErrorGradients(stepRight, y, X);
    
    // two step right
    gradientsTwoStepRight = sumSquaredErrorGradients(twoStepRight, y, X);
    
    // approximate hessian
    hessian.col(p) = arma::trans((gradientsTwoStepLeft - 
      8.0 * gradientsStepLeft + 
      8.0 * gradientsStepRight - 
      gradientsTwoStepRight)/(12.0 * eps));
    
    // reset
    stepLeft.at(p) += eps;
    twoStepLeft.at(p) += 2*eps;
    stepRight.at(p) -= eps;
    twoStepRight.at(p) -= 2*eps;
  }
  // make symmetric
  hessian = (hessian + arma::trans(hessian))/2.0;
  
  return(hessian);
}

// the first step when linking to lessSEM is to define a model class.
// The procedure is very similar to the ensmallen library, from which
// we have adapted the following approach (see https://ensmallen.org/)
// The model MUST inherit from the model class defined in lessSEM.
// this class lives in the lessSEM namespace and is accessed with lessSEM::model

class linearRegressionModel : public lessSEM::model{
  
public:
  // the lessSEM::model class has two methods: "fit" and "gradients".
  // Both of these methods must follow a fairly strict framework.
  // First: They must receive exactly two arguments: 
  //        1) an arma::rowvec with current parameter values
  //        2) an Rcpp::StringVector with current parameter labels 
  //          (NOTE: the lessSEM package currently does not make use of these labels.
  //                 This is just for future use. If you don"t want to use the labels,
  //                just pass any Rcpp::StringVector you want).
  // Second:
  //        1) fit must return a double (e.g., the -2-log-likelihood) 
  //        2) gradients must return an arma::rowvec with the gradients. It is
  //           important that the gradients are returned in the same order as the 
  //           parameters (i.e., don"t shuffle your gradients, lessSEM will assume
  //           that the first value in gradients corresponds to the derivative with
  //           respect to the first parameter passed to the function).
  
  double fit(arma::rowvec b, Rcpp::StringVector labels) override{
    //NOTE: In sumSquaredError we assumed that b was a column-vector. We
    // have to transpose b to make things work
    return(sumSquaredError(b.t(), y, X));
  }
  
  arma::rowvec gradients(arma::rowvec b, Rcpp::StringVector labels) override{
    //NOTE: In sumSquaredErrorGradients we assumed that b was a column-vector. We
    // have to transpose b to make things work
    return(sumSquaredErrorGradients(b.t(), y, X));
  }
  // IMPORTANT: Note that we used some arguments above which we did not pass to
  // the functions: y, and X. Without these arguments, we cannot use our
  // sumSquaredError and sumSquaredErrorGradients function! To make these accessible
  // to our functions, we have to define them:
  
  const arma::colvec y;
  const arma::mat X;
  
  // finally, we create a constructor for our class
  linearRegressionModel(arma::colvec y_, arma::mat X_):
    y(y_), X(X_){}
  
};

// Step 2: Now that we have defined our model, we can implement the functions
// which are necessary to optimize this model. We will use the elastic net with
// glmnet optimizer as an example.

// [[Rcpp::export]]
Rcpp::List elasticNet(
    const arma::colvec y, 
    arma::mat X,
    const arma::rowvec alpha,
    const arma::rowvec lambda
)
{
  
  // first, let"s add a column to X for our intercept
  X.insert_cols(0,1);
  X.col(0) += 1.0;
  
  // Now we define the parameter vector b
  // This must be an Rcpp::NumericVector with labels
  Rcpp::NumericVector b(X.n_cols,0);
  Rcpp::StringVector bNames;
  // now our vector needs names. This is a bit
  // cumbersome and it would be easier to just let
  // the user pass in a labeled R vector to the elasticNet
  // function. However, we want to make this as convenient
  // as possible for the user. Therefore, we have got to 
  // create these labels:
  for(int i = 0; i < b.length(); i++){
    bNames.push_back("b" + std::to_string(i));
  }
  // add the labels to the parameter vector:
  b.names() = bNames;
  
  // We also have to create a matrix which saves the parameter estimates
  // for all values of alpha and lambda. 
  Rcpp::NumericMatrix B(alpha.n_elem*lambda.n_elem, b.length());
  B.fill(NA_REAL);
  Rcpp::colnames(B) = bNames;
  // we also create a matrix to save the corresponding tuning parameter values
  Rcpp::NumericMatrix tpValues(alpha.n_elem*lambda.n_elem, 2);
  tpValues.fill(NA_REAL);
  Rcpp::colnames(tpValues) = Rcpp::StringVector{"alpha", "lambda"};
  // finally, let"s also return the fitting function value
  Rcpp::NumericVector loss(alpha.n_elem*lambda.n_elem);
  loss.fill(NA_REAL);
  
  // now, it is time to set up the model we defined above
  
  linearRegressionModel linReg(y,X);
  
  // next, we have to define the penalties we want to use.
  // The elastic net is a combination of a ridge penalty and 
  // a lasso penalty:
  lessSEM::penaltyLASSOGlmnet lasso;
  lessSEM::penaltyRidgeGlmnet ridge;
  // these penalties take tuning parameters of class tuningParametersEnetGlmnet
  lessSEM::tuningParametersEnetGlmnet tp;
  
  // finally, there is also the weights. The weights vector indicates, which
  // of the parameters is regularized (weight = 1) and which is unregularized 
  // (weight =0). It also allows for adaptive lasso weights (e.g., weight =.0123).
  // weights must be an arma::rowvec of the same length as our parameter vector.
  arma::rowvec weights(b.length());
  weights.fill(1.0); // we want to regularize all parameters
  weights.at(0) = 0.0; // except for the first one, which is our intercept.
  tp.weights = weights;
  
  // if we want to fine tune the optimizer, we can use the control
  // arguments. We will start with the default control elements and 
  // tweak some arguments to our liking:
  lessSEM::controlGLMNET control = lessSEM::controlGlmnetDefault();
  
  // now it is time to iterate over all lambda and alpha values:
  int it = 0;
  for(int a = 0; a < alpha.n_elem; a++){
    for(int l = 0; l < lambda.n_elem; l++){
      
      // set the tuning parameters
      tp.alpha = alpha.at(a);
      tp.lambda = lambda.at(l);
      
      tpValues(it,0) = alpha.at(a);
      tpValues(it,1) = lambda.at(a);
      
      // to optimize this model, we have to pass it to
      // one of the optimizers in lessSEM. These are 
      // glmnet and ista. We"ll use glmnet in the following. The optimizer will
      // return an object of class fitResults which will have the following fields:
      // - convergence: boolean indicating if the convergence criterion was met (true) or not (false)
      // - fit: double with fit values
      // - fits: a vector with the fits of all iterations
      // - parameterValues the final parameter values as an arma::rowvec
      // - Hessian: the BFGS Hessian approximation
      
      lessSEM::fitResults lmFit = lessSEM::glmnet(
        linReg, // the first argument is our model
        b, // the second are the parameters
        lasso, // the third is our lasso penalty
        ridge, // the fourth our ridge penalty
        tp, // the fifth is our tuning parameter 
        control // finally, let"s fine tune with the control
      );
      
      loss.at(it) = lmFit.fit;
      
      for(int i = 0; i < b.length(); i++){
        // let"s save the parameters
        B(it,i) = lmFit.parameterValues.at(i);
        // and also carry over the current estimates for the next iteration
        b.at(i) = lmFit.parameterValues.at(i);
      }
      
      // carry over Hessian for next iteration
      control.initialHessian = lmFit.Hessian;
      
      it++;
    }
  }
  
  Rcpp::List retList = Rcpp::List::create(
    Rcpp::Named("B") = B,
    Rcpp::Named("tuningParameters") = tpValues,
    Rcpp::Named("loss") = loss);
  return(
    retList
  );
}
'

sourceCpp(code = funCpp, rebuild = TRUE)

set.seed(123)
library(Matrix) # we will use the matrix package
# to print the matrices below in the same sparse format
# as glmnet or ncvreg

# let's first define a small print function to beautify our results
printCoefficients <- function(model){
  print(t(Matrix:::Matrix(model$B, sparse = TRUE)))
}

# first, we simulate data for our
# linear regression.
N <- 100 # number of persons
p <- 10 # number of predictors
X <- matrix(stats::rnorm(N*p),	nrow = N, ncol = p) # design matrix
b <- c(rep(1,4), 
       rep(0,6)) # true regression weights
y <- X%*%matrix(b,ncol = 1) + stats::rnorm(N,0,.2)

# define the tuning parameters
lambda = seq(1,0,length.out = 5)

lasso1 <- elasticNet(y = y,
                     X = X,
                     alpha = 1, # note: glmnet and lessSEM define 
                     # the elastic net differently (lessSEM follows lslx and regsem)
                     # Therefore, you will get different results if you change alpha
                     # when compared to glmnet
                     lambda = lambda
)


# For comparison, we will fit the model with the glmnet package:
library(glmnet)
lassoGlmnet <- glmnet(x = X, 
                      y = y, 
                      lambda = lambda,
                      standardize = FALSE)
testthat::expect_equal(all(abs(coef(lassoGlmnet) - t(lasso1$B)) < 1e-4), TRUE)
}, silent = TRUE)
warning("Not running test-default-control-glmnet. The function requires compilation and only works when running by hand. Make sure to test this function")

testthat::expect_equal(is(tmp, "try-error"), TRUE)

})
