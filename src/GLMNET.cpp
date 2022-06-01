#include <RcppArmadillo.h>

// [[Rcpp :: depends ( RcppArmadillo )]]

// [[Rcpp::export]]
arma::colvec innerGLMNET(const arma::colvec parameters, 
                         const int N,
                         const arma::colvec subGroupGradient, 
                         const arma::mat subGroupHessian, 
                         const double subGroupLambda,
                         const Rcpp::LogicalVector regularized,
                         const arma::colvec adaptiveLassoWeights,
                         const int maxIter, 
                         const double epsBreak,
                         const bool useMultipleConvergenceCriteria
){
  arma::colvec stepDirection(parameters.n_elem, arma::fill::zeros);
  arma::colvec z(parameters.n_elem, arma::fill::zeros), z_old(parameters.n_elem, arma::fill::zeros);
  arma::colvec newParameters(parameters.n_elem, arma::fill::zeros);
  arma::colvec hessianXdirection, HessTimesZ(subGroupHessian.n_rows, arma::fill::zeros);
  arma::mat HessDiag(subGroupHessian.n_rows, subGroupHessian.n_cols, arma::fill::zeros),
  dp_k(1,1, arma::fill::zeros), d2p_k(1,1, arma::fill::zeros), z_j(1,1, arma::fill::zeros), newParameter(1,1, arma::fill::zeros), zChange(1,1, arma::fill::zeros);
  
  for(int it = 0; it < maxIter; it++){
    
    // reset direction z
    z.fill(arma::fill::zeros); 
    z_old.fill(arma::fill::zeros);
    // iterate over parameters
    for(int p = 0; p < parameters.n_elem; p++){
      
      // compute derivative elements:
      hessianXdirection = subGroupHessian*stepDirection;
      dp_k = subGroupGradient.row(p) + hessianXdirection.row(p);
      d2p_k = subGroupHessian.row(p).col(p);
      
      // if the parameter is regularized:
      if(regularized.at(p)){
        newParameter = parameters.at(p);
        
        // adjust stepDirection for regularized parameters
        if((dp_k(0,0)-subGroupLambda*adaptiveLassoWeights.row(p)(0,0)) >= (d2p_k(0,0)*(newParameter(0,0)+stepDirection.row(p)(0,0)))){
          // condition 1
          z_j = -(dp_k-subGroupLambda*adaptiveLassoWeights.row(p))/(d2p_k);
          z.row(p) = z_j(0,0);
          stepDirection.row(p) += z_j(0,0);
        }else if((dp_k(0,0)+subGroupLambda*adaptiveLassoWeights.row(p)(0,0))<=(d2p_k(0,0)*(newParameter(0,0)+stepDirection.row(p)(0,0)))){
          // condition 2
          z_j = -(dp_k+subGroupLambda*adaptiveLassoWeights.row(p))/(d2p_k);
          z.row(p) = z_j(0,0);
          stepDirection.row(p) = stepDirection.row(p) + z_j(0,0);
        }else{
          // condition 3
          z_j = -(newParameter+stepDirection.row(p));
          z.row(p) = z_j(0,0);
          stepDirection.row(p) = stepDirection.row(p)+z_j(0,0);
        }
      }else{
        // if not regularized: coordinate descent with newton direction
        z_j = -dp_k/d2p_k;
        z.row(p) = z_j(0,0);
        stepDirection.row(p) = stepDirection.row(p)+z_j(0,0);
      }
    }
    
    // check inner stopping criterion:
    HessDiag.diag() = subGroupHessian.diag()/N;
    // scaling hessian with sample size as the hessian is based on the likelihood
    // which directly depends on the sample size -> it would make the
    // convergence criterion harder to hit for larger samples otherwise.
    HessTimesZ = HessDiag*pow(z,2);
    
    zChange = z*arma::trans(z_old);
    if(useMultipleConvergenceCriteria & (zChange(0,0) < epsBreak)){
      break;
    }
    
    if(HessTimesZ.max() < epsBreak){
      break;
    }
    z_old = z;
  }
  return(stepDirection);
}