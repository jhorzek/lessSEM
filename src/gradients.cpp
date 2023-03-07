#include <RcppArmadillo.h>
#include "SEM.h"
#include "multivariateNormal.h"
#include "scores.h"
#include "gradients.h"
#include "impliedDerivatives.h"

arma::rowvec gradientsByGroup(const SEMCpp& SEM, bool raw){
  
  // now that we have initialized all the derivatives of the implied means and
  // covariances, we can start computing the actual gradients
  
  const int numberOfMissingnessPatterns = SEM.data.nGroups;
  const std::vector<std::string>& uniqueParameterLabels = SEM.derivElements.uniqueLabels;
  
  const int nParameters = uniqueParameterLabels.size();
  
  arma::rowvec gradients(nParameters, arma::fill::zeros);
  
  for(int mp = 0; mp < numberOfMissingnessPatterns; mp++){
    
    // initialize group specific gradients for parallel computation:
    ParallelGradients parallelGradients(SEM, mp, raw);
    
    RcppParallel::parallelFor(0, 
                              nParameters, //number of parameters
                              parallelGradients
    );
    
    gradients += parallelGradients.groupGradients;
    
  }
  
  return(gradients);
  
}