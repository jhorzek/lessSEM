#' simulateExampleData
#' 
#' simulate data for a simple CFA model
#' @param N number of persons in the data set
#' @param loadings loadings of the latent variable on the manifest observations
#' @examples 
#' y <- lessSEM::simulateExampleData()
#' @export
simulateExampleData <- function(N = 100, # sample size
  loadings = c(rep(1,5), rep(.4,5), rep(0,5))
){
  
  f <- matrix(rnorm(N, 0, 1), ncol = 1) # latent factor
  L <- matrix(loadings, 
    nrow = 1) # loadings
  # covariances
  covs <- diag(max(L^2)+.2, length(loadings))
  
  y <- matrix(NA, nrow = N, ncol = ncol(L))
  
  for(i in 1:N){
    y[i,] <- L*f[i,] +  mvtnorm::rmvnorm(1, sigma = covs)
  }
  
  yNames <- paste0("y", 1:ncol(y))
  colnames(y) <- yNames
  
  return(y)
  
}
