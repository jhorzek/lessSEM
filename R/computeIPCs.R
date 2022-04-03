computeIPCs <- function(SEM, individualFitfunction, hessian, raw, stepSize = 1e-5, ...){
  
  parameters <- getParameters(SEM)
  N <- nrow(SEM$data$rawData)
  
  # step 1: compute scores
  scores <- computeScores(SEM = SEM, individualFitfunction = individualFitfunction, stepSize = stepSize, raw = raw, ... = ...)
  
  # step 2: compute IPCs
  IPCs <- scores
  IPCs[] <- NA
  
  hessianInverse <- solve(hessian)

  for(i in 1:nrow(scores)){
    IPCs[i,] <- parameters - N*matrix(scores[i,], nrow = 1)%*%hessianInverse
  }
  
  return(list("IPCs" = IPCs,
              "scores" = scores))
}