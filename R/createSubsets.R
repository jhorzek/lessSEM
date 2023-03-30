#' createSubsets
#' 
#' create subsets for cross-validation
#' @param N number of samples in the data set
#' @param k number of subsets to create
#' @returns matrix with subsets
#' @examples 
#' createSubsets(N=100, k = 5)
#' @export
createSubsets <- function(N,k){
  # build subgroups
  if(k < N){
    
    randomCases <- sample(1:N,N)
    
    subsets <- split(randomCases, sort(randomCases%%k)) # Harlan on https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks
    
  }else if(k == N){
    
    subsets <- vector("list",N)
    
    names(subsets) <- 1:N
    
    for(i in 1:N) subsets[[i]] <- i
    
  }else{
    
    stop(paste0("k must be <= ", N))
    
  }
  
  
  subsetMatrix <- matrix(NA, nrow = N, ncol = k,
                         dimnames = list(paste0("person",1:N), 
                                         paste0("testSet", 1:k)))
  for(s in 1:length(subsets)){
    
    subsetMatrix[,s] <- 1:N %in% subsets[[s]]
    
  }
  
  if(any(apply(subsetMatrix,1,sum) != 1)) stop("Error while splitting data in subsets: Some persons are in multiple or none of the subsets")
  
  return(subsetMatrix)
}