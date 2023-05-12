#include <Rcpp.h>


//' logicalMatch
//' 
//' Returns the rows for which all elements of a boolean matrix X are equal
//' to the elements in boolean vector x
//' 
//' @param X matrix with booleans
//' @param x vector of booleans
//' @return numerical vector with indices of matching rows
// [[Rcpp::export]]
Rcpp::NumericVector logicalMatch(const Rcpp::LogicalMatrix X, Rcpp::LogicalVector x) {
  
  if(X.ncol() != x.length())
    Rcpp::stop("Dimension mismatch");
  
  Rcpp::NumericVector rowsequal;
  bool equal = true;
  
  for(int i = 0; i < X.nrow(); i++){
    for(int j = 0; j < X.ncol(); j++){
      if(X(i,j) != x.at(j)){
        equal = false;
        break;
      }
    }
    if(equal)
      rowsequal.push_back(i+1); // c++ starts with 0
    equal = true;
  }
  return(rowsequal);
}