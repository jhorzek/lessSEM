#ifndef STRINGINVECTOR_H
#define STRINGINVECTOR_H

inline int findStringInVector(std::string what, std::vector<std::string> where, bool throwError){
  
  for (unsigned int i = 0; i < where.size(); i++){
    if(where.at(i).compare(what) == 0){
      return(i);
    }
  }
  if(throwError){
    Rcpp::stop("Could not find the label " + what);
  }else{
    return -1;
  }
}

inline int findStringInVector(std::string what, Rcpp::StringVector where, bool throwError){
  std::string currentString;
  for (int i = 0; i < where.size(); i++){
    currentString = Rcpp::as<std::string>(where.at(i));
    if(currentString.compare(what) == 0){
      return(i);
    }
  }
  if(throwError){
    Rcpp::stop("Could not find the label " + what);
  }else{
    return -1;
  }
}

# endif