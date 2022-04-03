computeExpected <- function(SEM){
  
  SEM$model$expected$covariance[] <- computeRAMExpectedCovariance(Fmatrix = SEM$model$matrices$Fmatrix, 
                                                          Amatrix = SEM$model$matrices$Amatrix, 
                                                          Smatrix = SEM$model$matrices$Smatrix);
  SEM$model$expected$means[] <- computeRAMExpectedMeans(Fmatrix = SEM$model$matrices$Fmatrix, 
                                                Amatrix = SEM$model$matrices$Amatrix,
                                                Mvector = SEM$model$matrices$Mvector);
  
  return(SEM)
}

