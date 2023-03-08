test_that("testing mixed penalty - 2", {
  testthat::skip_on_cran()
  
  library(lessSEM)
  set.seed(123)
  
  dataset <- simulateExampleData()
  
  lavaanSyntax <- "
f =~ l1*y1 + l2*y2 + l3*y3 + l4*y4 + l5*y5 +
     l6*y6 + l7*y7 + l8*y8 + l9*y9 + l10*y10 +
     l11*y11 + l12*y12 + l13*y13 + l14*y14 + l15*y15
f ~~ 1*f
"
  
  lavaanModel <- lavaan::sem(lavaanSyntax,
                             data = dataset,
                             meanstructure = TRUE,
                             std.lv = TRUE)
  
  
  tuningParameters <- expand.grid(lambdaLasso = seq(0,1,length.out=3),
                                  lambdaScad = seq(0,1,length.out=3),
                                  theta = 3.7)
  
  penaltyFunctionArguments <- list(
    eps = 1e-8,
    regularizedLasso = paste0("l", 6:10),
    regularizedScad = paste0("l", 11:15)
  )
  
  # Now, it is time to specify our custom penalty function:
  
  smoothLASSOXSCAD <- function(
    # here are our three arguments:
    parameters,
    tuningParameters,
    penaltyFunctionArguments
  ){
    lambdaLasso <- tuningParameters$lambdaLasso
    lambdaScad <- tuningParameters$lambdaScad
    theta <- tuningParameters$theta
    
    eps <- penaltyFunctionArguments$eps
    regularizedLasso <- penaltyFunctionArguments$regularizedLasso
    regularizedScad <- penaltyFunctionArguments$regularizedScad
    
    penaltyLasso <- lambdaLasso*sum(sqrt(parameters[regularizedLasso]^2 + eps))
    
    
    smoothAbs <- sqrt((parameters[regularizedScad])^2 + eps)
    
    penaltyScad <- 0
    
    for(p in 1:length(smoothAbs)){
      
      if(smoothAbs[p] <= lambdaScad){
        
        penaltyScad <- penaltyScad + lambdaScad * smoothAbs[p]
        
      }else if((lambdaScad < smoothAbs[p]) &&
               (smoothAbs[p] <= lambdaScad * theta)
      ){
        
        penaltyScad <- penaltyScad + (- smoothAbs[p]^2 + 
                                        2 * theta * lambdaScad * 
                                        smoothAbs[p] - 
                                        lambdaScad^2) /(
                                          2*(theta-1)
                                        ) 
        
      }else if(smoothAbs[p] > lambdaScad * theta){
        
        penaltyScad <- penaltyScad + .5*(theta + 1) * lambdaScad^2 
        
      }else{
        stop("impossible scad value")
      }
      
    }
    
    return(penaltyLasso + penaltyScad)
  }
  
  #### Now we are ready to optimize! ####
  regsemApprox <- lessSEM:::.regularizeSEMWithCustomPenaltyRsolnp(lavaanModel = lavaanModel,
                                                                  individualPenaltyFunction = smoothLASSOXSCAD,
                                                                  tuningParameters = tuningParameters,
                                                                  penaltyFunctionArguments = penaltyFunctionArguments)
  
  # let's compare the results to an exact optimization:
  
  exact <- lavaanModel |>
    mixedPenalty() |>
    addLasso(regularized = paste0("l", 6:10),
             lambdas = seq(0,1,length.out=3)) |>
    addScad(regularized = paste0("l", 11:15),
            lambdas = seq(0,1,length.out=3),
            thetas = 3.7) |>
    fit()
  
  
  testthat::expect_true(
    all(abs(exact@parameters[,exact@parameterLabels] -
              regsemApprox@parameters[,exact@parameterLabels]) < .1)
  )
  
  testthat::expect_true(
    all(abs(regsemApprox@fits$regM2LL - 
              exact@fits$regM2LL)/regsemApprox@fits$regM2LL < .01)
  )
  
})