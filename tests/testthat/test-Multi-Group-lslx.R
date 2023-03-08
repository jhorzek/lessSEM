test_that("multiplication works", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("lslx")
  library(lslx)
  library(lessSEM)
  set.seed(123)
  HS <- lavaan::HolzingerSwineford1939[complete.cases(lavaan::HolzingerSwineford1939),]
  
  lambdas <- seq(0,1,length.out=5)
  
  ## EXAMPLE: Multi-Group Factor Analysis ##
  # run `vignette("multi-group-analysis")` to see the vignette
  # specify multi-group factor analysis model
  model_mgfa <- "visual  :=> 1 * x1 + x2 + x3
               textual :=> 1 * x4 + x5 + x6
               speed   :=> 1 * x7 + x8 + x9"
  
  # "school" is set as group variable and "Pasteur" is specified as reference
  lslx_mgfa <- lslx$new(model = model_mgfa,
                        data = HS,
                        group_variable = "school",
                        reference_group = "Pasteur")
  
  # penalize increment components of loadings and intercepts in 'Grant-White'
  lslx_mgfa$penalize_heterogeneity(block = c("y<-1", "y<-f"), 
                                   group = "Grant-White")
  
  # fit with lasso penalty
  lslx_mgfa$fit_lasso(lambda_grid = lambdas)
  
  #### Fit the same model with lessSEM ####
  
  model <- "
visual  =~ l11_g * x1 + l21_g*x2 + l31_g*x3
textual =~ l42_g * x4 + l52_g*x5 + l62_g*x6
speed   =~ l73_g * x7 + l83_g*x8 + l93_g*x9

visual ~~ v11_g*visual + v12_g*textual + v13_g*speed
textual ~~ v22_g*textual + v23_g*speed
speed ~~ v33_g*speed

x1 ~~ mv11_g*x1
x2 ~~ mv22_g*x2
x3 ~~ mv33_g*x3
x4 ~~ mv44_g*x4
x5 ~~ mv55_g*x5
x6 ~~ mv66_g*x6
x7 ~~ mv77_g*x7
x8 ~~ mv88_g*x8
x9 ~~ mv99_g*x9

x1 ~ m1_g*1
x2 ~ m2_g*1
x3 ~ m3_g*1
x4 ~ m4_g*1
x5 ~ m5_g*1
x6 ~ m6_g*1
x7 ~ m7_g*1
x8 ~ m8_g*1
x9 ~ m9_g*1
"
  
  model_pasteur <- gsub(pattern = "_g", replacement = "_pasteur", x = model)
  
  fit_pasteur <- sem(model = model_pasteur,
                     data = HS[HS$school=="Pasteur",],
                     std.lv = FALSE
  )
  
  model_grand_white <- gsub(pattern = "_g", replacement = "_grand_white", x = model)
  
  fit_grand_white <- sem(model = model_grand_white,
                         data = HS[HS$school=="Grant-White",],
                         std.lv = FALSE
  )
  
  pasteur_coef <- unique(names(coef(fit_pasteur)))
  # don't regularize differences in variances:
  regularize <- ! (grepl(pattern = "^mv[0-9]*_", x = pasteur_coef) |
                     grepl(pattern = "^v[0-9]*_", x = pasteur_coef))
  pasteur_coef <- pasteur_coef[regularize]
  grand_white_coef <- gsub(pattern = "_pasteur", replacement = "_grand_white", x = pasteur_coef)
  delta_coef <- gsub(pattern = "_pasteur", replacement = "_delta", x = pasteur_coef)
  
  
  transformations <- paste0("parameters: ",paste0(c(pasteur_coef, grand_white_coef, delta_coef), collapse = ", "),"\n",
                            paste0(paste0(grand_white_coef, " = ", pasteur_coef, " + ", delta_coef, ";"), collapse = "\n")
  )
  cat(transformations)
  
  lsem <- lasso(lavaanModel = c(fit_grand_white, fit_pasteur), 
                regularized = delta_coef, 
                lambdas = lambdas, 
                modifyModel = modifyModel(transformations = transformations))
  
  #### Compare lslx and lessSEM ####
  
  ## for fit comparison:
  N <- nrow(HS)
  nvar <- 9
  saturated_pasteur <- -2*mvtnorm::dmvnorm(HS[HS$school == "Pasteur",paste0("x",1:9)], 
                                           apply(HS[HS$school == "Pasteur",paste0("x",1:9)],2,mean), 
                                           ((N-1)/N)* cov(HS[HS$school == "Pasteur",paste0("x",1:9)]), 
                                           log = TRUE)
  saturated_grand_white <- -2*mvtnorm::dmvnorm(HS[HS$school == "Grant-White",paste0("x",1:9)], 
                                               apply(HS[HS$school == "Grant-White",paste0("x",1:9)],2,mean), 
                                               ((N-1)/N)* cov(HS[HS$school == "Grant-White",paste0("x",1:9)]), 
                                               log = TRUE)
  
  lslx_regM2LL <- rep(NA, length(lambdas))
  
  for(l in 1:length(lambdas)){
    coef_lsem_l <- coef(lsem)@estimates[l,]
    coef_lslx_l <- lslx_mgfa$extract_coefficient(lambda = lambdas[l])
    
    # we have to reconstruct the fit measure used by lessSEM which is based
    # on the -2log-Likelihood
    lslx_objective_value <- lslx_mgfa$extract_numerical_condition(lambda = lambdas[l], 
                                                                  include_faulty = TRUE)["objective_value"]
    
    # Check differences in fit:
    
    lslx_regM2LL[l] <- N*lslx_objective_value + sum(saturated_pasteur) + sum(saturated_grand_white)
    
    testthat::expect_equal((lslx_regM2LL[l] - lsem@fits$regM2LL[l] > 0) & # lessSEM fits better
                             (lslx_regM2LL[l] - lsem@fits$regM2LL[l] < .2) # no large difference
                           , TRUE)
    
    # check loadings
    loadings_lsem_pasteur <- coef_lsem_l[grepl(pattern = "l[0-9][0-9]_pasteur", x = names(coef_lsem_l))]
    # rename to match lslx names:
    names(loadings_lsem_pasteur) <- names(loadings_lsem_pasteur) |>
      gsub(pattern = "^l", replacement = "x") |>
      gsub(pattern = "1_", replacement = "<-visual") |>
      gsub(pattern = "2_", replacement = "<-textual") |>
      gsub(pattern = "3_", replacement = "<-speed") |>
      gsub(pattern = "pasteur", replacement = "/Pasteur")
    
    
    loadings_lslx_pasteur <- coef_lslx_l[grepl(pattern = "x[0-9]<-[a-z]*/Pasteur", x = names(coef_lslx_l))]
    is_free <- loadings_lslx_pasteur != 1
    loadings_lslx_pasteur <- loadings_lslx_pasteur[is_free]
    
    # check difference
    testthat::expect_equal(max(abs(loadings_lslx_pasteur - loadings_lsem_pasteur[names(loadings_lslx_pasteur)])) < .03, TRUE)
    
    ### grand white
    loadings_lsem_delta <- coef_lsem_l[grepl(pattern = "l[0-9][0-9]_delta", x = names(coef_lsem_l))]
    # rename to match lslx names:
    names(loadings_lsem_delta) <- names(loadings_lsem_delta) |>
      gsub(pattern = "^l", replacement = "x") |>
      gsub(pattern = "1_", replacement = "<-visual") |>
      gsub(pattern = "2_", replacement = "<-textual") |>
      gsub(pattern = "3_", replacement = "<-speed") |>
      gsub(pattern = "delta", replacement = "/Grant-White")
    
    
    loadings_lslx_delta <- coef_lslx_l[grepl(pattern = "x[0-9]<-[a-z]*/Grant-White", x = names(coef_lslx_l))]
    loadings_lslx_delta <- loadings_lslx_delta[is_free]
    
    # check difference
    testthat::expect_equal(max(abs(loadings_lslx_delta - loadings_lsem_delta[names(loadings_lslx_delta)])) < .03, TRUE)
  }
  
})
