test_that("Multicore works", {
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
  
  # lasso
  # Note: Glmnet uses random update orders for the parameters -> we reset the
  # seed here to make sure that the orders are the same for single and multi-core
  # execution.
  set.seed(123)
  regsem1 <- lasso(
    lavaanModel = lavaanModel,
    regularized = paste0("l", 6:15),
    nLambdas = 5,
    control = controlGlmnet(nCores = 1))
  set.seed(123)
  regsem2 <- lasso(
    lavaanModel = lavaanModel,
    regularized = paste0("l", 6:15),
    nLambdas = 5,
    control = controlGlmnet(nCores = 2))
  
  testthat::expect_equal(all(abs(regsem1@parameters - regsem2@parameters) < 1e-8), TRUE)
  
  
  # adaptive lasso
  # Note: Glmnet uses random update orders for the parameters -> we reset the
  # seed here to make sure that the orders are the same for single and multi-core
  # execution.
  set.seed(123)
  regsem1 <- adaptiveLasso(
    lavaanModel = lavaanModel,
    regularized = paste0("l", 6:15),
    nLambdas = 5,
    control = controlGlmnet(nCores = 1))
  set.seed(123)
  regsem2 <- adaptiveLasso(
    lavaanModel = lavaanModel,
    regularized = paste0("l", 6:15),
    nLambdas = 5,
    control = controlGlmnet(nCores = 2))
  
  testthat::expect_equal(all(abs(regsem1@parameters - regsem2@parameters) < 1e-8), TRUE)
  
  # scad
  # Note: The update order is not random -> no need to reset the seed
  regsem1 <- scad(
    lavaanModel = lavaanModel,
    regularized = paste0("l", 6:15),
    lambdas = seq(0,.3,length.out = 3),
    thetas = c(2.1,2.7,3.5),
    control = controlIsta(nCores = 1))
  regsem2 <- scad(
    lavaanModel = lavaanModel,
    regularized = paste0("l", 6:15),
    lambdas = seq(0,.3,length.out = 3),
    thetas = c(2.1,2.7,3.5),
    control = controlIsta(nCores = 2))
  
  testthat::expect_equal(all(abs(regsem1@parameters - regsem2@parameters) < 1e-8), TRUE)
  
  # mcp
  # Note: The update order is not random -> no need to reset the seed
  regsem1 <- mcp(
    lavaanModel = lavaanModel,
    regularized = paste0("l", 6:15),
    lambdas = seq(0,.3,length.out = 3),
    thetas = c(2.1,2.7,3.5),
    control = controlIsta(nCores = 1))
  regsem2 <- mcp(
    lavaanModel = lavaanModel,
    regularized = paste0("l", 6:15),
    lambdas = seq(0,.3,length.out = 3),
    thetas = c(2.1,2.7,3.5),
    control = controlIsta(nCores = 2))
  
  testthat::expect_equal(all(abs(regsem1@parameters - regsem2@parameters) < 1e-8), TRUE)
  
  
  # cross-validation
  
  # Note: Glmnet uses random update orders for the parameters -> we reset the
  # seed here to make sure that the orders are the same for single and multi-core
  # execution.
  set.seed(123)
  regsem1 <- cvLasso(
    lavaanModel = lavaanModel,
    regularized = paste0("l", 6:15),
    lambdas = seq(0,1,length.out = 3),
    returnSubsetParameters = TRUE,
    control = controlGlmnet(nCores = 1))
  set.seed(123)
  regsem2 <- cvLasso(
    lavaanModel = lavaanModel,
    regularized = paste0("l", 6:15),
    lambdas = seq(0,1,length.out = 3),
    returnSubsetParameters = TRUE,
    control = controlGlmnet(nCores = 2))
  
  testthat::expect_equal(all(abs(regsem1@subsetParameters - regsem2@subsetParameters) < 1e-8), TRUE)
  
  # Note: The samples are random. Therefore we reset the seed to get the same samples
  # for single and multi-core execution.
  set.seed(123)
  regsem1 <- cvScad(
    lavaanModel = lavaanModel,
    regularized = paste0("l", 6:15),
    lambdas = seq(0,.3,length.out = 3),
    thetas = c(2.1,2.7,3.5),
    returnSubsetParameters = TRUE,
    control = controlIsta(nCores = 1))
  set.seed(123)
  regsem2 <- cvScad(
    lavaanModel = lavaanModel,
    regularized = paste0("l", 6:15),
    lambdas = seq(0,.3,length.out = 3),
    thetas = c(2.1,2.7,3.5),
    returnSubsetParameters = TRUE,
    control = controlIsta(nCores = 2))
  
  testthat::expect_equal(all(abs(regsem1@subsetParameters - regsem2@subsetParameters) < 1e-8), TRUE)
  
  #### multi-group ####
  set.seed(123)
  dataset1 <- simulateExampleData()
  dataset2 <- simulateExampleData()
  
  pt1 <- pt2 <- parTable(lavaanModel)
  pt1$label[pt1$label == ""] <- paste0(pt1$lhs[pt1$label == ""],
                                       pt1$op[pt1$label == ""],
                                       pt1$rhs[pt1$label == ""]
  )
  pt1$label <- paste0(pt1$label, "_g1")
  pt2$label[pt2$label == ""] <- paste0(pt2$lhs[pt2$label == ""],
                                       pt2$op[pt2$label == ""],
                                       pt2$rhs[pt2$label == ""]
  )
  pt2$label <- paste0(pt2$label, "_g2")
  
  lavaanModel1 <- lavaan::sem(pt1,
                              data = dataset1,
                              meanstructure = TRUE,
                              std.lv = TRUE)
  lavaanModel2 <- lavaan::sem(pt2,
                              data = dataset2,
                              meanstructure = TRUE,
                              std.lv = TRUE)
  lavaanModel <- c(lavaanModel1, lavaanModel2)
  
  # lasso
  # Note: Glmnet uses random update orders for the parameters -> we reset the
  # seed here to make sure that the orders are the same for single and multi-core
  # execution.
  regularized = c(paste0("l", 6:15, "_g1"),
                  paste0("l", 6:15, "_g2"))
  set.seed(123)
  regsem1 <- lasso(
    lavaanModel = lavaanModel,
    regularized = regularized,
    nLambdas = 5,
    control = controlGlmnet(nCores = 1))
  set.seed(123)
  regsem2 <- lasso(
    lavaanModel = lavaanModel,
    regularized = regularized,
    nLambdas = 5,
    control = controlGlmnet(nCores = 2))
  
  testthat::expect_equal(all(abs(regsem1@parameters - regsem2@parameters) < 1e-8), TRUE)
  
  
  # adaptive lasso
  # Note: Glmnet uses random update orders for the parameters -> we reset the
  # seed here to make sure that the orders are the same for single and multi-core
  # execution.
  set.seed(123)
  regsem1 <- adaptiveLasso(
    lavaanModel = lavaanModel,
    regularized = regularized,
    nLambdas = 5,
    control = controlGlmnet(nCores = 1))
  set.seed(123)
  regsem2 <- adaptiveLasso(
    lavaanModel = lavaanModel,
    regularized = regularized,
    nLambdas = 5,
    control = controlGlmnet(nCores = 2))
  
  testthat::expect_equal(all(abs(regsem1@parameters - regsem2@parameters) < 1e-8), TRUE)
  
  # scad
  # Note: The update order is not random -> no need to reset the seed
  regsem1 <- scad(
    lavaanModel = lavaanModel,
    regularized = regularized,
    lambdas = seq(0,.3,length.out = 3),
    thetas = c(2.1,2.7,3.5),
    control = controlIsta(nCores = 1))
  regsem2 <- scad(
    lavaanModel = lavaanModel,
    regularized = regularized,
    lambdas = seq(0,.3,length.out = 3),
    thetas = c(2.1,2.7,3.5),
    control = controlIsta(nCores = 2))
  
  testthat::expect_equal(all(abs(regsem1@parameters - regsem2@parameters) < 1e-8), TRUE)
  
  # mcp
  # Note: The update order is not random -> no need to reset the seed
  regsem1 <- mcp(
    lavaanModel = lavaanModel,
    regularized = regularized,
    lambdas = seq(0,.3,length.out = 3),
    thetas = c(2.1,2.7,3.5),
    control = controlIsta(nCores = 1))
  regsem2 <- mcp(
    lavaanModel = lavaanModel,
    regularized = regularized,
    lambdas = seq(0,.3,length.out = 3),
    thetas = c(2.1,2.7,3.5),
    control = controlIsta(nCores = 2))
  
  testthat::expect_equal(all(abs(regsem1@parameters - regsem2@parameters) < 1e-8), TRUE)
  
  #### transformations ####
  regularized = c(paste0("l", 6:15, "_g1"),
                  paste0("l", 6:15, "_g2"))
  
  transformation <- paste0("parameters: ", paste0(c(paste0("l", 6:15, "_g1"),
                                                    paste0("l", 6:15, "_g2"),
                                                    paste0("delta_l", 6:15)), collapse = ","),
                           "\n",
                           paste0(paste0(paste0("l", 6:15, "_g2"), " = ", paste0("l", 6:15, "_g1") , " + ", paste0("delta_l", 6:15), ";"),
                                  collapse = "\n")
  )
  
  # lasso
  # Note: Glmnet uses random update orders for the parameters -> we reset the
  # seed here to make sure that the orders are the same for single and multi-core
  # execution.
  regularized = paste0("delta_l", 6:15)
  
  set.seed(123)
  regsem1 <- lasso(
    lavaanModel = lavaanModel,
    regularized = regularized,
    nLambdas = 5,
    control = controlGlmnet(nCores = 1),
    modifyModel = modifyModel(transformations = transformation))
  set.seed(123)
  regsem2 <- lasso(
    lavaanModel = lavaanModel,
    regularized = regularized,
    nLambdas = 5,
    control = controlGlmnet(nCores = 2),
    modifyModel = modifyModel(transformations = transformation))
  
  testthat::expect_equal(all(abs(regsem1@parameters - regsem2@parameters) < 1e-8), TRUE)
  
  
  # adaptive lasso
  # Note: Glmnet uses random update orders for the parameters -> we reset the
  # seed here to make sure that the orders are the same for single and multi-core
  # execution.
  set.seed(123)
  regsem1 <- adaptiveLasso(
    lavaanModel = lavaanModel,
    regularized = regularized,
    nLambdas = 5,
    control = controlGlmnet(nCores = 1),
    modifyModel = modifyModel(transformations = transformation))
  set.seed(123)
  regsem2 <- adaptiveLasso(
    lavaanModel = lavaanModel,
    regularized = regularized,
    nLambdas = 5,
    control = controlGlmnet(nCores = 2),
    modifyModel = modifyModel(transformations = transformation))
  
  testthat::expect_equal(all(abs(regsem1@parameters - regsem2@parameters) < 1e-8), TRUE)
  
  # scad
  # Note: The update order is not random -> no need to reset the seed
  regsem1 <- scad(
    lavaanModel = lavaanModel,
    regularized = regularized,
    lambdas = seq(0,.3,length.out = 3),
    thetas = c(2.1,2.7,3.5),
    control = controlIsta(nCores = 1),
    modifyModel = modifyModel(transformations = transformation))
  regsem2 <- scad(
    lavaanModel = lavaanModel,
    regularized = regularized,
    lambdas = seq(0,.3,length.out = 3),
    thetas = c(2.1,2.7,3.5),
    control = controlIsta(nCores = 2),
    modifyModel = modifyModel(transformations = transformation))
  
  testthat::expect_equal(all(abs(regsem1@parameters - regsem2@parameters) < 1e-8), TRUE)
  
  # mcp
  # Note: The update order is not random -> no need to reset the seed
  regsem1 <- mcp(
    lavaanModel = lavaanModel,
    regularized = regularized,
    lambdas = seq(0,.3,length.out = 3),
    thetas = c(2.1,2.7,3.5),
    control = controlIsta(nCores = 1),
    modifyModel = modifyModel(transformations = transformation))
  regsem2 <- mcp(
    lavaanModel = lavaanModel,
    regularized = regularized,
    lambdas = seq(0,.3,length.out = 3),
    thetas = c(2.1,2.7,3.5),
    control = controlIsta(nCores = 2),
    modifyModel = modifyModel(transformations = transformation))
  
  testthat::expect_equal(all(abs(regsem1@parameters - regsem2@parameters) < 1e-8), TRUE)
  
})
