## ----install_BNrich, eval=FALSE-----------------------------------------------
#  install.packages("BNrich_0.1.0.tar.gz", type="source", repos=NULL)
#  library("BNrich")

## ----start files, eval=FALSE--------------------------------------------------
#  destfile = tempfile("files", fileext = ".rda")
#  files <- fetch_data_file()
#  load(destfile)

## ----example_of_dataset, eval=FALSE-------------------------------------------
#  Data <- system.file("extdata", "Test_DATA.RData", package = "BNrich", mustWork = TRUE)
#  load(Data)
#  head(dataH)

## ----example of dataset, eval=FALSE-------------------------------------------
#  head(dataD)

## ----unify_path, eval=FALSE---------------------------------------------------
#  unify_results <- unify_path(dataH, dataD, mapkG, pathway.id)

## ----the_number_of_original_pathways, eval=FALSE------------------------------
#  length(mapkG)

## ----the_number_of_unified_pathways, eval=FALSE-------------------------------
#  mapkG1 <- unify_results$mapkG1
#  length(mapkG1)

## ----the_first_original_pathwayID, eval=FALSE---------------------------------
#  pathway.id[1]

## ----the_first_original_pathway_information, eval=FALSE-----------------------
#  mapkG[[1]]

## ----the_first_unified_pathwayID, eval=FALSE----------------------------------
#  pathway.id1 <- unify_results$pathway.id1
#  pathway.id1[1]

## ----the_first_unified_pathway_information, eval=FALSE------------------------
#  mapkG1[[1]]

## ----construct_BN, eval=FALSE-------------------------------------------------
#  BN <- BN_struct(unify_results$mapkG1)

## ----LASSO_function, eval=FALSE-----------------------------------------------
#  data_h <- unify_results$data_h
#  data_d <- unify_results$data_d
#  LASSO_results <- LASSO_BN(BN, data_h, data_d)

## ----the_number_of_edge_in_first_initial_BN, eval=FALSE-----------------------
#  nrow(arcs(BN[[1]]))

## ----the_number_of_edge_in_first_control_BN, eval=FALSE-----------------------
#  nrow(arcs(LASSO_results$BN_H[[1]]))

## ----the_number_of_edge_in_first_disease_BN, eval=FALSE-----------------------
#  nrow(arcs(LASSO_results$BN_D[[1]]))

## ----estimate_parameters, eval=FALSE------------------------------------------
#  BN_H <- LASSO_results$BN_H
#  BN_D <- LASSO_results$BN_D
#  esti_results <- esti_par(BN_H, BN_D, data_h, data_d)

## ----parameters_information_in_BN_h, eval=FALSE-------------------------------
#  esti_results$BNs_h[[1]]$` hsa:1978`

## ----parameters_information_in_BN_d, eval=FALSE-------------------------------
#  esti_results$BNs_d[[1]]$`hsa:1978`

## ----variance_matrix_function, eval=FALSE-------------------------------------
#  BN_h <- esti_results$BNs_h
#  BN_d <- esti_results$BNs_d
#  coef_h <- esti_results$coef_h
#  coef_d <- esti_results$coef_d
#  var_mat_results<- var_mat (data_h, coef_h, BN_h, data_d, coef_d, BN_d)

## ----variance-covariance_matrixe_for_the_fifth_node_in_first_BNh, eval=FALSE----
#  (var_mat_results$var_mat_Bh[[1]])[5]

## ----variance-covariance_matrixe_for_the_fifth_node_in_first_BNd, eval=FALSE----
#  (var_mat_results$var_mat_Bd[[1]])[5]

## ----parm_Ttest_function, eval=FALSE------------------------------------------
#  var_mat_Bh <- var_mat_results $var_mat_Bh
#  var_mat_Bd <- var_mat_results $var_mat_Bd
#  Ttest_results <- parm_Ttest(data_h, coef_h, BN_h, data_d, coef_d, BN_d, var_mat_Bh, var_mat_Bd, pathway.id1)
#  head(Ttest_results)

## ----BNrich_function, eval=FALSE----------------------------------------------
#  BNrich_results <- BNrich(Ttest_results, pathway.id1, PathName_final, fdr.value = 0.05)
#  head(BNrich_results)

