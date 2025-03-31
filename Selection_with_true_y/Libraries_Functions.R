# Libraries

library(doParallel)
library(tictoc)

library(knockoff)
library(glmnet)

library(stringr)
library(tidyverse)
library(ggplot2)

library(reticulate)


# LSCIP method

create.LSCIP_parallel <- function(X){
  
  # X: data matrix
  
  p <- ncol(X)
  n <- nrow(X)
  mat_residuals <- matrix(nrow = n, ncol = p)
  mat_lasso_x <- matrix(nrow = n, ncol = p)
  
  cl <- parallel::makeCluster(12)
  doParallel::registerDoParallel(cl)
  
  result <- foreach(j = 1:p, .combine = "cbind",.packages = c("glmnet")) %dopar% {
    y = X[, j]
    X_ = X[, -j]
    modelj <- cv.glmnet(X_, y)
    betaj <- as.matrix(coef(modelj, s = "lambda.min"))
    X_j <- cbind(rep(1, n), X_)
    residuals_j = y - X_j %*% betaj
    lasso_x_j = X_j %*% betaj
    
    list(residuals_j, lasso_x_j)
  }
  stopCluster(cl)
  # Combine the results in a matrix
  for (i in 1:p) {
    mat_residuals[, i] <- result[[2*i-1]]
    mat_lasso_x[, i] <- result[[2*i]]
  }
  
  
  X_KO <- mat_lasso_x + mat_residuals[, sample(p)] # Randomly assign residuals
  return(X_KO)
}



create_KO_stats <- function(X, list_X_k, y, num_cores = 5){
  
  # X: data matrix (nxp)
  # y: otcomes (vector of size n)
  # list_X_k: list of D pre-computed KO matrices
  
  Draws <- length(list_X_k)
  
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  KO_stats_list <- foreach(X_k = list_X_k, .combine = rbind, .packages = "knockoff") %dopar% { # Do parallel function seems more robust to high dimension tha mclaplly
    stat.glmnet_coefdiff(X, X_k, y, nfolds = 10, family = "binomial")
  }
  
  stopCluster(cl)
  
  return(KO_stats_list)
}


# Apply the KOPI scheme to perform selection 
# get table with selected features for different target FDR 

KOPI_selection <- function(X, y, list_X_k, KO_stats = NULL){
  
  # X: data matrix (nxp)
  # y: otcomes (vector of size n)
  # list_X_k: list of D pre-computed KO matrices
  # KO_stats: Computed matrix of statistics (with create_KO_stats)
  
  # cross-validation
  n <- nrow(X)
  p <- ncol(X)
  
  table.coeff.KOPI <- data.frame(matrix(0, ncol = 4  , nrow = p ))
  
  
  # KOPI
  
  Draws <- length(list_X_k)
  
  if (is.null(KO_stats)){
    KO_stats <- create_KO_stats(X, list_X_k, y, num_cores = 5)
  }
  
  
  KO_stats_py <- r_to_py(KO_stats)
  
  # get pi0_hmean and learned_tpl_hmean
  
  result <- create_hmean_p_values_learned_tpl(p = as.integer(p), draws = as.integer(Draws), B = as.integer(1000), n_jobs = as.integer(10))
  
  pi0_hmean <- result[[1]]
  learned_tpl_hmean <- result[[2]]
  
  # select features
  
  set_selected_0.1 <- KOPI_given_all(ko_stats = KO_stats_py, draws = Draws, alpha = 0.1, target_fdp = 0.1, learned_tpl_hmean = learned_tpl_hmean, pi0_hmean = pi0_hmean, k_max = as.integer(round(p/50)))
  selected_KOPI_0.1 <- set_selected_0.1[[1]] + 1
  
  set_selected_0.2 <- KOPI_given_all(ko_stats = KO_stats_py, draws = Draws, alpha = 0.1, target_fdp = 0.2, learned_tpl_hmean = learned_tpl_hmean, pi0_hmean = pi0_hmean, k_max = as.integer(round(p/50)))
  selected_KOPI_0.2 <- set_selected_0.2[[1]] + 1
  
  set_selected_0.3 <- KOPI_given_all(ko_stats = KO_stats_py, draws = Draws, alpha = 0.1, target_fdp = 0.3, learned_tpl_hmean = learned_tpl_hmean, pi0_hmean = pi0_hmean, k_max = as.integer(round(p/50)))
  selected_KOPI_0.3 <- set_selected_0.3[[1]] + 1
  
  set_selected_0.5 <- KOPI_given_all(ko_stats = KO_stats_py, draws = Draws, alpha = 0.1, target_fdp = 0.5, learned_tpl_hmean = learned_tpl_hmean, pi0_hmean = pi0_hmean, k_max = as.integer(round(p/50)))
  selected_KOPI_0.5 <- set_selected_0.5[[1]] + 1
  
  table.coeff.KOPI[selected_KOPI_0.1, 1] = 1
  table.coeff.KOPI[selected_KOPI_0.2, 2] = 1
  table.coeff.KOPI[selected_KOPI_0.3, 3] = 1
  table.coeff.KOPI[selected_KOPI_0.5, 4] = 1
  
  colnames(table.coeff.KOPI) <- c("q = 0.1"," q= 0.2", "q = 0.3", "q = 0.5")
  rownames(table.coeff.KOPI) <- colnames(X)
  
  return(table.coeff.KOPI)
  
}
