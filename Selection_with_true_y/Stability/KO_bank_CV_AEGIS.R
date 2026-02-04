# This script has been use to create a bank of KO matrices: (10 x 10-fold CV). Given the number of required KO matrices, we compute this bank using computing cluster
# Hence this script does not depend on other scripts

#library(knockoff)
library(glmnet)
library(doParallel)
library(stringr)


# Settings
# get X_CRUKPAP

RNA.data <- load("/cluster/CBIO/data1/jcartier/data/AEGIS/AEGIS_data.R")
X_real = X_BC
n = nrow(X_real) # Number of samples
p = ncol(X_real) # Number of variables


# Function to create KO matrices (LSCIP)

create.LSCIP_parallel <- function(X) {
  # original data matrix
  
  p = ncol(X)
  n = nrow(X)
  mat_residuals <- matrix(nrow = n, ncol = p)
  mat_lasso_x <- matrix(nrow = n, ncol = p)
  
  cl <- parallel::makeCluster(10)
  doParallel::registerDoParallel(cl)
  
  result <- foreach(j = 1:p, .combine = "cbind",.packages = c("glmnet")) %dopar% {
    y = X[, j]
    X_ = X[, -j]
    modelj <- cv.glmnet(X_, y)
    betaj = as.matrix(coef(modelj, s = "lambda.min"))
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


list_KO_cv_LSCIP <- list()
list_y_cv_LSCIP <- list()

X_scaled = scale(X_real)

i <- 1

vec <- sample(1:n)
split_vec <- split(vec, cut(seq_along(vec), 10, labels = FALSE)) # split the vector for one cross validation
save(split_vec, file = str_glue("/cluster/CBIO/data1/jcartier/Knockoffs/KO_Matrices_Bank/KO_bank_CV_sets_AEGIS/CV_indices/CV_",i,".RData")) # save the indices (for comparison with other methods)
for (j in 1:10){
  X_cv = X_scaled[-split_vec[[j]], ]
    
  list_KO_cv_100 <- list()
    
  for (k in 1:100){# given a data matrix with n = 0.9*n samples and p features create 100 KO matrices for KO aggregation schemes
      
    print(str_glue("i",i,", j",j," ,k", k))
    X_k = create.LSCIP_parallel(X_cv)
      
    list_KO_cv_100[[k]] <- X_k
    list_KO_cv_LSCIP[[length(list_KO_cv_LSCIP) + 1]] <- X_k
  }
  save(list_KO_cv_100, file = str_glue("/cluster/CBIO/data1/jcartier/Knockoffs/KO_Matrices_Bank/KO_bank_CV_sets_AEGIS/CV_",i,"_fold_",j,"_X_k_100.RData")) # save into different files
  }

