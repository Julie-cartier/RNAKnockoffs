#library(knockoff)
library(glmnet)
library(doParallel)
library(stringr)

# This script has been use to create a bank of KO matrices: (100 to test stability) x (Draws for aggregation methods = 100)-> 10000 KO matrices. Given the number of required KO matrices, we compute this bank using computing cluster
# This script does not depend on other scripts
### Settings

RNA.data <- read.table(file = "/home/julie/Documents/Paper_codes/Data/CRUKPAP/dat_vst_NS_no_Inel_selected.txt", head = T) 
X_real <- as.matrix(RNA.data)
n = 369 # Number of samples
p = 749 # Number of variables

create.LSCIP_parallel <- function(X){
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

list_KO_real_LSCIP <- list()



for (i in (1:10000)){
  print(i)
  X_k <- create.LSCIP_parallel(scale(X_real))
  list_KO_real_LSCIP[[i]] <- X_k
}


save(list_KO_real_LSCIP, file = "/cluster/CBIO/data1/jcartier/Knockoffs/KO_Matrices_Bank/KO_real_LSCIP_bank.R")


load("/cluster/CBIO/data1/jcartier/Knockoffs/KO_Matrices_Bank/KO_real_LSCIP_bank.R")

# Separate the library of 10000 matrices into smaller matrices (easier to load in an R script)

for (i in 1:100){
  print(i)
  list_X_k_100 <- list_KO_real_LSCIP[(100*(i-1)+1):(100*i)]
  save(list_X_k_100, file = str_glue("/media/julie/T5 EVO/Data_cluster/KO_bank_under_sets/X_k_",i,"_100.R"))
}
