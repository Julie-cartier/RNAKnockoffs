library(doParallel)
library(glmnet)

# This script has been use to create a bank of 100 KO matrices for the breast cancer data matrix :. Given the dimension of the data matrix, we compute this bank using computing cluster


load("/home/julie/Documents/Paper_codes/Selection_with_true_y/R_files/BC_data.R")


create.LSCIP_parallel <- function(X) {
  p = ncol(X)
  n = nrow(X)
  mat_residuals <- matrix(nrow = n, ncol = p)
  mat_lasso_x <- matrix(nrow = n, ncol = p)
  
  # Create a foreach loop for parallelization
  cl <- parallel::makeCluster(12)
  doParallel::registerDoParallel(cl)
  
  result <- foreach(j = 1:p, .combine = "cbind",.packages = c("glmnet")) %dopar% {
    y = X[, j]
    X_ = X[, -j]
    modelj <- cv.glmnet(X_, y)
    betaj = as.matrix(coef(modelj, s = "lambda.min"))
    X_j <- cbind(rep(1, n), X_)
    residuals_j = y - X_j %*% betaj
    lasso_x_j = X_j %*% betaj
    
    # Return the results for this iteration
    list(residuals_j, lasso_x_j)
  }
  stopCluster(cl)
  # Combine the results into matrices
  for (i in 1:p) {
    mat_residuals[, i] <- result[[2*i-1]]
    mat_lasso_x[, i] <- result[[2*i]]
  }
  
  
  X_KO <- mat_lasso_x + mat_residuals[, sample(p)]
  return(X_KO)
}

D = 100

list.KO.BC <- list()

for (i in 1:D){
  print(i)
  list.KO.BC[[i]] <- create.LSCIP_parallel(X_BC)
}

save(list.KO.BC, file = "/home/julie/Documents/Paper_codes/Data/BC/KO_BC_bank.R")



