beta.sample <- function(p, k, amplitude){ 
  
  # p: vector size (number of features)
  # k: number of non null coefficients (non-null features). To simulate interacting settings, k is assumed to be even.
  # amplitude: amplitude of the non null coefficients 
  
  nonzero <- sample(p, k) # randomly choose the indices of the k non-null coefficients 
  beta <- rep(0,p)
  beta[nonzero] <- amplitude * sample(c(1,-1), k, replace = T) # Randomly multiply amplitude by +1 or -1
  return(beta) 
}

y.sample <- function(X, beta, method, scale = FALSE, threshold = 0.5){ 
  
  # X: data matrix of size nxp
  # beta: vector of size p (coefficients)
  # method: either "linear" or "interaction" 
  # scale: If TRUE data are scaled (Default to FALSE)
  # Threshold : used to assign classes 
  
  
  n <- nrow(X)
  p <- ncol(X)
  
  if (scale == TRUE){X = scale(X)}
  
  if (method == "linear"){ # y is simulated with a linear logistic model
    y_prob <- 1 / (1 + exp(- (X%*%beta))) 
    y <- 1*(y_prob>=threshold) # assign class
  }
  else if (method == "interaction"){ # y is simulated with a logistic model and with pairwise interactions between features
    
    
    beta.ind <- which(beta!=0) # get the indices of non-null features
    X.part <- X[, sample(beta.ind)]
    k <- length(beta.ind) # number of non null features
    amplitude <- abs(beta[beta.ind][1])# get amplitude
    
    choice <- sample(c(1, -1), 1)
    
    # Define the amplitude 
    
    # select amplitude sign
    
    if (choice == 1) { # if k/2 is odd
      vector <- c(rep(1, round(k/4)), rep(-1, k/2-round(k/4)))
    } else {
      vector <- c(rep(1, k/2-round(k/4)), rep(-1, round(k/4)))
    }
    
    shuffled.vector <- sample(vector)
    
    # get non-null coefficients
    
    coeff <- amplitude * shuffled.vector
    
    int.product = 0
    
    for (j in 1:(k/2)){ # pairwise interaction
      int.product <- int.product + coeff[j]*X.part[,2*j-1]*X.part[,2*j]
    }
    y_prob <- 1 / (1 + exp(- int.product ))
    y <- 1*(y_prob>=threshold)
  }
  if (sum(y)<=(n*0.45) | sum(y)>(n*0.55)){
    warning(str_glue("unbalanced classes : proportion of 1 : ",sum(y)/n))
  }
  return(y)
}

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

create_KO_stats <- function(X, list_X_k, y, num_cores = 10) {
  
  # X: data matrix
  # list_X_k: list of KO matrices
  # y: outcomes (vector of size n)
  
  Draws <- length(list_X_k)
  
  # Use mclapply for parallel computation
  KO_stats_list <- mclapply(list_X_k, function(X_k) {
    stat.glmnet_coefdiff(X, X_k, y, nfolds = 10, family = "binomial")
  }, mc.cores = num_cores)
  
  # Combine the list of results into a matrix
  KO_stats_mat <- do.call(rbind, KO_stats_list)
  
  return(KO_stats_mat)
}



# These functions return tables with, for each method and for each target FDP/target FDR/lambda parameter the set of selected features. In addition, for the LASSO, tables contain AUCs. 

KOPI_KO_stability <- function(X, y, list_X_k){
  
  # X: data matrix
  # y: outcomes (vector of size n)
  # list_X_k: list of KO matrices


  n <- nrow(X)
  p <- ncol(X)
  
  table.coeff.Vanilla <- data.frame(matrix(0, ncol = 4  , nrow = p ))
  table.coeff.KOPI <- data.frame(matrix(0, ncol = 4, nrow = p ))
  
  ## Vanilla KO
  
  X_k <- list_X_k[[1]]
  
  W_LCD <- stat.glmnet_coefdiff(X, X_k, y, nfolds=10, family="binomial")
  
  
  thres_LCD_0.1 <- knockoff.threshold(W_LCD, fdr=0.1, offset=1)
  selected_LCD_0.1 <- which(W_LCD >= thres_LCD_0.1)
  
  thres_LCD_0.2 <- knockoff.threshold(W_LCD, fdr=0.2, offset=1)
  selected_LCD_0.2 <- which(W_LCD >= thres_LCD_0.2)
  
  thres_LCD_0.3 <- knockoff.threshold(W_LCD, fdr=0.3, offset=1)
  selected_LCD_0.3 <- which(W_LCD >= thres_LCD_0.3)
  
  thres_LCD_0.5 <- knockoff.threshold(W_LCD, fdr=0.5, offset=1)
  selected_LCD_0.5 <- which(W_LCD >= thres_LCD_0.5)
  
  
  table.coeff.Vanilla[selected_LCD_0.1, 1] <- 1
  table.coeff.Vanilla[selected_LCD_0.2, 2] <- 1
  table.coeff.Vanilla[selected_LCD_0.3, 3] <- 1
  table.coeff.Vanilla[selected_LCD_0.5, 4] <- 1
  
  
  ## KOPI
  
  Draws <- length(list_X_k)
  
  KO_stats <- create_KO_stats(X, list_X_k, y, num_cores = 10)
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
  
  table.coeff.KOPI[selected_KOPI_0.1, 1] <- 1
  table.coeff.KOPI[selected_KOPI_0.2, 2] <- 1
  table.coeff.KOPI[selected_KOPI_0.3, 3] <- 1
  table.coeff.KOPI[selected_KOPI_0.5, 4] <- 1
  
  return(list(table.coeff.Vanilla, table.coeff.KOPI))
  
}

LPLR_stability <- function(X_init, y_init, vec.ind){
  
  # X_init: data matrix (before subsampling)
  # y_init: outcomes (vector of size n)
  # vec.ind: vector of indices. (Indices of samples removed from the data matrix) 
  
  X <- X_init[-vec.ind, ]
  X_test <- X_init[vec.ind, ]
  
  y <- y_init[-vec.ind]
  y_test <- y_init[vec.ind]
  
  n <- nrow(X)
  p <- ncol(X)
  
  table.coeff.LPLR <- data.frame(matrix(0, ncol = 1, nrow = 2+p))
  
  ## LPLR
  
  cv_fitlog <- cv.glmnet(X, y, type.measure = "deviance", family = "binomial", nfolds = 10)
  coef_lasso <- coef(cv_fitlog, s = "lambda.min")
  table.coeff.LPLR[which(coef_lasso[-1]!=0)+2, 1] <- 1
  
  # AUC
  
  mat_i <- matrix(data = NA, nrow = nrow(X_test), ncol = length(coef_lasso))
  mat_i[,1] <- rep(1, nrow(X_test))
  mat_i[,2:length(coef_lasso)] <- X_test
  
  vec.prod <- mat_i %*% coef_lasso
  score <- 1/(exp(-vec.prod)+1)
  
  
  ROC <- roc.curve(scores.class0 = score, weights.class0 = y_test, curve = FALSE)
  PR <- pr.curve(scores.class0 = score, weights.class0 = y_test, curve = FALSE)
  table.coeff.LPLR[1, ] <- ROC$auc
  table.coeff.LPLR[2, ] <- PR$auc.integral
  
  return(table.coeff.LPLR)
  
}

LPLR_stability_oracle <- function(X_init, y_init, vec.ind, beta){
  
  # X_init: data matrix (before subsampling)
  # y_init: outcomes (vector of size n)
  # vec.ind: vector of indices. (Indices of samples removed from the data matrix) 
  # beta : vector used to compute y 
  
  X <- X_init[-vec.ind, ] # subsampling
  X_test <- X_init[vec.ind, ]
  
  y <- y_init[-vec.ind] # subsampling
  y_test <- y_init[vec.ind]
  
  n <- nrow(X)
  p <- ncol(X)
  
  table.coeff.LPLR <- data.frame(matrix(0, ncol = 1, nrow = 4 + p))
  
  ## get lambda oracle
  
  vec.lambda <-  seq(1, 0, -0.001) # sequence of lambda used in the LASSO
  vec.ratio.TD.FD <- c()
  
  # fit the logistic regression
  
  fitlog <- glmnet(X, y, type.measure = "deviance", family = "binomial", lambda = vec.lambda, nfolds = 10)
  
  for (lambda.value in vec.lambda){ # get the ratio TD/FD for all lambda
    
    coef_lasso <- coef(fitlog, s = lambda.value)
    selected_lasso <- which(coef_lasso[-1]!=0)
    vec.ratio.TD.FD <- c(vec.ratio.TD.FD, sum(beta[selected_lasso] != 0)/sum(beta[selected_lasso] == 0))
    
  }
  
  lambda.oracle <- vec.lambda[which.max(vec.ratio.TD.FD)] # get the lambda that maximizes the TD/FD  ratio
  coef_lasso <- coef(fitlog, s = lambda.oracle)
  table.coeff.LPLR[which(coef_lasso[-1]!=0) + 4, 1] <- 1
  selected_lasso <- which(coef_lasso[-1]!=0) # get the selected features (faetures with a non null coefficient)
  
  # Compute AUCs (Safety check)
  
  mat_i <- matrix(data = NA, nrow = nrow(X_test), ncol = length(coef_lasso))
  mat_i[,1] <- rep(1, nrow(X_test))
  mat_i[,2:length(coef_lasso)] <- X_test
  
  vec.prod <- mat_i %*% coef_lasso
  score <- 1/(exp(-vec.prod)+1)
  
  
  ROC <- roc.curve(scores.class0 = score, weights.class0 = y_test, curve = FALSE)
  PR <- pr.curve(scores.class0 = score, weights.class0 = y_test, curve = FALSE)
  table.coeff.LPLR[1, ] <- ROC$auc
  table.coeff.LPLR[2, ] <- PR$auc.integral
  table.coeff.LPLR[3, ] <- lambda.oracle
  table.coeff.LPLR[4, ] <- max(vec.ratio.TD.FD, na.rm = TRUE) 
  return(table.coeff.LPLR)
  
}
