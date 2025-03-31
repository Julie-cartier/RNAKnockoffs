# data simulation

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

# Performances metrics 

fdp <- function(selected, beta_coeff) + return(sum(beta_coeff[selected] == 0) / max(1, length(selected)))

power <- function(selected, beta_coeff) + return(length(which(beta_coeff[selected]!=0))/max(1, length(which(beta_coeff != 0))))


# Implementation of LSCIP method

create.LSCIP_parallel <- function(X) {
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

########################################################################################################################################
#These functions return table with, for each method/statistics and for each target FDR levels, the power, the FDP and the set of selected features.

# Methods comparison 

Knockoff_comp_tot <- function(X, X_k = NULL, y = NULL, beta = NULL, scale = NULL, k = NULL, amplitude = NULL, list_target_fdr, method = NULL){# k: number of non-null features, methods: "linear" or "interaction" 
  
  # X: data matrix of size nxp
  # X_k: X KO matrix
  # y: outcomes (vector of size n)
  # beta: vector of size p (coefficients)
  # scale: If TRUE data are scaled (Default to FALSE)
  # k: number of non null coefficients (non-null features). To simulate interacting settings, k is assumed to be even.
  # amplitude: amplitude of the non null coefficients 
  # list_target_fdr : vector of target FDR for the selection with the KO framework
  # method: either "linear" or "interaction"
  
  n <- nrow(X)
  p <- ncol(X)
  y.copy <- y
  
  mat.features <- data.frame(matrix(0, ncol = 4*length(list_target_fdr), nrow = p+2 ))
  
  if (is.null(X_k)){
    X_k <- create.LSCIP_parallel(X)
  }
  
  if (is.null(y)){
    beta <- beta.sample(p, k, amplitude)
    y <- y.sample(X, beta, method, scale, threshold = 0.5)
  }
  
  # get the knockoff matrices for  different methods
  
  X_LSCIP <- create.LSCIP_parallel(X)
  
  knock_sampler <- knockpy$knockoffs$GaussianSampler(X, method = 'mvr')
  X_MVR <- knock_sampler$sample_knockoffs()
  X_MVR <- py_to_r(X_MVR)
  
  knock_sampler <- knockpy$knockoffs$GaussianSampler(X, method = 'ci')
  X_CI <- knock_sampler$sample_knockoffs()
  X_CI <- py_to_r(X_CI)
  
  knock_sampler <- knockpy$knockoffs$GaussianSampler(X, method = 'sdp')
  X_SDP <- knock_sampler$sample_knockoffs()
  X_SDP <- py_to_r(X_SDP)
  
  
  # Compute the LCD statistic
  
  W_MVR <- stat.glmnet_coefdiff(X, X_MVR, y, nfolds=10, family="binomial")
  W_SDP <- stat.glmnet_coefdiff(X, X_SDP, y, nfolds=10, family="binomial")
  W_CI <- stat.glmnet_coefdiff(X, X_CI, y, nfolds=10, family="binomial")
  W_LSCIP <- stat.glmnet_coefdiff(X, X_LSCIP, y, nfolds=10, family="binomial")
  
  
  for (i in 1:length(list_target_fdr)){
    
    target_fdr <- list_target_fdr[i]
    
    thres_SDP <- knockoff.threshold(W_SDP, fdr=target_fdr, offset=1) 
    selected_SDP <- which(W_SDP >= thres_SDP)
    
    thres_MVR <- knockoff.threshold(W_MVR, fdr=target_fdr, offset=1) 
    selected_MVR <- which(W_MVR >= thres_MVR)
    
    thres_LSCIP <- knockoff.threshold(W_LSCIP, fdr=target_fdr, offset=1) 
    selected_LSCIP <- which(W_LSCIP >= thres_LSCIP)
    
    
    thres_CI <- knockoff.threshold(W_CI, fdr=target_fdr, offset=1) 
    selected_CI <- which(W_CI >= thres_CI)
    
    
    mat.features[1, (1+4*(i-1)):(4+4*(i-1))] <- c(fdp(selected_SDP, beta), fdp(selected_MVR, beta), fdp(selected_LSCIP, beta), fdp(selected_CI, beta))
    mat.features[2, (1+4*(i-1)):(4+4*(i-1))] <- c(power(selected_SDP, beta), power(selected_MVR, beta), power(selected_LSCIP, beta), power(selected_CI, beta))
    mat.features[selected_SDP + 2, 1+4*(i-1)] <- 1
    mat.features[selected_MVR + 2, 2+4*(i-1)] <- 1
    mat.features[selected_LSCIP + 2, 3+4*(i-1)] <- 1
    mat.features[selected_CI + 2, 4+4*(i-1)] <- 1
    
    
    colnames(mat.features)[c(1+4*(i-1), 2+4*(i-1), 3+4*(i-1), 4+4*(i-1)) ] <- c(str_glue("SDP", target_fdr), str_glue("MVR", target_fdr), str_glue("LSCIP", target_fdr), str_glue("CI", target_fdr))
    
  }
  
  rownames(mat.features) <- c("FDP", "power", c(1:p))
  
  if (is.null(y.copy)){# get the beta vector in the case where we did not provide it
    return(list(mat.features, beta))
  }
  else{
    return(mat.features)
  }
}

# Statistics comparison
Knockoff_comp_stat_LSCIP <- function(X, X_k = NULL, y = NULL, beta = NULL, scale = NULL, k = NULL, amplitude = NULL, list_target_fdr, list_lambda, method = NULL){# list lambda is usef for comparison with the LASSO
  
  # X: data matrix of size nxp
  # X_k: X KO matrix
  # y: outcomes (vector of size n)
  # beta: vector of size p (coefficients)
  # scale: If TRUE data are scaled (Default to FALSE)
  # k: number of non null coefficients (non-null features). To simulate interacting settings, k is assumed to be even.
  # amplitude: amplitude of the non null coefficients 
  # list_target_fdr : vector of target FDR for the selection with the KO framework
  # list_lambda : vector of penalization values (lambda) for the selection with LASSO
  # method: either "linear" or "interaction"
  
  n <- nrow(X)
  p <- ncol(X)
  
  y.copy = y
  
  mat.features <- data.frame(matrix(0, ncol = 5*length(list_target_fdr), nrow = p+2 ))
  
  # create LSCIP KO
  
  if (is.null(X_k)){
    X_k <- create.LSCIP_parallel(X)
  }
  
  # simulate data 
  
  if (is.null(y)){
    beta <- beta.sample(p, k, amplitude)
    y <- y.sample(X, beta, method, scale, threshold = 0.5)
  }
  
  # Compute statistics
  
  # EN 
  
  model_lasso <- cv.glmnet(as.matrix(data.frame(X,X_k)), y, alpha = 0.5, type.measure = "deviance", family = "binomial", nfolds = 10)
  coef_lasso <- coef(model_lasso, s = "lambda.min")
  coef <- coef_lasso[-1]
  orig <- 1:p
  W_EN <- abs(coef[orig])-abs(coef[orig+p])
  
  # Lasso
  
  W_LCD <- stat.glmnet_coefdiff(X, X_k, y, nfolds=10, family="binomial")
  
  # MLR
  
  X_py <- numpy$ascontiguousarray(X)
  X_k_py <- numpy$ascontiguousarray(X_k)
  y_py <- numpy$ascontiguousarray(numpy$ravel(y))
  
  KO_MLR <- knockpy$mlr$MLR_Spikeslab()
  KO_MLR$fit(X_py, X_k_py, y_py, groups = numpy$arange(1, p + 1, 1))
  W_MLR <- KO_MLR$W
  W_MLR <- py_to_r(W_MLR)
  
  # RF
  
  y.factor <- as.factor(y)
  W_RF <- stat.random_forest(X, X_k, y.factor)
  
  # select features for different target fdr levels
  
  for (i in 1:length(list_target_fdr)){
    target_fdr <- list_target_fdr[i]
    lambda <- list_lambda[i]
    
    thres_EN <- knockoff.threshold(W_EN, fdr=target_fdr, offset=1) 
    selected_EN <- which(W_EN>= thres_EN)
    
    
    thres_LCD <- knockoff.threshold(W_LCD, fdr=target_fdr, offset=1) 
    selected_LCD <- which(W_LCD >= thres_LCD)
    
    thres_MLR <- knockoff.threshold(W_MLR, fdr=target_fdr, offset=1) 
    selected_MLR <- which(W_MLR >= thres_MLR)
    
    thres_RF <- knockoff.threshold(W_RF, fdr=target_fdr, offset=1) 
    selected_RF <- which(W_RF >= thres_RF)
    
    # get the fdp and the power obtained with the Lasso penalized logistic regression for different penalization parameters
    
    fitlog <- glmnet(X, y, type.measure = "deviance", family = "binomial", lambda = sort(list_lambda, decreasing = TRUE))
    coef_lasso <- coef(fitlog, s = lambda)
    selected_lasso <- which(coef_lasso[-1]!=0)
    
    predicted_y  <- predict(fitlog, s = lambda, X, type='response')
    predicted_y[predicted_y>0.5]  <- 1
    predicted_y[predicted_y<=0.5] <- 0
    
    mat.features[1, (1+5*(i-1)):(5+5*(i-1))] <- c(fdp(selected_LCD, beta), fdp(selected_EN, beta), fdp(selected_MLR, beta), fdp(selected_RF, beta), fdp(selected_lasso, beta))
    mat.features[2, (1+5*(i-1)):(5+5*(i-1))] <- c(power(selected_LCD, beta), power(selected_EN, beta), power(selected_MLR, beta), power(selected_RF, beta), power(selected_lasso, beta))
    mat.features[selected_LCD + 2, 1+5*(i-1)] = 1
    mat.features[selected_EN + 2, 2+5*(i-1)] = 1
    mat.features[selected_MLR + 2, 3+5*(i-1)] = 1
    mat.features[selected_RF + 2, 4+5*(i-1)] = 1
    mat.features[selected_lasso + 2, 5+5*(i-1)] = 1
    
    colnames(mat.features)[c(1+5*(i-1), 2+5*(i-1), 3+5*(i-1), 4+5*(i-1), 5+5*(i-1)) ] <- c(str_glue("LCD", target_fdr), str_glue("Elastic_net", target_fdr), str_glue("MLR", target_fdr), str_glue("RF", target_fdr), str_glue("lasso_", lambda))
    
  }
  
  rownames(mat.features) <- c("FDP", "power", c(1:p))
  
  if (is.null(y.copy)){# get the beta vector in the case where we did not provide it
    return(list(mat.features, beta))
  }
  else{
    return(mat.features)
  }
  
}

