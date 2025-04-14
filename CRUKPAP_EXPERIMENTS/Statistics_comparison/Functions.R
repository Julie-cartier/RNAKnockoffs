# Simulate data

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

# Compute performance metrics

fdp <- function(selected, beta_coeff) + return(sum(beta_coeff[selected] == 0) / max(1, length(selected)))

power <- function(selected, beta_coeff) + return(length(which(beta_coeff[selected]!=0))/max(1, length(which(beta_coeff != 0))))

# Create LSCIP KO

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

# Compute RRB statistic 

# From the original stat.mboost_varimp2 function provided in the code of the paper : Shen A, Fu H, He K, Jiang H (2019) False Discovery Rate Control in Cancer Biomarker Selection using Knockoffs (https://jhui2014.github.io/knockoff/)
W.mboost_varimp2 <- function(X, X_k, y, max.mstop = 200, bl = c("btree"), k_fold = 10, family = Binomial()){ 
  
  # X: data matrix
  # X_k: X corresponding KO matrix
  # y: outcomes (vector of size n)
  # max.stop: parameter of the mboost algorithm ("an integer giving the number of initial boosting iterations") 
  # bl : base learners used in the mboost algorithm (mboost function parameter)
  # k_fold: number of fold for the cross validation
  # family: type of prediction (depend on y) (mboost function parameter) 
  
  Xaug <- data.frame(X, X_k, y)
  
  Xaug$y <- as.factor(Xaug$y)
  
  model <- mboost(y ~ ., data = Xaug, control = boost_control(mstop = max.mstop), baselearner = bl, family = family) # use the first base learner of the list 
  
  cv.kf <- cv(model.weights(model), type="kfold", B = k_fold) # function that separates the data between train data and validation data for cross validation
  
  cvm <- cvrisk(model, folds = cv.kf, papply = mclapply) ### used to choose the number of boosting iterations : hyperparameter chosen via cross validation
  
  best.model <- model[mstop(cvm)]
  Z <- as.numeric(varimp(best.model))
  
  orig <- 1:p
  W <- Z[orig] - Z[orig+p]
  
  return(W)
}

# Comparison functions (to create tables with saved results and performances). THe goal is to compare the performance of the KO methods for different statistics
# These functions returns tables with, for each method and for each target FDR and/or lambda value, the power, the FDP and the set of selected features 

Knockoff_comp_stat_LSCIP <- function(X, X_k = NULL, y = NULL, beta = NULL, scale = NULL, k = NULL, amplitude = NULL, list_target_fdr, list_lambda, method = NULL){
  
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
  
  if (is.null(X_k)){ #create KO matrx
    X_k <- create.LSCIP_parallel(X)
  }
  
  if (is.null(y)){# Simulate y
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
  
  # LCD
  
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
  
  # Select features for different target fdr levels
  
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
    
    # Compute LASSO for different penalization parameter
    
    fitlog <- glmnet(X, y, type.measure = "deviance", family = "binomial", lambda = sort(list_lambda, decreasing = TRUE))
    coef_lasso <- coef(fitlog, s = lambda)
    selected_lasso <- which(coef_lasso[-1]!=0)
    
    # LASSO prediction
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
  
  if (is.null(y.copy)){
    return(list(mat.features, beta))
  }
  else{
    return(mat.features)
  }
  
}

Knockoff_LassoBasedMethod <- function(X, X_k = NULL, y = NULL, beta = NULL, scale = NULL, k = NULL, amplitude = NULL, list_target_fdr, method = NULL){
  
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
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  y.copy = y
  
  mat.features <- data.frame(matrix(0, ncol = 3*length(list_target_fdr), nrow = p+2 ))
  
  if (is.null(X_k)){# Compute KO matrix 
    X_k <- create.LSCIP_parallel(X)
  }
  
  if (is.null(y)){# Simulate y
    beta <- beta.sample(p, k, amplitude)
    y <- y.sample(X, beta, method, scale, threshold = 0.5)
  }
  
  # Compute statistic
  
  #LCD
  
  W_LCD <- stat.glmnet_coefdiff(X, X_k, y, family="binomial")
  
  # LSM
  
  W_LSM <- stat.glmnet_lambdasmax(X, X_k, y, family="binomial")
  
  # LSM
  
  W_LM <- stat.glmnet_lambdadiff(X, X_k, y, family="binomial")
  
  # Fetaure selection for different target FDR levels
  
  for (i in 1:length(list_target_fdr)){
    target_fdr <- list_target_fdr[i]
    
    thres_LCD <- knockoff.threshold(W_LCD, fdr=target_fdr, offset=1) 
    selected_LCD <- which(W_LCD >= thres_LCD)
    
    thres_LSM <- knockoff.threshold(W_LSM, fdr=target_fdr, offset=1) 
    selected_LSM <- which(W_LSM >= thres_LSM)
    
    thres_LM <- knockoff.threshold(W_LM, fdr=target_fdr, offset=1) 
    selected_LM <- which(W_LM >= thres_LM)
    
    mat.features[1, (1+3*(i-1)):(3+3*(i-1))] <- c(fdp(selected_LSM, beta), fdp(selected_LCD, beta), fdp(selected_LM, beta))
    mat.features[2, (1+3*(i-1)):(3+3*(i-1))] <- c(power(selected_LSM, beta), power(selected_LCD, beta), power(selected_LM, beta))
    mat.features[selected_LSM + 2, 1+3*(i-1)] = 1
    mat.features[selected_LCD + 2, 2+3*(i-1)] = 1
    mat.features[selected_LM + 2, 3+3*(i-1)] = 1
    
    
    colnames(mat.features)[c(1+3*(i-1), 2+3*(i-1), 3+3*(i-1))] <- c(str_glue("LSM", target_fdr), str_glue("LCD", target_fdr), str_glue("LM", target_fdr))
    
  }
  
  rownames(mat.features) <- c("FDP", "power", c(1:p))
  
  if (is.null(y.copy)){
    return(list(mat.features, beta))
  } else {
    return(mat.features)
  }
  
}

Knockoff_TreeBasedMethod.synthesis <- function(X, X_k = NULL, y = NULL, beta = NULL, scale = NULL, k = NULL, amplitude = NULL, list_target_fdr, method = NULL){
  
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
  
  mat.features <- data.frame(matrix(0, ncol = 4*length(list_target_fdr), nrow = p+2 ))
  
  if (is.null(X_k)){# Copute the KO matrix
    X_k <- create.LSCIP_parallel(X)
  }
  
  if (is.null(y)){# Simulate the outcome
    beta <- beta.sample(p, k, amplitude)
    y <- y.sample(X, beta, method, scale, threshold = 0.5)
  }
  
  
  X_py <- numpy$ascontiguousarray(X)
  X_k_py <- numpy$ascontiguousarray(X_k)
  y_py <- numpy$ascontiguousarray(numpy$ravel(y))
  
  
  # GINI
  y.factor <- as.factor(y)
  W_GINI <- stat.random_forest(X, X_k, y.factor)
  
  
  # swap
  KO_swap <- knockpy$knockoff_stats$RandomForestStatistic()
  KO_swap$fit(X_py, X_k_py, y_py, feature_importance ="swap")
  W_swap <- KO_swap$W
  W_swap <- py_to_r(W_swap)
  
  # swapint
  KO_swapint <- knockpy$knockoff_stats$RandomForestStatistic()
  KO_swapint$fit(X_py, X_k_py, y_py, feature_importance ="swapint")
  W_swapint <- KO_swapint$W
  W_swapint <- py_to_r(W_swapint)
  
  # RRB
  W_RRB <- W.mboost_varimp2(X, X_k, y, max.mstop = 200, bl = c("btree"), k_fold = 10, family = Binomial())
  
  # SElect features for a sequence of target fdr levels
  for (i in 1:length(list_target_fdr)){
    
    target_fdr <- list_target_fdr[i]
    
    thres_GINI <- knockoff.threshold(W_GINI, fdr=target_fdr, offset=1) 
    selected_GINI <- which(W_GINI >= thres_GINI)
    
    thres_swap <- knockoff.threshold(W_swap, fdr=target_fdr, offset=1) 
    selected_swap <- which(W_swap >= thres_swap)
    
    thres_swapint <- knockoff.threshold(W_swapint, fdr=target_fdr, offset=1) 
    selected_swapint <- which(W_swapint >= thres_swapint)
    
    thres_RRB <- knockoff.threshold(W_RRB, fdr=target_fdr, offset=1) 
    selected_RRB <- which(W_RRB >= thres_RRB)
    
    
    
    mat.features[1, (1+4*(i-1)):(4+4*(i-1))] <- c(fdp(selected_swap, beta), fdp(selected_GINI, beta), fdp(selected_swapint, beta), fdp(selected_RRB, beta))
    mat.features[2, (1+4*(i-1)):(4+4*(i-1))] <- c(power(selected_swap, beta), power(selected_GINI, beta), power(selected_swapint, beta), power(selected_RRB, beta))
    mat.features[selected_swap + 2, 1+4*(i-1)] = 1
    mat.features[selected_GINI + 2, 2+4*(i-1)] = 1
    mat.features[selected_swapint + 2, 3+4*(i-1)] = 1
    mat.features[selected_RRB + 2, 4+4*(i-1)] = 1
    
    
    colnames(mat.features)[c(1+4*(i-1), 2+4*(i-1), 3+4*(i-1), 4+4*(i-1))] <- c(str_glue("swap", target_fdr), str_glue("GINI", target_fdr), str_glue("swapint", target_fdr), str_glue("RRB", target_fdr))
    
  }
  
  rownames(mat.features) <- c("FDP", "power", c(1:p))
  
  if (is.null(y.copy)){
    return(list(mat.features, beta))
  } else {
    return(mat.features)
  }
  
}

Knockoff_DeepBasedMethod <- function(X, X_k = NULL, y = NULL, beta = NULL, scale = NULL, k = NULL, amplitude = NULL, list_target_fdr, method = NULL){
  
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
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  y.copy = y
  
  mat.features <- data.frame(matrix(0, ncol = 4*length(list_target_fdr), nrow = p+2 ))
  
  if (is.null(X_k)){# Create KO matrix
    X_k <- create.LSCIP_parallel(X)
  }
  
  if (is.null(y)){# Simulate y
    beta <- beta.sample(p, k, amplitude)
    y <- y.sample(X, beta, method, scale, threshold = 0.5)
  }
  
  # Compute statistics
  
  # pre-processing to use py functions 
  X_py <- numpy$ascontiguousarray(X)
  X_k_py <- numpy$ascontiguousarray(X_k)
  y_py <- numpy$ascontiguousarray(numpy$ravel(y))
  
  # swap
  KO_swap <- knockpy$knockoff_stats$DeepPinkStatistic()
  KO_swap$fit(X_py, X_k_py, y_py, feature_importance ="swap")
  W_swap <- KO_swap$W
  W_swap <- py_to_r(W_swap)
  
  # swapint
  KO_swapint <- knockpy$knockoff_stats$DeepPinkStatistic()
  KO_swapint$fit(X_py, X_k_py, y_py, feature_importance ="swapint")
  W_swapint <- KO_swapint$W
  W_swapint <- py_to_r(W_swapint)
  
  
  # deeppink
  KO_deeppink <- knockpy$knockoff_stats$DeepPinkStatistic()
  KO_deeppink$fit(X_py, X_k_py, y_py, feature_importance ="deeppink")
  W_deeppink <- KO_deeppink$W
  W_deeppink <- py_to_r(W_deeppink)
  
  # unweighted
  KO_VI <- knockpy$knockoff_stats$DeepPinkStatistic()
  KO_VI$fit(X_py, X_k_py, y_py, feature_importance ="unweighted")
  W_VI <- KO_VI$W
  W_VI <- py_to_r(W_VI)
  

  
  # Features selection for different target fdr levels
  
  for (i in 1:length(list_target_fdr)){
    
    target_fdr <- list_target_fdr[i]
    
    thres_deeppink <- knockoff.threshold(W_deeppink, fdr=target_fdr, offset=1) 
    selected_deeppink <- which(W_deeppink >= thres_deeppink)
    
    thres_swap <- knockoff.threshold(W_swap, fdr=target_fdr, offset=1) 
    selected_swap <- which(W_swap >= thres_swap)
    
    thres_swapint <- knockoff.threshold(W_swapint, fdr=target_fdr, offset=1) 
    selected_swapint <- which(W_swapint >= thres_swapint)
    
    thres_VI <- knockoff.threshold(W_VI, fdr=target_fdr, offset=1) 
    selected_VI <- which(W_VI>= thres_VI)
    
    
    
    mat.features[1, (1+4*(i-1)):(4+4*(i-1))] <- c(fdp(selected_swap, beta), fdp(selected_deeppink, beta), fdp(selected_swapint, beta), fdp(selected_VI, beta))
    mat.features[2, (1+4*(i-1)):(4+4*(i-1))] <- c(power(selected_swap, beta), power(selected_deeppink, beta), power(selected_swapint, beta), power(selected_VI, beta))
    mat.features[selected_swap + 2, 1+4*(i-1)] = 1
    mat.features[selected_deeppink + 2, 2+4*(i-1)] = 1
    mat.features[selected_swapint + 2, 3+4*(i-1)] = 1
    mat.features[selected_VI + 2, 4+4*(i-1)] = 1
    
    
    colnames(mat.features)[c(1+4*(i-1), 2+4*(i-1), 3+4*(i-1), 4+4*(i-1))] <- c(str_glue("swap", target_fdr), str_glue("deeppink", target_fdr), str_glue("swapint", target_fdr), str_glue("VI", target_fdr))
    
  }
  
  rownames(mat.features) <- c("FDP", "power", c(1:p))
  
  if (is.null(y.copy)){
    return(list(mat.features, beta))
  } else {
    return(mat.features)
  }
  
}


########################################################################################################################################


                                       ### Functions for graphical representations ###


df_perf_maker <- function(table.perf, list_target_fdr, stats_names, n.rep, list_lambda=NULL){
  
  # table.perf: performance tables with FDP and power (rows 1 and 2)
  # list_target_fdr: target fdr used in the experiments to generate table.perf
  # stats_names: names of the different statistics used to create table.perf
  # number of times the comparison has been done 
  # list_lambda: lambda values used to perform variable selection with the LASSO
  
  n.stats = length(stats_names)
  compteur = 0
  
  len_target_fdr = length(list_target_fdr)
  
  if ("LPLR" %in% stats_names){ # Account for LASSO in the comparison (not use in the paper)
    index = which(stats_names=="LPLR" | stats_names=="LASSO")
    target_fdr = rep(list_target_fdr, n.stats)
    target_fdr[((index-1)*len_target_fdr + 1) : (index*len_target_fdr)] = list_lambda
  }else{target_fdr = rep(list_target_fdr, n.stats)}
  
  df_power_fdr <- data.frame(power = rep(NA, n.stats*len_target_fdr), target_fdr = target_fdr, std = rep(NA, n.stats*len_target_fdr), KO = rep(stats_names, each = len_target_fdr))
  df_power_fdr$KO <- as.vector(df_power_fdr$KO)
  
  df_fdp_fdr <- data.frame(fdp = rep(NA, n.stats*len_target_fdr), target_fdr = target_fdr , std = rep(NA, n.stats*len_target_fdr), KO = rep(stats_names, each = len_target_fdr))
  df_fdp_fdr$KO <- as.vector(df_fdp_fdr$KO)
  
  
  list_power_bin = c()
  list_power_sd_bin = c()
  list_FDP_bin = c()
  list_rep_bin = c()
  
  
  
  for (i in 1:n.stats){
    compteur = compteur+1
    
    
    ### creating df for performances
    
    # fdr curves
    
    table = table.perf
    list_FDP_perf = as.numeric(table[1, seq(i, len_target_fdr*n.stats*n.rep, n.stats)])
    list_power_perf = as.numeric(table[2, seq(i, len_target_fdr*n.stats*n.rep, n.stats)])
    
    list_power_mean = c()
    list_FDP_mean = c()
    
    for  (j in 1:len_target_fdr){
      list_power_mean = c(list_power_mean, mean(list_power_perf[seq(j,len_target_fdr*n.rep,len_target_fdr)]))
      list_FDP_mean = c(list_FDP_mean, mean(list_FDP_perf[seq(j,len_target_fdr*n.rep,len_target_fdr)]))
    }
    
    list_power_sd = c()
    list_FDP_sd = c()
    
    for (j in 1:len_target_fdr){
      list_power_sd = c(list_power_sd, sd(list_power_perf[seq(j,len_target_fdr*n.rep,len_target_fdr)]))
      list_FDP_sd = c(list_FDP_sd, sd(list_FDP_perf[seq(j,len_target_fdr*n.rep,len_target_fdr)]))
    }
    
    df_power_fdr$power[((compteur-1)*len_target_fdr + 1) : (compteur*len_target_fdr)] = list_power_mean
    df_power_fdr$std[((compteur-1)*len_target_fdr + 1) : (compteur*len_target_fdr)] = list_power_sd
    
    df_fdp_fdr$fdp[((compteur-1)*len_target_fdr + 1) : (compteur*len_target_fdr)] = list_FDP_mean
    df_fdp_fdr$std[((compteur-1)*len_target_fdr + 1) : (compteur*len_target_fdr)] = list_FDP_sd
    
    # FDP curve
    
    list_perf_m = bin(list_FDP_perf, list_power_perf, 0.05) # mean power over an interval of FDP values

    list_sd = list_perf_m[[3]]
    list_sd[which(is.na(list_sd))]=0

    list_power_bin = c(list_power_bin, list_perf_m[[2]])
    list_power_sd_bin = c(list_power_sd_bin, list_sd)
    list_FDP_bin = c(list_FDP_bin, list_perf_m[[1]])
    list_rep_bin = c(list_rep_bin, length(list_perf_m[[1]]))
    
  }
  
  
  df_power_fdp = data.frame(power = list_power_bin, fdp = list_FDP_bin, std = list_power_sd_bin, KO = rep(stats_names, list_rep_bin))
  df_power_fdp$KO <- as.factor(df_power_fdp$KO)
  
  
  
  return(list(df_power_fdr, df_fdp_fdr, df_power_fdp))
}
bin <- function(list_x, list_y, length_bin){# create bins for Power against FDP curves
  
  # list_x : list of FDP values
  # list_y : list of power values
  # length_bin : FDP intervals used to compute the average power
  
  n <- length(list_x)
  l_min <- min(list_x)
  l_max <- max(list_x)
  
  if (l_min == l_max){bin_number = 1}
  else{
    bin_number = ceiling((l_max-l_min)/length_bin)
  }
  list_x_m = c()
  list_y_m = c()
  list_y_sd = c()
  for (i in 1:bin_number){ # defining bin
    bound_inf = l_min+(i-1)*length_bin
    bound_sup = l_min+i*length_bin
    
    if (bound_sup == l_max){# accounting for bound issues
      l_y = which(list_x <= bound_sup & list_x >= bound_inf)
    }
    else{l_y = which(list_x < bound_sup & list_x >= bound_inf)}
    
    list_y_m = c(list_y_m, mean(list_y[l_y]))
    if (length(list_y[l_y])==1){
      list_y_sd = c(list_y_sd, 0)
    }
    else{
      list_y_sd = c(list_y_sd, sd(list_y[l_y]))
    }
    if (bound_sup > l_max){bound_sup = l_max}
    
    list_x_m = c(list_x_m,bound_inf+(bound_sup-bound_inf)/2)
  }
  i = intersect(which(is.na(list_y_m)),which(is.na(list_y_sd)))
  if (length(intersect(which(is.na(list_y_m)),which(is.na(list_y_sd))))!=0){ # avoid empty bin issues
    list_x_m = list_x_m[-i]
    list_y_m = list_y_m[-i]
    list_y_sd = list_y_sd[-i] 
  }
  return(list(list_x_m, list_y_m, list_y_sd))
}
