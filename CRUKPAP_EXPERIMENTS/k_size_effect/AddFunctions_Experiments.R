#get settings and usual functions

source("/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Statistics_comparison/Settings.R")
source("/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Statistics_comparison/Functions.R")

source("/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Wilcoxon_vs_KO_vs_LASSO_vs_KOPI/Functions.R")

# Additional setting 

# Target FDR for the BH procedure (post wilcoxon test)
exponents <- -10:0
vec <- c(10**exponents, 5*10**exponents)
list_target_fdr_wilcoxon <- vec[-length(vec)]
list_target_fdr_wilcoxon <- list_target_fdr_wilcoxon[order(list_target_fdr_wilcoxon)]

# Additional function

Knockoff_LPLR <- function(X, X_k = NULL, y = NULL, beta = NULL, scale = NULL, k = NULL, amplitude = NULL, list_target_fdr, list_lambda, method = NULL){
  
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
  
  mat.features <- data.frame(matrix(0, ncol = 2*length(list_target_fdr), nrow = p+2 ))
  
  if (is.null(X_k)){ # compute KO matrix
    X_k <- create.LSCIP_parallel(X)
  }
  
  if (is.null(y)){ # simulate y
    beta <- beta.sample(p, k, amplitude)
    y <- y.sample(X, beta, method, scale, threshold = 0.5)
  }
  
  # LCD statistics
  if (method == "linear"){
    W_KO <- stat.glmnet_coefdiff(X, X_k, y, nfolds=10, family="binomial")
  }
  else if(method == "interaction"){
    # GINI
    y.factor <- as.factor(y)
    W_KO <- stat.random_forest(X, X_k, y.factor)
  }
  
  
  # Select features
  
  for (i in 1:length(list_target_fdr)){
    
    target_fdr <- list_target_fdr[i]
    lambda <- list_lambda[i]
    
    # KO
    
    thres_KO <- knockoff.threshold(W_KO, fdr=target_fdr, offset=1) 
    selected_KO <- which(W_KO >= thres_KO)
    
    #LASSO
    
    fitlog <- glmnet(X, y, type.measure = "deviance", family = "binomial", lambda = sort(list_lambda, decreasing = TRUE))
    coef_lasso <- coef(fitlog, s = lambda)
    selected_lasso <- which(coef_lasso[-1]!=0)
    
    predicted_y  <- predict(fitlog, s = lambda, X, type='response')
    predicted_y[predicted_y>0.5]  <- 1
    predicted_y[predicted_y<=0.5] <- 0
    
    mat.features[1, (1+2*(i-1)):(2+2*(i-1))] <- c(fdp(selected_KO, beta), fdp(selected_lasso, beta))
    mat.features[2, (1+2*(i-1)):(2+2*(i-1))] <- c(power(selected_KO, beta), power(selected_lasso, beta))
    mat.features[selected_KO + 2, 1+2*(i-1)] = 1
    mat.features[selected_lasso + 2, 2+2*(i-1)] = 1
    
    if (method == "linear"){
      colnames(mat.features)[c(1+2*(i-1), 2+2*(i-1)) ] <- c(str_glue("LCD", target_fdr), str_glue("lasso_", lambda))
    }
    else if (method == "interaction"){
      colnames(mat.features)[c(1+2*(i-1), 2+2*(i-1)) ] <- c(str_glue("RF", target_fdr), str_glue("lasso_", lambda))
    }
  }
  
  rownames(mat.features) = c("FDP", "power", c(1:p))
  
  # get the FDP and power for lambda.min : 
  
  cv.fitlog <- cv.glmnet(X, y, type.measure = "deviance", family = "binomial")
  beta.lambda.min <- coef(cv.fitlog, s = "lambda.min")
  selected_lambda.min <- which(beta.lambda.min[-1]!=0)
  
  vec.lambda.min <- c(cv.fitlog$lambda.min, fdp(selected_lambda.min, beta), power(selected_lambda.min, beta))
  
  
  if (is.null(y.copy)){
    return(list(mat.features, beta, vec.lambda.min))
  }
  else{
    return(list(mat.features, vec.lambda.min))
  }
  
}



KO_ksize_effect <- function(X, X_k = NULL, vec_k, amplitude, list_target_fdr, list_lambda,  method){
  
  # X: data matrix of size nxp
  # X_k: X KO matrix
  # vec_k : vector of k values (number of non null coefficients (non-null features)). To simulate interacting settings, k is assumed to be even.
  # amplitude: amplitude of the non null coefficients 
  # list_target_fdr : vector of target FDR for the selection with the KO framework
  # method: either "linear" or "interaction"
  
  n <- nrow(X)
  p <- ncol(X)
  
  
  k_length <- length(vec_k)
  mat.features <- data.frame(matrix(0, ncol = k_length*length(list_target_fdr), nrow = p+2 ))
  mat.features.LPLR <- data.frame(matrix(0, ncol = k_length*length(list_target_fdr), nrow = p+2 ))
  
  if (is.null(X_k)){
    X_k <- create.LSCIP_parallel(X)
  }
  
  compteur = 0
  for (k in vec_k){
    
    compteur = compteur + 1
    
    #simulate y 
    beta <- beta.sample(p, k, amplitude)
    y <- y.sample(X, beta, method, scale = FALSE, threshold = 0.5)
    
    if (k == 30){
      beta.k30 = beta
      y.k30 = y
    }
    # compute the LCD statistic
    W_LCD <- stat.glmnet_coefdiff(X, X_k, y, nfolds=10, family="binomial")
    
    # get sets of selected features for different targets FDR level
    for (i in 1:length(list_target_fdr)){
      target_fdr <- list_target_fdr[i]
      lambda <- list_lambda[i]
      
      thres_LCD <- knockoff.threshold(W_LCD, fdr=target_fdr, offset=1) 
      selected_LCD <- which(W_LCD >= thres_LCD)
      
      #LASSO
      
      fitlog <- glmnet(X, y, type.measure = "deviance", family = "binomial", lambda = sort(list_lambda, decreasing = TRUE))
      coef_lasso <- coef(fitlog, s = lambda)
      selected_lasso <- which(coef_lasso[-1]!=0)
      
      mat.features[1, compteur + k_length*(i-1)] <- fdp(selected_LCD, beta)
      mat.features[2, compteur + k_length*(i-1)] <- power(selected_LCD, beta)
      mat.features[selected_LCD + 2, compteur + k_length*(i-1)] <- 1
      
      mat.features.LPLR[1, compteur + k_length*(i-1)] <- fdp(selected_lasso, beta)
      mat.features.LPLR[2, compteur + k_length*(i-1)] <- power(selected_lasso, beta)
      mat.features.LPLR[selected_lasso + 2, compteur + k_length*(i-1)] <- 1
      
      colnames(mat.features)[compteur + k_length*(i-1)] <- c(str_glue("k_", k, "_tFDR_", target_fdr,"_rep_"))
      colnames(mat.features.LPLR)[compteur + k_length*(i-1)] <- c(str_glue("k_", k, "_lambda_", lambda, "_rep_"))
      
    }
  }
  return(list(mat.features, mat.features.LPLR, beta.k30, y.k30))
}

# Experiments 


table.perf.linear.ksize.KO <- data.frame(rep(NA, ncol(X_real)+2))
table.perf.linear.ksize.LPLR <- data.frame(rep(NA, ncol(X_real)+2))
table.perf.linear.k30.wilcoxon <- data.frame(rep(NA, ncol(X_real)+2))



tic()
for (i in 1:10){
  print(i)
  table.perf.part <- KO_ksize_effect(scale(X_real), vec_k = c(10,30,50), amplitude = 10, list_target_fdr = list_target_fdr, list_lambda = list_lambda,  method = "linear")  
  table.perf.linear.ksize.KO <- cbind(table.perf.linear.ksize.KO, table.perf.part[[1]])
  table.perf.linear.ksize.LPLR <- cbind(table.perf.linear.ksize.LPLR, table.perf.part[[2]])
  # # For the comparison with k = 30 of all methods
  table.perf.part.wilcoxon <- DE_Wilcoxon_BH(scale(X_real), y = table.perf.part[[4]], beta = table.perf.part[[3]], list_target_fdr = list_target_fdr_wilcoxon)
  table.perf.linear.k30.wilcoxon <- cbind(table.perf.linear.k30.wilcoxon, table.perf.part.wilcoxon)
}

toc()

table.perf.linear.ksize.KO <- table.perf.linear.ksize.KO[, -1]
table.perf.linear.ksize.LPLR <- table.perf.linear.ksize.LPLR[, -1]
table.perf.linear.k30.wilcoxon <- table.perf.linear.k30.wilcoxon[, -1]

save(table.perf.linear.ksize.KO, table.perf.linear.ksize.LPLR, table.perf.linear.k30.wilcoxon, file = "/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/k_size_effect/R_files/ksize_trueX_fakeY_tables_perfx10.R")

