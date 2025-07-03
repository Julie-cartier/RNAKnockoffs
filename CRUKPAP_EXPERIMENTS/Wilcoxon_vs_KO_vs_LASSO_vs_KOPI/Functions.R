# y simulation

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

# Performance metrics

fdp <- function(selected, beta_coeff) + return(sum(beta_coeff[selected] == 0) / max(1, length(selected)))

power <- function(selected, beta_coeff) + return(length(which(beta_coeff[selected]!=0))/max(1, length(which(beta_coeff != 0))))

# Create KO

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

# BH procedure for a vector of p_values

BH_procedure <- function(pvalues, target_fdr = 0.05){
  
  # pvalues : vector of pvalues
  # target_fdr : target fdr 
  
  p <- length(pvalues)
  
  # Order the p-values 
  
  p_order <- sort(pvalues, index.return = TRUE)
  p_sorted <- p_order$x      # Sorted p-values
  p_indices <- p_order$ix     # Original indices of the sorted p-values
  
  # Calculate the threshold for each p-value
  
  BH_threshold <- (1:p) * target_fdr / p 
  
  # Select features
  selected <- which(p_sorted <= BH_threshold)
  
  if (length(selected) > 0) {
    
    max_selected <- max(selected)
    return(p_indices[1:max_selected])
    
  } else {
    return(integer(0)) 
  }
}

# Function to compute KOPI

create_KO_stats <- function(X, list_X_k, y, num_cores = 10) {# Save D statistics vector for aggregation  (p x draws matrix)
  
  # X : data matrix
  # list_X_k : list of D computed KO matrices 
  # y : outcomes
  
  Draws <- length(list_X_k)

  KO_stats_list <- mclapply(list_X_k, function(X_k) {
    stat.glmnet_coefdiff(X, X_k, y, nfolds = 10, family = "binomial")
  }, mc.cores = num_cores)
  
  # Combine the list of results into a matrix
  KO_stats_mat <- do.call(rbind, KO_stats_list)
  
  return(KO_stats_mat)
}

# The following functions return tables with, for each method (LPLR/KO/KOPI) and for each target FDR/FDP/lambda, the power, the FDP and the set of selected features. 

# Selection with the KO framework and the lasso penalized logistic regression

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

# Differential expression analysis with Wilcoxon rank-sum test 

DE_Wilcoxon_BH <- function(X, y = NULL, beta = NULL, scale = FALSE, k = NULL, amplitude = NULL, list_target_fdr, method = NULL){
  
  # X: data matrix of size nxp
  # y: outcomes (vector of size n)
  # beta: vector of size p (coefficients)
  # scale: If TRUE data are scaled (Default to FALSE)
  # k: number of non null coefficients (non-null features). 
  # amplitude: amplitude of the non null coefficients 
  # list_target_fdr : vector of target FDR for the selection with the BH procedure
  # method: either "linear" or "interaction" (to generate y)
  
  p <- ncol(X)
  n <- nrow(X)
  y.copy <- y
  
  if (is.null(y)){
    beta = beta.sample(p, k, amplitude)
    y = y.sample(X, beta, method, scale, threshold = 0.5)
  }
  
  mat.features <- data.frame(matrix(0, ncol = length(list_target_fdr)*1, nrow = p+2 ))
  
  ### get p values with the wilcoxon test ###
  
  pvalues <- sapply(1:ncol(X),function(i){
    data <- cbind.data.frame(gene = as.numeric(X[,i]), conditions = factor(y))
    p = wilcox.test(gene~conditions, data)$p.value
  })
  
  compteur = 0
  # Selection for different target FDR
  
  for (alpha in list_target_fdr){
    
    compteur = compteur + 1
    
    # BH 
    
    p_adjust <- p.adjust(pvalues, method = "fdr")
    selected_BH <- which(p_adjust < alpha)
    
    mat.features[1, compteur] <- c(fdp(selected_BH, beta))
    mat.features[2, compteur] <- c(power(selected_BH, beta))
    mat.features[selected_BH + 2, compteur] <- 1
 
    
    colnames(mat.features)[compteur] <- str_glue("BH", alpha)
  }
  
  rownames(mat.features) <- c("FDP", "power", c(1:p))
  
  if (is.null(y.copy)){
    return(list(mat.features, beta))
  }
  else{
    return(mat.features)
  }
}

# Selection with the KO framework and the lasso penalized logistic regression where the penalization parameter is lambda oracle

Add_KOPI_LPLR_oracle <- function(X, D_draws_X_k, y, beta, list_target_FDP, k_max, alpha = 0.1,  B = 1000){
  
  # X: data matrix of size nxp
  # D_draw_X_k : number of KO matrices used in the aggregation scheme
  # y: outcomes vector 
  # beta: vector of size p (coefficients)
  # list_target_FDP : FDP targets for KOPI
  # k_max : KOPI parameter (number of thresholds in the threshold family (KOPI))
  # alpha : level of confidence (KOPI)
  # B : KOPI parameter (Number of Monte Carlo simulations (KOPI))
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  mat.features <- data.frame(matrix(0, ncol = length(list_target_FDP), nrow = p+2 ))
  
  ### KOPI Settings ###
  
  Draws <- length(D_draws_X_k)
  
  # create KO_stats (matrix of size D x p)
  
  KO_stats <- create_KO_stats(X, D_draws_X_k, y, num_cores = 10)
  KO_stats_py <- r_to_py(KO_stats)
  
  # get pi0_hmean and learned_tpl_hmean
  
  result <- create_hmean_p_values_learned_tpl(p = as.integer(p), draws = as.integer(Draws), B = as.integer(B), n_jobs = as.integer(10))
  
  pi0_hmean <- result[[1]]
  learned_tpl_hmean <- result[[2]]
  
  
  for (i in 1:length(list_target_FDP)){
    
    q <- list_target_FDP[i]
    
    set_selected <- KOPI_given_all(ko_stats = KO_stats_py, draws = Draws, alpha = alpha, target_fdp = q, learned_tpl_hmean = learned_tpl_hmean, pi0_hmean = pi0_hmean, k_max = as.integer(k_max))
    selected_KOPI <- set_selected[[1]] + 1
    
    
    # Data frame creator
    
    mat.features[1, i] <- fdp(selected_KOPI, beta)
    mat.features[2, i] <- power(selected_KOPI, beta)
    mat.features[selected_KOPI + 2, i] <- 1
    
    
    
    colnames(mat.features)[i] =  str_glue("KOPI_q_", q ,"_")
    
  }
  
  rownames(mat.features) = c("FDP", "power", c(1:p))
  

  ## get lambda oracle 
  
  new_lambda_sequence <- seq(1, 0, -0.001)
  fitlog_oracle <- glmnet(X, y, type.measure = "deviance", family = "binomial", lambda = seq(1, 0, -0.001))
  
  vec.ratio.TD.FD = c()
  
  # Find the best lambda
  
  for (lambda in new_lambda_sequence){
    coef_lasso <- coef(fitlog_oracle, s = lambda)
    selected_lasso <- which(coef_lasso[-1]!=0)
    vec.ratio.TD.FD <- c(vec.ratio.TD.FD, sum(beta[selected_lasso] != 0)/sum(beta[selected_lasso] == 0))
  }
  
  # get the power and FDP corresponding to this lambda 
  
  lambda.oracle <- new_lambda_sequence[which.max(vec.ratio.TD.FD)]
  beta.lambda.oracle <- coef(fitlog_oracle, s = lambda.oracle)
  selected_lambda.oracle <- which(beta.lambda.oracle[-1]!=0)
  vec.lambda.oracle <- c(lambda.oracle, fdp(selected_lambda.oracle, beta), power(selected_lambda.oracle, beta))
  
  
  return(list(mat.features, vec.lambda.oracle))
  
}

# Performance comparison: KO vs KOPI

KOPI_KO_comp <- function(X, D_draws_X_k, y = NULL, beta = NULL, scale = FALSE, k = NULL, amplitude = NULL, alpha = 0.1, list_target_fdr_FDP, method = NULL, B = 100, k_max = NULL){
  
  # X: data matrix of size nxp
  # D_draw_X_k : number of KO matrices used in KOPI
  # y: outcomes (vector of size n)
  # beta: vector of size p (coefficients)
  # scale: If TRUE data are scaled (Default to FALSE)
  # k: number of non null coefficients (non-null features). To simulate interacting settings, k is assumed to be even.
  # amplitude: amplitude of the non null coefficients 
  # alpha : level of confidence (KOPI)
  # list_target_fdr_FDP : vector of target FDR/FDP for the selection with the KO framework/KOPI
  # method: either "linear" or "interaction"
  # D_draw_X_k : number of KO matrices 
  # list_target_fdr_FDP : FDP targets for KOPI, and fdr targets for KO
  # k_max : KOPI parameter (number of thresholds in the threshold family (KOPI))
  # B : KOPI parameter (NUmber of Monte Carlo simulations (KOPI))
  
  p <- ncol(X)
  n <- nrow(X)
  y.copy = y
  
  if (is.null(y)){#simulate y
    beta <- beta.sample(p, k, amplitude)
    y <- y.sample(X, beta, method, scale, threshold = 0.5)
  }
  
  mat.features <- data.frame(matrix(0, ncol = length(list_target_fdr_FDP), nrow = p+2 ))
  
  
  
  ### KOPI Settings ###
  
  Draws <- length(D_draws_X_k)
  
  # create KO_stats (matrix of size D x p)
  
  KO_stats <- create_KO_stats(X, D_draws_X_k, y, num_cores = 10)
  KO_stats_py <- r_to_py(KO_stats)
  
  # get pi0_hmean and learned_tpl_hmean
  
  result <- create_hmean_p_values_learned_tpl(p = as.integer(p), draws = as.integer(Draws), B = as.integer(B), n_jobs = as.integer(10))
  
  pi0_hmean <- result[[1]]
  learned_tpl_hmean <- result[[2]]
  
  
  
  ###  Vanilla KO ###
  
  X_LSCIP <- D_draws_X_k[[1]]
  W_LSCIP <- stat.glmnet_coefdiff(X, X_LSCIP, y, nfolds=10, family="binomial")
  
  
  for (i in 1:length(list_target_fdr_FDP)){
    
    #  KO selection 
    
    target_fdr <- list_target_fdr_FDP[i]
    thres_LSCIP <- knockoff.threshold(W_LSCIP, fdr=target_fdr, offset=1) 
    selected_LSCIP <- which(W_LSCIP >= thres_LSCIP)
    
    # KOPI Selection
    
    q <- list_target_fdr_FDP[i]
    
    
    set_selected <- KOPI_given_all(ko_stats = KO_stats_py, draws = Draws, alpha = alpha, target_fdp = q, learned_tpl_hmean = learned_tpl_hmean, pi0_hmean = pi0_hmean, k_max = as.integer(k_max))
    selected_KOPI <- set_selected[[1]] + 1
    
    
    # Data frame creator
    
    mat.features[1, (1+2*(i-1)):(2+2*(i-1))] <- c(fdp(selected_LSCIP, beta), fdp(selected_KOPI, beta))
    mat.features[2, (1+2*(i-1)):(2+2*(i-1))] <- c(power(selected_LSCIP, beta), power(selected_KOPI, beta))
    mat.features[selected_LSCIP + 2, 1+2*(i-1)] <- 1
    mat.features[selected_KOPI + 2, 2+2*(i-1)] <- 1
    
    
    
    colnames(mat.features)[c(1+2*(i-1), 2+2*(i-1)) ] <- c(str_glue("Vanilla KO", target_fdr), str_glue("KOPI", target_fdr))
    
  }
  
  rownames(mat.features) <- c("FDP", "power", c(1:p))
  
  if (is.null(y.copy)){
    return(list(mat.features, beta))
  }
  else{
    return(mat.features)
  }
}

################################################################################################################
                          #### Functions for graphical representation ####

df_perf_maker <- function(table.perf, list_target_fdr, stats_names, n.rep, list_lambda=NULL){
  
  # table.perf : data frame where the two first rows give the FDP and the power of the selection (experiment result)
  # list_target_fdr : vector of targets fdr used to create table.perf
  # stats_names : give the names of the methods compared (must be given in the same order as in table.perf)
  # n.rep : number of repetition of the experiment
  # list_lambda : vector with lambda values used in the experiment

  n.stats = length(stats_names)
  compteur = 0
  
  len_target_fdr = length(list_target_fdr)
  
  if ("LASSO" %in% stats_names | "LPLR" %in% stats_names){
    index = which(stats_names=="LASSO" | stats_names == "LPLR")
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
    
    # FDP curves 
    
    list_perf_m = bin(list_FDP_perf, list_power_perf, 0.05)
    
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



########################################################################################################################################


                                        ### Functions for graphical representations ###



bin <- function(list_x, list_y, length_bin){ # For the power vs FDP curves
  
  # list_x : list of FDP values
  # list_y : list of power values
  # length_bin : FDP intervals used to compute the average power
  
  n = length(list_x)
  l_min = min(list_x)
  l_max = max(list_x)
  if (l_min == l_max){bin_number = 1}
  else{
    bin_number = ceiling((l_max-l_min)/length_bin)
  }
  list_x_m = c()
  list_y_m = c()
  list_y_sd = c()
  for (i in 1:bin_number){
    bound_inf = l_min+(i-1)*length_bin
    bound_sup = l_min+i*length_bin
    
    if (bound_sup == l_max){
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
  if (length(intersect(which(is.na(list_y_m)),which(is.na(list_y_sd))))!=0){
    list_x_m = list_x_m[-i]
    list_y_m = list_y_m[-i]
    list_y_sd = list_y_sd[-i] 
  }
  return(list(list_x_m, list_y_m, list_y_sd))
}

df_perf_maker_boxplot <- function(table.perf, list_target_fdr, stats_names, n.rep, length_bin = 0.05, FDP.limits = c(0,1)){ 
  
  # table.perf : data frame where the two first rows give the FDP and the power of the selection
  # list_target_fdr : vector with sequence of target fdr for KO selection
  # stats_names : give the names of the methods compared (must be given in the same order as the table.perf)
  # n.rep : number of repetition of the comparison 
  # length_bin : bin size for the boxplots (corresponding to FDP intervals)
  # FDP.limits : Interval of FD values that we consider
  
  n.stats = length(stats_names)
  compteur = 0
  len_target_fdr  = length(list_target_fdr)
  
  df_power_fdp_list <- list()
  
  
  for (i in 1:n.stats){
    compteur = compteur+1
    
    
    
    table = table.perf
    list_FDP_perf = as.numeric(table[1, seq(i, len_target_fdr*n.stats*n.rep, n.stats)])
    list_power_perf = as.numeric(table[2, seq(i, len_target_fdr*n.stats*n.rep, n.stats)])
    
    breaks <- seq(min(FDP.limits), max(FDP.limits), length_bin)
    FDP_levels <- cut(list_FDP_perf, breaks = breaks, right = FALSE, include.lowest = TRUE)
    
    
    df_part <- data.frame(FDP = list_FDP_perf, Power = list_power_perf, FDP_levels = as.character(FDP_levels), KO = stats_names[i])
    
    df_power_fdp_list[[compteur]] <- df_part
  }
  df_power_fdp <- do.call(rbind, df_power_fdp_list)
  df_power_fdp$KO <- as.factor(df_power_fdp$KO)
  
  
  
  return(df_power_fdp)
}


