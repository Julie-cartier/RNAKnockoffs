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
  # method: either "linear" or "interaction". 
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


# Linear variant of the sequential conditional independant pairs

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

# Comparison between : LSCIP/MVR/CI/SDP
# It returns a table with, for each method and for each target FDR, the power, the FDP and the set of selected features. 

Knockoff_comp_tot <- function(X, scale = FALSE, k, amplitude, list_target_fdr, method){
  
  # X: data matrix of size nxp
  # scale: If TRUE data are scaled (Default to FALSE)
  # k: number of non null coefficients (non-null features). To simulate interacting settings, k is assumed to be even.
  # amplitude: amplitude of the non null coefficients 
  # list_target_fdr : vector of target FDR for the selection with the KO framework
  # method: either "linear" or "interaction" (to generate y)
  
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  mat.features <- data.frame(matrix(0, ncol = 4*length(list_target_fdr), nrow = p+2 ))
  
  # Simulate linear outcomes 
  
  beta <- beta.sample(p, k, amplitude)
  y <- y.sample(X, beta, method, scale, threshold = 0.5)
  
  # create LSCIP KO
  
  X_LSCIP <- create.LSCIP_parallel(X)
  
  # create MVR KO
  
  knock_sampler <- knockpy$knockoffs$GaussianSampler(X, method = 'mvr')
  X_MVR <- knock_sampler$sample_knockoffs()
  X_MVR <- py_to_r(X_MVR)
  
  #create CI KO
  
  knock_sampler <- knockpy$knockoffs$GaussianSampler(X, method = 'ci')
  X_CI <- knock_sampler$sample_knockoffs()
  X_CI <- py_to_r(X_CI)
  
  # create SDP KO 
  knock_sampler <- knockpy$knockoffs$GaussianSampler(X, method = 'sdp')
  X_SDP <- knock_sampler$sample_knockoffs()
  X_SDP <- py_to_r(X_SDP)
  
  # Compute LCD statistics 
  
  W_MVR <- stat.glmnet_coefdiff(X, X_MVR, y, nfolds=10, family="binomial")
  W_SDP <- stat.glmnet_coefdiff(X, X_SDP, y, nfolds=10, family="binomial")
  W_CI <- stat.glmnet_coefdiff(X, X_CI, y, nfolds=10, family="binomial")
  W_LSCIP <- stat.glmnet_coefdiff(X, X_LSCIP, y, nfolds=10, family="binomial")
  
  # Selection for different target fdr 
  
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
    mat.features[selected_SDP + 2, 1+4*(i-1)] = 1
    mat.features[selected_MVR + 2, 2+4*(i-1)] = 1
    mat.features[selected_LSCIP + 2, 3+4*(i-1)] = 1
    mat.features[selected_CI + 2, 4+4*(i-1)] = 1
    
    
    colnames(mat.features)[c(1+4*(i-1), 2+4*(i-1), 3+4*(i-1), 4+4*(i-1)) ] = c(str_glue("SDP", target_fdr), str_glue("MVR", target_fdr), str_glue("LSCIP", target_fdr), str_glue("CI", target_fdr))
    
  }
  
  rownames(mat.features) = c("FDP", "power", c(1:p))
  
  
  return(list(mat.features, beta))
}

# Comparison of the performance obtained on gaussian data (with known or estimated covariance matrix)

Knockoff_comp_KO_COV <- function(X, Sigma.sim, scale = FALSE, k, amplitude, list_target_fdr, method){
  
  # X : data matrix
  # Sigma.sim : X true covariance matrix 
  # scale : if TRUE X is standardized (Default to FALSE)
  # k : number of non-null features
  # amplitude: amplitude of the non null coefficients 
  # list_target_fdr : vector of target FDR 
  # method: either "linear" or "interaction" (to generate y)
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  mat.features <- data.frame(matrix(0, ncol = 8*length(list_target_fdr), nrow = p+2))
  
  # Simulate outcomes
  
  beta <- beta.sample(p, k, amplitude)
  y <- y.sample(X, beta, method, scale, threshold = 0.5)
  
  # create KO 
  
  # SDP 
  
  knock_sampler <- knockpy$knockoffs$GaussianSampler(X, method = 'sdp') # Estimated covariance matrix
  X_1 <- knock_sampler$sample_knockoffs()
  X_1 <- py_to_r(X_1)
  knock_sampler <- knockpy$knockoffs$GaussianSampler(X, Sigma = Sigma.sim, method = 'sdp') # True covariance matrix
  X_2 <- knock_sampler$sample_knockoffs()
  X_2 <- py_to_r(X_2)

 # MVR
  
  knock_sampler <- knockpy$knockoffs$GaussianSampler(X, method = 'mvr') # Estimated covariance matrix
  X_3 <- knock_sampler$sample_knockoffs()
  X_3 <- py_to_r(X_3)
  knock_sampler <- knockpy$knockoffs$GaussianSampler(X, Sigma = Sigma.sim, method = 'mvr')  # True covariance matrix
  X_4 <- knock_sampler$sample_knockoffs()
  X_4 <- py_to_r(X_4)
  
  # CI
  
  knock_sampler <- knockpy$knockoffs$GaussianSampler(X, method = 'ci') # Estimated covariance matrix
  X_5 <- knock_sampler$sample_knockoffs()
  X_5 <- py_to_r(X_5)
  knock_sampler <- knockpy$knockoffs$GaussianSampler(X, Sigma = Sigma.sim, method = 'ci')  # True covariance matrix
  X_6 <- knock_sampler$sample_knockoffs()
  X_6 <- py_to_r(X_6)
  
  # Not use in the paper ( safety check : the results are not caused by the algorithm used to estimate the covariance matrix)
  # SDP (with R package)
  
  X_7 <- create.second_order(X, method = "sdp", shrink=T) # Estimated covariance matrix
  X_8 <- create.gaussian(X, mu = rep(0,p), Sigma = Sigma.sim, method = 'sdp')# True covariance matrix
  
  # Compute LCD statistics for each case
  
  W_SDP <- stat.glmnet_coefdiff(X, X_1, y, nfolds=10, family="binomial")
  W_SDP_T <- stat.glmnet_coefdiff(X, X_2, y, nfolds=10, family="binomial")
  W_MVR <- stat.glmnet_coefdiff(X, X_3, y, nfolds=10, family="binomial")
  W_MVR_T <- stat.glmnet_coefdiff(X, X_4, y, nfolds=10, family="binomial")
  W_CI <- stat.glmnet_coefdiff(X, X_5, y, nfolds=10, family="binomial")
  W_CI_T <- stat.glmnet_coefdiff(X, X_6, y, nfolds=10, family="binomial")
  W_SDPR <- stat.glmnet_coefdiff(X, X_7, y, nfolds=10, family="binomial")
  W_SDPR_T <- stat.glmnet_coefdiff(X, X_8, y, nfolds=10, family="binomial")
  
  # feature selection
  
  for (i in 1:length(list_target_fdr)){
    target_fdr <- list_target_fdr[i]
    
    thres_SDP <- knockoff.threshold(W_SDP, fdr=target_fdr, offset=1) 
    selected_SDP <- which(W_SDP>= thres_SDP)
    
    thres_SDP_T <- knockoff.threshold(W_SDP_T, fdr=target_fdr, offset=1) 
    selected_SDP_T <- which(W_SDP_T>= thres_SDP_T)
    
    thres_MVR= knockoff.threshold(W_MVR, fdr=target_fdr, offset=1) 
    selected_MVR <- which(W_MVR >= thres_MVR)
    
    thres_MVR_T <- knockoff.threshold(W_MVR_T, fdr=target_fdr, offset=1) 
    selected_MVR_T <- which(W_MVR_T >= thres_MVR_T)
    
    thres_CI <- knockoff.threshold(W_CI, fdr=target_fdr, offset=1) 
    selected_CI <- which(W_CI>= thres_CI)
    
    thres_CI_T <- knockoff.threshold(W_CI_T, fdr=target_fdr, offset=1) 
    selected_CI_T <- which(W_CI_T>= thres_CI_T)
    
    thres_SDPR <- knockoff.threshold(W_SDPR, fdr=target_fdr, offset=1) 
    selected_SDPR <- which(W_SDPR >= thres_SDPR)
    
    thres_SDPR_T <- knockoff.threshold(W_SDPR_T, fdr=target_fdr, offset=1) 
    selected_SDPR_T <- which(W_SDPR_T >= thres_SDPR_T)
    
    
    mat.features[1, (1+8*(i-1)):(8+8*(i-1))] <- c(fdp(selected_SDP, beta), fdp(selected_SDP_T, beta), fdp(selected_MVR, beta), fdp(selected_MVR_T, beta), fdp(selected_CI, beta), fdp(selected_CI_T, beta), fdp(selected_SDPR, beta), fdp(selected_SDPR_T, beta))
    mat.features[2, (1+8*(i-1)):(8+8*(i-1))] <- c(power(selected_SDP, beta), power(selected_SDP_T, beta), power(selected_MVR, beta), power(selected_MVR_T, beta), power(selected_CI, beta), power(selected_CI_T, beta), power(selected_SDPR, beta), power(selected_SDPR_T, beta))
    mat.features[selected_SDP + 2, 1+8*(i-1)] = 1
    mat.features[selected_SDP_T + 2, 2+8*(i-1)] = 1
    mat.features[selected_MVR + 2, 3+8*(i-1)] = 1
    mat.features[selected_MVR_T + 2, 4+8*(i-1)] = 1
    mat.features[selected_CI + 2, 5+8*(i-1)] = 1
    mat.features[selected_CI_T + 2, 6+8*(i-1)] = 1
    mat.features[selected_SDPR + 2, 7+8*(i-1)] = 1
    mat.features[selected_SDPR_T + 2, 8+8*(i-1)] = 1
    
    colnames(mat.features)[c(1+8*(i-1), 2+8*(i-1), 3+8*(i-1), 4+8*(i-1), 5+8*(i-1), 6+8*(i-1), 7+8*(i-1), 8+8*(i-1))] <- c(str_glue("SDP", target_fdr), str_glue("SDP_T", target_fdr), str_glue("MVR", target_fdr), str_glue("MVR_T", target_fdr),  str_glue("CI", target_fdr), str_glue("CI_T", target_fdr),  str_glue("SDPR", target_fdr), str_glue("SDPR_T", target_fdr))
    
  }
  
  rownames(mat.features) <- c("FDP", "power", c(1:p))
  
  
  return(mat.features)
}

# Comparison between MRC methods: MVR and ME 

Knockoff_comp_KO_MRC <- function(X, scale = FALSE, k, amplitude, list_target_fdr, method){
  
  # X : data matrix
  # scale : if TRUE X is standardized (Default to FALSE)
  # k : number of non-null features
  # amplitude: amplitude of the non null coefficients 
  # list_target_fdr : vector of target FDR 
  # method: either "linear" or "interaction" (to generate y)
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  mat.features <- data.frame(matrix(0, ncol = 2*length(list_target_fdr), nrow = p+2 ))
  
  beta <- beta.sample(p, k, amplitude)
  y <- y.sample(X, beta, method, scale, threshold = 0.5)
  
  knock_sampler <- knockpy$knockoffs$GaussianSampler(X, method = 'mvr')
  X_1 <- knock_sampler$sample_knockoffs()
  X_1 <-  py_to_r(X_1)
  knock_sampler <- knockpy$knockoffs$GaussianSampler(X, method = 'maxent')
  X_2 <- knock_sampler$sample_knockoffs()
  X_2 <- py_to_r(X_2)

  
  W_MVR <- stat.glmnet_coefdiff(X, X_1, y, nfolds=10, family="binomial")
  W_maxent <- stat.glmnet_coefdiff(X, X_2, y, nfolds=10, family="binomial")
  
  
  for (i in 1:length(list_target_fdr)){
    target_fdr <- list_target_fdr[i]
    
    thres_MVR <- knockoff.threshold(W_MVR, fdr=target_fdr, offset=1) 
    selected_MVR <- which(W_MVR>= thres_MVR)
    
    thres_maxent <- knockoff.threshold(W_maxent, fdr=target_fdr, offset=1) 
    selected_maxent <- which(W_maxent>= thres_maxent)
    
    mat.features[1, (1+2*(i-1)):(2+2*(i-1))] <- c(fdp(selected_MVR, beta), fdp(selected_maxent, beta))
    mat.features[2, (1+2*(i-1)):(2+2*(i-1))] <- c(power(selected_MVR, beta), power(selected_maxent, beta))
    mat.features[selected_MVR + 2, 1+2*(i-1)] = 1
    mat.features[selected_maxent + 2, 2+2*(i-1)] = 1

    
    
    colnames(mat.features)[c(1+2*(i-1), 2+2*(i-1)) ] <- c(str_glue("MVR", target_fdr), str_glue("maxent", target_fdr))
    
  }
  
  rownames(mat.features) <- c("FDP", "power", c(1:p))
  
  
  return(mat.features)
}



########################################################################################################################################


                                        ### Functions for graphical representations ###


df_perf_maker <- function(table.perf, list_target_fdr, KO_meth_names, n.rep){ # create data frames for graphical representation
  
  # table.perf: performance tables with FDP and power (rows 1 and 2)
  # list_target_fdr: target fdr used in the experiments to generate table.perf
  # stats_names: names of the different methods used to create table.perf
  # n.rep: number of times the comparison has been done 
  # list_lambda: lambda values used to perform variable selection with the LASSO
  
  n.methods = length(KO_meth_names)
  compteur = 0
  
  len_target_fdr = length(list_target_fdr)
  
  target_fdr = rep(list_target_fdr, n.methods)
  
  df_power_fdr <- data.frame(power = rep(NA, n.methods*len_target_fdr), target_fdr = target_fdr, std = rep(NA, n.methods*len_target_fdr), KO = rep(KO_meth_names, each = len_target_fdr))
  df_power_fdr$KO <- as.vector(df_power_fdr$KO)
  
  df_fdp_fdr <- data.frame(fdp = rep(NA, n.methods*len_target_fdr), target_fdr = target_fdr , std = rep(NA, n.methods*len_target_fdr), KO = rep(KO_meth_names, each = len_target_fdr))
  df_fdp_fdr$KO <- as.vector(df_fdp_fdr$KO)
  
  
  list_power_bin = c()
  list_power_sd_bin = c()
  list_FDP_bin = c()
  list_rep_bin = c()
  
  
  
  for (i in 1:n.methods){
    compteur = compteur+1
    
    
    ### creating df for performances
    
    # fdr curves
    
    table = table.perf
    list_FDP_perf = as.numeric(table[1, seq(i, len_target_fdr*n.methods*n.rep, n.methods)])
    list_power_perf = as.numeric(table[2, seq(i, len_target_fdr*n.methods*n.rep, n.methods)])
    
    
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
    
    
    #fdp curves
    
    list_perf_m = bin(list_FDP_perf, list_power_perf, 0.05)
    
    list_sd = list_perf_m[[3]]
    list_sd[which(is.na(list_sd))]=0
    
    
    
    list_power_bin = c(list_power_bin, list_perf_m[[2]])
    list_power_sd_bin = c(list_power_sd_bin, list_sd)
    list_FDP_bin = c(list_FDP_bin, list_perf_m[[1]])
    list_rep_bin = c(list_rep_bin, length(list_perf_m[[1]]))
    
    
  }
  
  
  df_power_fdp = data.frame(power = list_power_bin, fdp = list_FDP_bin, std = list_power_sd_bin, KO = rep(KO_meth_names, list_rep_bin))
  df_power_fdp$KO <- as.factor(df_power_fdp$KO)
  
  
  
  return(list(df_power_fdr, df_fdp_fdr, df_power_fdp))
}

bin <- function(list_x, list_y, length_bin){ # used for Power vs FDP curves
  
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
  
  for (i in 1:bin_number){
    bound_inf = l_min+(i-1)*length_bin
    bound_sup = l_min+i*length_bin
    
    # Gives the indices of the variables contained in the intervals defined by the size of the bins.
    
    if (bound_sup == l_max){
      l_y = which(list_x <= bound_sup & list_x >= bound_inf) # accounting for the upper bound of the list
    }
    else{l_y = which(list_x < bound_sup & list_x >= bound_inf)}
    
    # Get the mean 
    list_y_m = c(list_y_m, mean(list_y[l_y]))
    
    # Get the standard deviation
    
    if (length(list_y[l_y])==1){ # when there is only one value in the bin
      list_y_sd = c(list_y_sd, 0)
    }
    else{
      list_y_sd = c(list_y_sd, sd(list_y[l_y]))
    }
    
    # determined the value on the x-axis (middle of the bin)
    
    if (bound_sup > l_max){bound_sup = l_max}# accounting for last bin when the size of the bin might be reduced
    
    list_x_m = c(list_x_m, bound_inf+(bound_sup-bound_inf)/2)
  }
  
  i = intersect(which(is.na(list_y_m)),which(is.na(list_y_sd))) #remove empty bins
  if (length(intersect(which(is.na(list_y_m)),which(is.na(list_y_sd))))!=0){
    list_x_m = list_x_m[-i]
    list_y_m = list_y_m[-i]
    list_y_sd = list_y_sd[-i] 
  }
  return(list(list_x_m, list_y_m, list_y_sd))
}
