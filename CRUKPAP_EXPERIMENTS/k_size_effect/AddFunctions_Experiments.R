#get settings and usual functions

source("/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Statistics_comparison/Settings.R")
source("/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Statistics_comparison/Functions.R")

# Additional function

KO_ksize_effect <- function(X, X_k = NULL, vec_k, amplitude, list_target_fdr,  method){
  
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
  
  if (is.null(X_k)){
    X_k <- create.LSCIP_parallel(X)
  }
  
  compteur = 0
  for (k in vec_k){
    
    compteur = compteur + 1
    
    #simulate y 
    beta <- beta.sample(p, k, amplitude)
    y <- y.sample(X, beta, method, scale = FALSE, threshold = 0.5)
    
    # compute the LCD statistic
    W_LCD <- stat.glmnet_coefdiff(X, X_k, y, nfolds=10, family="binomial")
    
    # get sets of selected features for different targets FDR level
    for (i in 1:length(list_target_fdr)){
      target_fdr <- list_target_fdr[i]
      
      thres_LCD <- knockoff.threshold(W_LCD, fdr=target_fdr, offset=1) 
      selected_LCD <- which(W_LCD >= thres_LCD)
      
      
      mat.features[1, compteur + k_length*(i-1)] <- fdp(selected_LCD, beta)
      mat.features[2, compteur + k_length*(i-1)] <- power(selected_LCD, beta)
      mat.features[selected_LCD + 2, compteur + k_length*(i-1)] <- 1
      
      colnames(mat.features)[compteur + k_length*(i-1)] <- c(str_glue("k_", k, "_tFDR_", i))
      
    }
  }
  return(mat.features)
}

# Experiments 


table.perf.linear.ksize <- data.frame(rep(NA, ncol(X_real)+2))



tic()
for (i in 1:1){
  print(i)
  table.perf.part <- KO_ksize_effect(scale(X_real), vec_k = c(10,30,50), amplitude = 10, list_target_fdr = list_target_fdr,  method = "linear")  
  table.perf.linear.ksize <- cbind(table.perf.linear.ksize, table.perf.part)
}

toc()

table.perf.linear.ksize <- table.perf.linear.ksize[, -1]


save(table.perf.linear.ksize, file = "/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/k_size_effect/R_files/ksize_trueX_fakeY_tables_perfx10.R")
