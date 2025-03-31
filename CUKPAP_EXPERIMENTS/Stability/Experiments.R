source("/home/julie/Documents/Paper_codes/CUKPAP_EXPERIMENTS/Stability/Settings.R")
source("/home/julie/Documents/Paper_codes/CUKPAP_EXPERIMENTS/Stability/Functions.R")

library(reticulate)

# Set the Python environment
Sys.setenv(RETICULATE_PYTHON = "/home/julie/anaconda3/envs/py-environment/bin/python")

# Import necessary Python modules
knockpy <- import("knockpy", convert = FALSE)
sklearn <- import("sklearn", convert = FALSE)
numpy <- import("numpy", convert = FALSE)

# Source the home-made Python file to perform KOPI
source_python("/home/julie/Documents/Paper_codes/CUKPAP_EXPERIMENTS/Wilcoxon_vs_KO_vs_LASSO_vs_KOPI/Functions.py")


### !!! DO not forget to rerun the List_indices_maker.R script if the script in the cluster has been launched again !!! 

# Get CV indices for the comparison

load("/home/julie/Documents/Paper_codes/CUKPAP_EXPERIMENTS/Stability/R_files/list_cv_indices.R")



################################################################################################################################################
###################################### True X, fake Y, 10 simulated y ##########################################################################
######################################## Stability 10 x 10-fold CV #############################################################################
################################# Comparison : LPLR (CV, orcale), KOPI, KO #####################################################################
################################################################################################################################################
# Repeat : (10 times 10 - fold cross validation) x (iteration = 10)
# For the ten iterations, only beta and y change. The subset of the 10 times 10 - fold cross validation remain the same, so does the corresponding KO matrices

y.rep = 10

for (rep in 1:y.rep){# Repeat the experiment 10 times with different outcome y 
  
  beta <- beta.sample(p, 10, amplitude = 10)
  y <- y.sample(scale(X_real), beta, method = "linear", scale = TRUE)
  X_scaled <- scale(X_real)
  
  
  # LASSO Penalized Logistic Regression (LPLR), cv lambda 
  
  tic()
  
  cl <- parallel::makeCluster(5)
  doParallel::registerDoParallel(cl)
  table.selected_genes.perf.LPLR <- foreach(j = 1:100, .combine = 'cbind', .packages = c("glmnet", "stringr","PRROC")) %dopar% {LPLR_stability(scale(X_real), y, list.cv.ind[[j]])}
  parallel::stopCluster(cl)
  
  toc()
  
  table.aucs.LPLR <- table.selected_genes.perf.LPLR[1:2,]
  table.selected.genes.LPLR <- table.selected_genes.perf.LPLR[3:751,]
  rownames(table.selected.genes.LPLR) <- colnames(X_real)
  
  save(y, beta, table.selected.genes.LPLR, table.aucs.LPLR, file = str_glue("/home/julie/Documents/Paper_codes/CUKPAP_EXPERIMENTS/Stability/R_files/y_rep/table_selected_genes_LPLR_",rep,".R"))
  
  
  # LPLR, oracle lambda 
  
  tic()
  
  cl <- parallel::makeCluster(5)
  doParallel::registerDoParallel(cl)
  table.selected_genes.perf.LPLR.oracle <- foreach(j = 1:100, .combine = 'cbind', .packages = c("glmnet", "stringr","PRROC")) %dopar% {LPLR_stability_oracle(scale(X_real), y, list.cv.ind[[j]], beta)}
  parallel::stopCluster(cl)
  
  toc()
  
  table.aucs.LPLR.oracle <- table.selected_genes.perf.LPLR.oracle[1:4,]
  table.selected.genes.LPLR.oracle <- table.selected_genes.perf.LPLR.oracle[5:753,]
  rownames(table.selected.genes.LPLR.oracle) <- c(colnames(X_real))
  
  save(y, beta, table.selected.genes.LPLR.oracle, table.aucs.LPLR.oracle, file = str_glue("/home/julie/Documents/Paper_codes/CUKPAP_EXPERIMENTS/Stability/R_files/y_rep/table_selected_genes_LPLR_oracle_",rep,".R"))
  
  
  # KO + KOPI
  
  table.selected_genes.Vanilla_0.1 <- data.frame(rep(NA, ncol(X_real)))
  table.selected_genes.Vanilla_0.2 <- data.frame(rep(NA, ncol(X_real)))
  table.selected_genes.Vanilla_0.3 <- data.frame(rep(NA, ncol(X_real)))
  table.selected_genes.Vanilla_0.5 <- data.frame(rep(NA, ncol(X_real)))
  
  
  table.selected_genes.KOPI_0.1 <- data.frame(rep(NA, ncol(X_real)))
  table.selected_genes.KOPI_0.2 <- data.frame(rep(NA, ncol(X_real)))
  table.selected_genes.KOPI_0.3 <- data.frame(rep(NA, ncol(X_real)))
  table.selected_genes.KOPI_0.5 <- data.frame(rep(NA, ncol(X_real)))
  
  tic()
  
  for (CV in 1:10){

    load(str_glue('/media/julie/T5 EVO/Data_cluster/KO_bank_CV_sets/CV_indices/CV_',CV,'.R')) # get indices
    
    
    for (fold in 1:10){
      ind.vec <- split_vec[[fold]]
      load(str_glue('/media/julie/T5 EVO/Data_cluster/KO_bank_CV_sets/CV_', CV, '_fold_', fold,  '_X_k_100.R')) # get pre-computed KO matrices (KO_banks_CV_cluster_script.R)
      
      
      print(str_glue('CV : ',CV, ', fold : ',fold))
      
      table.selected_genes <- KOPI_KO_stability(X_scaled[-ind.vec,], y[-ind.vec], list_KO_cv_100)
      
      table.selected_genes.Vanilla <- table.selected_genes[[1]]
      table.selected_genes.KOPI <- table.selected_genes[[2]]
      
      table.selected_genes.Vanilla_0.1 <- cbind(table.selected_genes.Vanilla_0.1, table.selected_genes.Vanilla[,1])
      table.selected_genes.Vanilla_0.2 <- cbind(table.selected_genes.Vanilla_0.2, table.selected_genes.Vanilla[,2])
      table.selected_genes.Vanilla_0.3 <- cbind(table.selected_genes.Vanilla_0.3, table.selected_genes.Vanilla[,3])
      table.selected_genes.Vanilla_0.5 <- cbind(table.selected_genes.Vanilla_0.5, table.selected_genes.Vanilla[,4])
      
      table.selected_genes.KOPI_0.1 <- cbind(table.selected_genes.KOPI_0.1, table.selected_genes.KOPI[,1])
      table.selected_genes.KOPI_0.2 <- cbind(table.selected_genes.KOPI_0.2, table.selected_genes.KOPI[,2])
      table.selected_genes.KOPI_0.3 <- cbind(table.selected_genes.KOPI_0.3, table.selected_genes.KOPI[,3])
      table.selected_genes.KOPI_0.5 <- cbind(table.selected_genes.KOPI_0.5, table.selected_genes.KOPI[,4])
    }
  }
  
  toc()
  
  table.selected_genes.KOPI_0.1 <- table.selected_genes.KOPI_0.1[,-1]
  table.selected_genes.KOPI_0.2 <- table.selected_genes.KOPI_0.2[,-1]
  table.selected_genes.KOPI_0.3 <- table.selected_genes.KOPI_0.3[,-1]
  table.selected_genes.KOPI_0.5 <- table.selected_genes.KOPI_0.5[,-1]
  
  table.selected_genes.Vanilla_0.1 <- table.selected_genes.Vanilla_0.1[,-1]
  table.selected_genes.Vanilla_0.2 <- table.selected_genes.Vanilla_0.2[,-1]
  table.selected_genes.Vanilla_0.3 <- table.selected_genes.Vanilla_0.3[,-1]
  table.selected_genes.Vanilla_0.5 <- table.selected_genes.Vanilla_0.5[,-1]

  save(beta, y, table.selected_genes.Vanilla_0.1, table.selected_genes.Vanilla_0.2, table.selected_genes.Vanilla_0.3, table.selected_genes.Vanilla_0.5, file = str_glue("/home/julie/Documents/Paper_codes/CUKPAP_EXPERIMENTS/Stability/R_files/y_rep/table.selected_genes_cv_Vanilla_",rep,".R"))
  save(beta, y, table.selected_genes.KOPI_0.1, table.selected_genes.KOPI_0.2, table.selected_genes.KOPI_0.3, table.selected_genes.KOPI_0.5, file = str_glue("/home/julie/Documents/Paper_codes/CUKPAP_EXPERIMENTS/Stability/R_files/y_rep/table.selected_genes_cv_KOPI_",rep,".R"))
  
}



################################################################################################################################################
###################################### True X, fake Y, KOPI vs KO  #############################################################################
######################## Stability over 100 iterations (study KOPI impact on stability) ########################################################
################################# Comparison :   KOPI, KO ######################################################################################
################################################################################################################################################
# Repeat : 100 times with the same y
# At each iteration : 


beta = beta.sample(p, 10, amplitude = 10)
y = y.sample(X_real, beta, method = "linear", scale = TRUE)

table.selected_genes.Vanilla_0.1 <- data.frame(rep(NA, ncol(X_real)))
table.selected_genes.Vanilla_0.2 <- data.frame(rep(NA, ncol(X_real)))
table.selected_genes.Vanilla_0.3 <- data.frame(rep(NA, ncol(X_real)))
table.selected_genes.Vanilla_0.5 <- data.frame(rep(NA, ncol(X_real)))


table.selected_genes.KOPI_0.1 <- data.frame(rep(NA, ncol(X_real)))
table.selected_genes.KOPI_0.2 <- data.frame(rep(NA, ncol(X_real)))
table.selected_genes.KOPI_0.3 <- data.frame(rep(NA, ncol(X_real)))
table.selected_genes.KOPI_0.5 <- data.frame(rep(NA, ncol(X_real)))

for (i in 1:100){
  
  # load KO matrices from the bank of KO matrices generated with the KO_banks_cluster_script.R (/home/julie/Documents/Paper_codes/CUKPAP_EXPERIMENTS/Wilcoxon_vs_KO_vs_LASSO_vs_KOPI/)
  load(str_glue("/media/julie/T5 EVO/Data_cluster/KO_bank_under_sets/X_k_",i,"_100.R"))
  
  print(i)
  
  table.selected_genes <- KOPI_KO_stability(scale(X_real), y, list_X_k_100)
  
  table.selected_genes.Vanilla <- table.selected_genes[[1]]
  table.selected_genes.KOPI <- table.selected_genes[[2]]
  
  table.selected_genes.Vanilla_0.1 <- cbind(table.selected_genes.Vanilla_0.1, table.selected_genes.Vanilla[,1])
  table.selected_genes.Vanilla_0.2 <- cbind(table.selected_genes.Vanilla_0.2, table.selected_genes.Vanilla[,2])
  table.selected_genes.Vanilla_0.3 <- cbind(table.selected_genes.Vanilla_0.3, table.selected_genes.Vanilla[,3])
  table.selected_genes.Vanilla_0.5 <- cbind(table.selected_genes.Vanilla_0.5, table.selected_genes.Vanilla[,4])
  
  table.selected_genes.KOPI_0.1 <- cbind(table.selected_genes.KOPI_0.1, table.selected_genes.KOPI[,1])
  table.selected_genes.KOPI_0.2 <- cbind(table.selected_genes.KOPI_0.2, table.selected_genes.KOPI[,2])
  table.selected_genes.KOPI_0.3 <- cbind(table.selected_genes.KOPI_0.3, table.selected_genes.KOPI[,3])
  table.selected_genes.KOPI_0.5 <- cbind(table.selected_genes.KOPI_0.5, table.selected_genes.KOPI[,4])
}


table.selected_genes.KOPI_0.1 <- table.selected_genes.KOPI_0.1[,-1]
table.selected_genes.KOPI_0.2 <- table.selected_genes.KOPI_0.2[,-1]
table.selected_genes.KOPI_0.3 <- table.selected_genes.KOPI_0.3[,-1]
table.selected_genes.KOPI_0.5 <- table.selected_genes.KOPI_0.5[,-1]

table.selected_genes.Vanilla_0.1 <- table.selected_genes.Vanilla_0.1[,-1]
table.selected_genes.Vanilla_0.2 <- table.selected_genes.Vanilla_0.2[,-1]
table.selected_genes.Vanilla_0.3<- table.selected_genes.Vanilla_0.3[,-1]
table.selected_genes.Vanilla_0.5 <- table.selected_genes.Vanilla_0.5[,-1]


save(beta, y, table.selected_genes.Vanilla_0.1, table.selected_genes.Vanilla_0.2, table.selected_genes.Vanilla_0.3, table.selected_genes.Vanilla_0.5, file = "/home/julie/Documents/Paper_codes/CUKPAP_EXPERIMENTS/Stability/R_files/table.selected_genes_no_cv_Vanilla.R")
save(beta, y, table.selected_genes.KOPI_0.1, table.selected_genes.KOPI_0.2, table.selected_genes.KOPI_0.3, table.selected_genes.KOPI_0.5, file = "/home/julie/Documents/Paper_codes/CUKPAP_EXPERIMENTS/Stability/R_files/table.selected_genes_no_cv_KOPI.R")


