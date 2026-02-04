# This script is used to study feature selection stability using the Stability Selection methods
# on the CRUKPAP (with only 749 genes)

library(sharp)
library(stabs)

#source("Setup/settings.R")
#source("Setup/functions.R")

source("/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Stability/Settings.R")
source("/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Stability/Functions.R")

list_beta <- list()
for (rep in 1:10){
  load(str_glue("/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Stability/R_files/y_rep/table.selected_genes_cv_KOPI_",rep,".R"))
  list_beta[[rep]]<- beta
}


# Load the CRUKPAP dataset (with only the identified 749 risk genes)
#X_real <- read.table("/Users/youmnaayadi/Documents/PhD_RNA_knockoffs/Data/CRUKPAP/CRUKPAP_risk_genes_only.txt", header = TRUE)
X_scaled <- scale(X_real)

n <- nrow(X_scaled)
p <- ncol(X_scaled)

k <- 10 # Number of positive / non-null features

y_rep <- 10 # Repeat the experiment 10 times with different outcomes y 
n_cv <- 10
n_folds <- 10

  
for (rep in 1:y_rep){
  print(paste0("Rep ", rep))
  
  beta <- list_beta[[rep]]
  y <- y.sample(X_scaled, beta, scale = FALSE, method = "linear")
    
    # Initialize list of empty data frames that will store the selected genes for the different methods
    selected_genes.Consensus <- data.frame(matrix(0, nrow = p, ncol = (n_cv * n_folds)))
    selected_genes.MB <- data.frame(matrix(0, nrow = p, ncol = (n_cv * n_folds)))
    selected_genes.SS <- data.frame(matrix(0, nrow = p, ncol = (n_cv * n_folds)))
    
    fcount <- 1 # fold count
    
    for (i in 1:n_cv){
      print(paste0("---- CV ", i))
      
      vec <- sample(1:n, replace = FALSE)
      split_vec <- split(vec, cut(seq_along(vec), n_folds, labels = FALSE)) # split the vector for one cross validation
      
      for (j in 1:n_folds){
        print(paste0("-------- Fold ", j))
        fold_index <- split_vec[[j]] # Indices of the left-out fold
        
        # Real X subsampled and scaled
        X_cv <- X_scaled[-fold_index, ]
        
        # Feature selection
        ## "Consensus" method
        stab <- VariableSelection(xdata = X_cv, ydata = y[-fold_index], family = "binomial", K = 100, seed = sample(1:(y_rep*n_cv*n_folds), 1), verbose = FALSE)
        selected_genes.Consensus[, fcount] <- SelectedVariables(stab)
        
        ## Meinshausen & Buehlmann (MB) method
        stab.MB <- stabsel(x = X_cv, y = y[-fold_index], fitfun = glmnet.lasso, cutoff = 0.75, sampling.type = "MB",PFER = 1)
        selected_genes.MB[stab.MB$selected, fcount] <- 1 
        
        
        ## SS method
        stab.SS <- stabsel(x = X_cv, y = y[-fold_index], fitfun = glmnet.lasso, cutoff = 0.75, sampling.type = "SS", PFER = 1)
        selected_genes.SS[stab.SS$selected, fcount] <- 1 
        
        fcount <- fcount + 1
      }
    }
    
    # Storing the selected genes for a particular p/n ratio and value of y into .Rdata files
    save(beta, y, selected_genes.Consensus, selected_genes.MB, selected_genes.SS, file = paste0("/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Stability/R_files/y_rep/SS_selected_features_rep_", rep, ".Rdata"))
    
}



