source('/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Wilcoxon_vs_KO_vs_LASSO_vs_KOPI/Settings.R')
source('/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Wilcoxon_vs_KO_vs_LASSO_vs_KOPI/Functions.R')


# Set the Python environment
Sys.setenv(RETICULATE_PYTHON = "/home/julie/anaconda3/envs/py-environment/bin/python")

# Import necessary Python modules
knockpy <- import("knockpy", convert = FALSE)
sklearn <- import("sklearn", convert = FALSE)
numpy <- import("numpy", convert = FALSE)

# Source the home-made Python file
source_python("/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Wilcoxon_vs_KO_vs_LASSO_vs_KOPI/Functions.py")




################################################################################################################
############################# True X (CRUKPAP), simulated linear y ####################################################


# pre-computed bank of 100 KO matrices obtained with the LSCIP method and X_real
# using the KO_banks_cluster_script script (with less iterations)

load("/home/julie/Documents/Paper_codes/Data/CRUKPAP/KO_partial_bank_100_CRUKPAP_LSCIP/KO_real_LSCIP_bank_part.R")
#load("/media/julie/T5 EVO/Copie ordinateur/Paper_codes/Data/CRUKPAP/KO_partial_bank_100_CRUKPAP_LSCIP/KO_real_LSCIP_bank_part.R")

list.beta <- list()
table.perf.linear.LPLRKO <- data.frame(rep(NA, ncol(X_real)+2))
table.perf.linear.wilcoxon <- data.frame(rep(NA, ncol(X_real)+2))
table.perf.linear.LPLR.lambda.min <- data.frame(matrix(ncol = 100, nrow = 3))

tic()

for (i in 1:100){
  
  print(i)
  X_k <- list_KO_real_LSCIP[[i]]
  
  beta <- beta.sample(ncol(X_real), k = 10, amplitude = 10)
  y <- y.sample(scale(X_real), beta, method = "linear", scale = FALSE, threshold = 0.5)
  
  table.perf.part <- Knockoff_LPLR(scale(X_real, center = TRUE, scale = TRUE), X_k = X_k, y = y, beta = beta, list_target_fdr = list_target_fdr, list_lambda = list_lambda, method = "linear")
  table.perf.linear.LPLRKO <- cbind(table.perf.linear.LPLRKO, table.perf.part[[1]])
  
  table.perf.linear.LPLR.lambda.min[,i] <- table.perf.part[[2]]
  
  table.perf.part <- DE_Wilcoxon_BH(scale(X_real), y = y, beta = beta, list_target_fdr = list_target_fdr_wilcoxon, method = "linear")
  table.perf.linear.wilcoxon <- cbind(table.perf.linear.wilcoxon, table.perf.part)
  
  list.beta[[i]] <- beta
}

toc()

table.perf.linear.LPLRKO<- table.perf.linear.LPLRKO[, -1]
table.perf.linear.wilcoxon <- table.perf.linear.wilcoxon[, -1]
rownames(table.perf.linear.LPLR.lambda.min) <- c("lambda.min", "FDP", "Power")


save(list.beta, table.perf.linear.LPLRKO, table.perf.linear.wilcoxon, table.perf.linear.LPLR.lambda.min, file = "/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Wilcoxon_vs_KO_vs_LASSO_vs_KOPI/R_files/LPLR_KO_BH_linear_table_perfx10.R")

###################### Additional comparison (KOPI and oracle LASSO) ###########################################

# load data to get the list of beta for further comparison

load("/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Wilcoxon_vs_KO_vs_LASSO_vs_KOPI/R_files/LPLR_KO_BH_linear_table_perfx10.R")

## Get the results with KOPI and The oracle lambda ##

table.perf.linear.KOPI <- data.frame(rep(NA, ncol(X_real)+2))
table.perf.linear.LPLR.lambda.oracle <- data.frame(matrix(ncol = 100, nrow = 3))


tic()
for (i in 1:100){
  # pre-computed bank of 10 000 KO matrices obtained with the LSCIP method and X_real (KO_banks_cluster_script)
  load(str_glue("/media/julie/T5 EVO/Data_cluster/KO_bank_under_sets/X_k_",i,"_100.R"))
  
  print(i)
  
  beta <- list.beta[[i]]
  y <- y.sample(scale(X_real), beta, method = "linear", scale = FALSE, threshold = 0.5)
  
  table.perf.part <- Add_KOPI_LPLR_oracle(scale(X_real), list_X_k_100, y = y, beta = beta, list_target_FDP = list_target_fdr, k_max = round(p/50))
  table.perf.linear.KOPI <- cbind(table.perf.linear.KOPI, table.perf.part[[1]])
  
  table.perf.linear.LPLR.lambda.oracle[,i] <- table.perf.part[[2]]
}

toc()

table.perf.linear.KOPI <- table.perf.linear.KOPI[, -1]
rownames(table.perf.linear.LPLR.lambda.oracle) <- c("lambda.oracle", "FDP", "Power")


save(table.perf.linear.KOPI, table.perf.linear.LPLR.lambda.oracle, file = "/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Wilcoxon_vs_KO_vs_LASSO_vs_KOPI/R_files/KOPI_KO_BH_add_linear_table_perfx10.R")




################################################################################################################
############################# For KOPI vs KO comparison, True X (CRUKAP), linear y  ############################
################################################################################################################


list.beta <- list()
table.perf.linear.KOPIKO <- data.frame(rep(NA, ncol(X_real)+2))
table.perf.linear.wilcoxon <- data.frame(rep(NA, ncol(X_real)+2))


tic()
for (i in 1:100){
  # pre-computed bank of 10 000 KO matrices obtained with the LSCIP method and X_real
  load(str_glue("/media/julie/T5 EVO/Data_cluster/KO_bank_under_sets/X_k_",i,"_100.R"))
  
  print(i)
  #D_draws_X_k <- list_KO_real_LSCIP[(Draws*(i-1)+1):(Draws*i)]
  
  beta <- beta.sample(ncol(X_real), k = 10, amplitude = 10)
  y <- y.sample(scale(X_real), beta, method = "linear", scale = FALSE, threshold = 0.5)
  
  table.perf.part <- KOPI_KO_comp(scale(X_real), list_X_k_100, y = y, beta = beta, list_target_fdr_FDP = list_target_fdr, alpha = 0.1, B = 1000, k_max = round(p/50))
  table.perf.linear.KOPIKO <- cbind(table.perf.linear.KOPIKO, table.perf.part)

  
  list.beta[[i]] <- beta
}

toc()

table.perf.linear.KOPIKO <- table.perf.linear.KOPIKO[, -1]


save(list.beta, table.perf.linear.KOPIKO, file = "/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Wilcoxon_vs_KO_vs_LASSO_vs_KOPI/R_files/KOPI_KO_linear_table_perfx10.R")


################################################################################################################
############################# True X (CRUKPAP), simulated y with interactions ##################################


# pre-computed bank of 100 KO matrices obtained with the LSCIP method and X_real
# using the KO_banks_cluster_script script (with less iterations)

#load("/home/julie/Documents/Paper_codes/Data/CRUKPAP/KO_partial_bank_100_CRUKPAP_LSCIP/KO_real_LSCIP_bank_part.R")
load("/media/julie/T5 EVO/Copie ordinateur/Paper_codes/Data/CRUKPAP/KO_partial_bank_100_CRUKPAP_LSCIP/KO_real_LSCIP_bank_part.R")

load("/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Statistics_comparison/R_files/list_beta_interaction.R")

list.beta <- list()
table.perf.interaction.LPLRKO <- data.frame(rep(NA, ncol(X_real)+2))
table.perf.interaction.wilcoxon <- data.frame(rep(NA, ncol(X_real)+2))
table.perf.interaction.LPLR.lambda.min <- data.frame(matrix(ncol = 100, nrow = 3))

tic()

for (i in 1:10){
  
  print(i)
  X_k <- list_KO_real_LSCIP[[i]]
  
  y <- list_y_interaction[[i]]
  beta <- list_beta_interaction[[i]]
  
  table.perf.part <- Knockoff_LPLR(scale(X_real, center = TRUE, scale = TRUE), X_k = X_k, y = y, beta = beta, list_target_fdr = list_target_fdr, list_lambda = list_lambda, method = "interaction")
  table.perf.interaction.LPLRKO <- cbind(table.perf.interaction.LPLRKO, table.perf.part[[1]])
  
  table.perf.interaction.LPLR.lambda.min[,i] <- table.perf.part[[2]]
  
  table.perf.part <- DE_Wilcoxon_BH(scale(X_real), y = y, beta = beta, list_target_fdr = list_target_fdr_wilcoxon, method = "interaction")
  table.perf.interaction.wilcoxon <- cbind(table.perf.interaction.wilcoxon, table.perf.part)
  
  list.beta[[i]] <- beta
}

toc()

table.perf.interaction.LPLRKO<- table.perf.interaction.LPLRKO[, -1]
table.perf.interaction.wilcoxon <- table.perf.interaction.wilcoxon[, -1]
rownames(table.perf.interaction.LPLR.lambda.min) <- c("lambda.min", "FDP", "Power")


save(list.beta, table.perf.interaction.LPLRKO, table.perf.interaction.wilcoxon, table.perf.interaction.LPLR.lambda.min, file = "/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Wilcoxon_vs_KO_vs_LASSO_vs_KOPI/R_files/LPLR_KO_BH_interaction_table_perfx10.R")


