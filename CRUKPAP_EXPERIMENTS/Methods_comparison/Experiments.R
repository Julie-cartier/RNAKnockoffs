source("/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Methods_comparison/Settings.R")
source("/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Methods_comparison/Functions.R")

#use python functions 
require(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/home/julie/anaconda3/envs/py-environment/bin/python")
knockpy <- import("knockpy", convert = FALSE)
sklearn <- import("sklearn", convert = FALSE)
numpy <- import("numpy", convert = FALSE)



############################ TRUE X, linear Y,  LSCIP VS MVR vs SDP vs CI#############################################################
######################################################################################################################################

table.perf.linear.comp_tot <- data.frame(rep(NA, ncol(X_real)+2))
list_beta <- list()


tic()
for (i in 1:100){
  print(i)
  table.perf.part <- Knockoff_comp_tot(scale(X_real), scale = FALSE, k = 10, amplitude = 10, list_target_fdr, method = "linear")
  table.perf.linear.comp_tot <- cbind(table.perf.linear.comp_tot, table.perf.part[[1]])
  list_beta[[i]] <- table.perf.part[[2]] # get beta vector
}
toc()

table.perf.linear.comp_tot <- table.perf.linear.comp_tot[, -1]

save(list_beta, table.perf.linear.comp_tot, file = "/home/julie/Bureau/Methods/Knockoffs/Knockoff_generation_synthesis/R_files/Comp_tot/table_perfx10.R")


######################################################################################################################################
######################################################################################################################################




############################ TRUE X, linear y,  MVR vs ME KO #########################################################################
######################################################################################################################################

table.perf.linear.MRC <- data.frame(rep(NA, ncol(X_real)+2))


tic()
for (i in 1:10){
  print(i)
  table.perf.part <- Knockoff_comp_KO_MRC(scale(X_real), scale = FALSE, k = 10, amplitude = 10, list_target_fdr, method="linear")
  table.perf.linear.MRC <- cbind(table.perf.linear.MRC, table.perf.part)
}
toc()

table.perf.linear.MRC <- table.perf.linear.MRC[, -1]

save(table.perf.linear.MRC, file = "/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Methods_comparison/R_files/MRC_table_perfx10.R")

######################################################################################################################################
######################################################################################################################################




####################### Simulated X (known covariance structure), linear y, CI/MVR/SDP ###############################################
######################################################################################################################################

table.perf.linear.cov <- data.frame(rep(NA, ncol(X_real)+2))


tic()
for (i in 1:10){
  print(i)
  table.perf.part <- Knockoff_comp_KO_COV(X.sim2, Sigma.sim = S.sim2, scale = FALSE, k = 10, amplitude = 10, list_target_fdr, method = "linear")
  table.perf.linear.cov <- cbind(table.perf.linear.cov, table.perf.part)
}
toc()

table.perf.linear.cov <- table.perf.linear.cov[, -1]

save(table.perf.linear.cov, file = "/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Methods_comparison/R_files/table_perf_extx10.R")


######################################################################################################################################
######################################################################################################################################










