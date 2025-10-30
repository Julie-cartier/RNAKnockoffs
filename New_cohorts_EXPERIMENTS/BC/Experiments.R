# get access to the settings
source("/home/julie/Documents/Paper_codes/New_cohorts_EXPERIMENTS/BC/Settings.R")

# get access to the functions
source("/home/julie/Documents/Paper_codes/New_cohorts_EXPERIMENTS/Functions.R")
source('/home/julie/Documents/Paper_codes/New_cohorts_EXPERIMENTS/BC/Add_functions.R')

#get access to python packages 

require(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/home/julie/anaconda3/envs/py-environment/bin/python")
knockpy <- import("knockpy", convert = FALSE)
sklearn <- import("sklearn", convert = FALSE)
numpy <- import("numpy", convert = FALSE)


############################ X_BC (experiments), linear y,  statistics comparison ####################################################
######################################################################################################################################

list.beta <- list()
table.perf.linear.comp <- data.frame(rep(NA, ncol(X_BC)+2))

#  Perform the comparison a 100 times 
tic()
for (i in 1:100){
  print(i)
  table.perf.part <- Knockoff_comp_stat_LSCIP(scale(X_BC), scale = FALSE, k = 10, amplitude = 10, list_target_fdr = list_target_fdr, list_lambda = list_lambda, method = "linear")
  table.perf.linear.comp <- cbind(table.perf.linear.comp, table.perf.part[[1]])
  list.beta[[i]] <- table.perf.part[[2]]
}

toc()

table.perf.linear.comp <- table.perf.linear.comp[, -1]

# save the table in an R file to use later to make plots

save(table.perf.linear.comp, list.beta, file = "/home/julie/Documents/Paper_codes/New_cohorts_EXPERIMENTS/BC/R_files/table_perf_statistic_comparaison_3e_cohorte.R")

######################################################################################################################################
######################################################################################################################################



############################ X_BC (experiments), linear y,  methods comparison #######################################################
######################################################################################################################################

table.perf.linear.comp.meth <- data.frame(rep(NA, ncol(X_BC)+2))
list.beta = list()


#  Perform the comparison a 100 times 

tic()

for (i in 1:100){
  print(i)
  table.perf.part <- Knockoff_comp_tot(scale(X_BC), scale = FALSE, k = 10, amplitude = 10, list_target_fdr = list_target_fdr, method = "linear") 
  table.perf.linear.comp.meth <- cbind(table.perf.linear.comp.meth, table.perf.part[[1]])
  list.beta[[i]] <- table.perf.part[[2]]
  
}

toc()

table.perf.linear.comp.meth <- table.perf.linear.comp.meth[, -1]

# save the table in an R file to use later to make plots

save(table.perf.linear.comp.meth, list.beta, file = "/home/julie/Documents/Paper_codes/New_cohorts_EXPERIMENTS/BC/R_files/table_perf_method_comparaison_3e_cohorte.R")

######################################################################################################################################
######################################################################################################################################

######################################################################################################################################
######################################################################################################################################

######################################################################################################################################
######################################################################################################################################


############################### KO, wilcoxon and Lasso comparison ####################################################################
######################################################################################################################################


load("/home/julie/Documents/Paper_codes/New_cohorts_EXPERIMENTS/BC/R_files/table_perf_method_comparaison_3e_cohorte.R")

table.perf.linear.KO <- table.perf.linear.comp.meth[, grepl("LSCIP", names(table.perf.linear.comp.meth))]
table.perf.linear.LPLR <- data.frame(rep(NA, ncol(X_BC)+2))
table.perf.linear.wilcoxon <- data.frame(rep(NA, ncol(X_BC)+2))
table.perf.linear.LPLR.lambda.min <- data.frame(matrix(ncol = 100, nrow = 3))
table.perf.linear.LPLR.lambda.oracle <- data.frame(matrix(ncol = 100, nrow = 3))

tic()

for (i in 1:100){
  print(i)
  beta <- list.beta[[i]]
  y <- y.sample(scale(X_BC), beta, method = "linear", scale = FALSE, threshold = 0.5)
  
  table.perf.part <- LPLR_min_oracle(scale(X_BC, center = TRUE, scale = TRUE), y = y, list_lambda = list_lambda)
  table.perf.linear.LPLR <- cbind(table.perf.linear.LPLR, table.perf.part[[1]])
  table.perf.linear.LPLR.lambda.min[,i] <- table.perf.part[[2]]
  table.perf.linear.LPLR.lambda.oracle[,i] <- table.perf.part[[3]]
  
  table.perf.part <- DE_Wilcoxon_BH(scale(X_BC), y = y, beta = beta, list_target_fdr = list_target_fdr_wilcoxon, method = "linear")
  table.perf.linear.wilcoxon <- cbind(table.perf.linear.wilcoxon, table.perf.part)
  
}

toc()

table.perf.linear.LPLR<- table.perf.linear.LPLR[, -1]
table.perf.linear.wilcoxon <- table.perf.linear.wilcoxon[, -1]
rownames(table.perf.linear.LPLR.lambda.min) <- c("lambda.min", "FDP", "Power")
rownames(table.perf.linear.LPLR.lambda.oracle) <- c("lambda.oracle", "FDP", "Power")


save(table.perf.linear.KO, table.perf.linear.LPLR, table.perf.linear.wilcoxon, table.perf.linear.LPLR.lambda.min, table.perf.linear.LPLR.lambda.oracle, file = "/home/julie/Documents/Paper_codes/New_cohorts_EXPERIMENTS/BC/R_files/LPLR_KO_BH_linear_table_perfx10.R")

