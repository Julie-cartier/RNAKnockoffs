# get access to the settings
source("/home/julie/Documents/Paper_codes/New_cohorts_EXPERIMENTS/AEGIS/Settings.R")

# get access to the functions
source("/home/julie/Documents/Paper_codes/New_cohorts_EXPERIMENTS/Functions.R")

#get access to python packages 

require(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/home/julie/anaconda3/envs/py-environment/bin/python")
knockpy <- import("knockpy", convert = FALSE)
sklearn <- import("sklearn", convert = FALSE)
numpy <- import("numpy", convert = FALSE)


############################ X_AEGIS (experiments), linear y,  statistics comparison #################################################
######################################################################################################################################

# Linear 
list.beta <- list()
table.perf.linear.comp <- data.frame(rep(NA, ncol(X_AEGIS_exp)+2))

#perform the comparison a 100 times 

tic()
for (i in 1:100){
  print(i)
  table.perf.part <- Knockoff_comp_stat_LSCIP(scale(X_AEGIS_exp), scale = FALSE, k = 10, amplitude = 10, list_target_fdr = list_target_fdr, list_lambda = list_lambda, method = "linear")
  table.perf.linear.comp <- cbind(table.perf.linear.comp, table.perf.part[[1]])
  list.beta[[i]] <- table.perf.part[[2]]
}

toc()

table.perf.linear.comp <- table.perf.linear.comp[, -1]

# save the table in an R file you can use later to make the plots

save(table.perf.linear.comp, list.beta, file = "/home/julie/Documents/Paper_codes/New_cohorts_EXPERIMENTS/AEGIS/R_files/table_perf_statistic_comparaison_2e_cohorte.R")

######################################################################################################################################
######################################################################################################################################

############################ X_AEGIS (experiments), linear y,  methods comparison ####################################################
######################################################################################################################################

table.perf.linear.comp.meth <- data.frame(rep(NA, ncol(X_AEGIS_exp)+2))
list.beta <- list()



tic()

for (i in 1:100){
  print(i)
  table.perf.part <- Knockoff_comp_tot(scale(X_AEGIS_exp), scale = FALSE, k = 10, amplitude = 10, list_target_fdr = list_target_fdr, method = "linear") 
  table.perf.linear.comp.meth <- cbind(table.perf.linear.comp.meth, table.perf.part[[1]])
  list.beta[[i]] <- table.perf.part[[2]]
  
}

toc()

table.perf.linear.comp.meth <- table.perf.linear.comp.meth[, -1]

# save the table in an R file you can use later to make the plots

save(table.perf.linear.comp.meth, list_beta, file = "/home/julie/Documents/Paper_codes/New_cohorts_EXPERIMENTS/AEGIS/R_files/table_perf_method_comparaison_2e_cohorte.R")

######################################################################################################################################
######################################################################################################################################
