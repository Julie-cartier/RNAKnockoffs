source("/home/julie/Documents/Paper_codes/CUKPAP_EXPERIMENTS/Statistics_comparison/Settings.R")
source("/home/julie/Documents/Paper_codes/CUKPAP_EXPERIMENTS/Statistics_comparison/Functions.R")

load("/home/julie/Documents/Paper_codes/CUKPAP_EXPERIMENTS/Statistics_comparison/R_files/list_beta_interaction.R")

### Experiment for Tree_based_statistics files 

require(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/home/julie/anaconda3/envs/py-environment/bin/python")
knockpy <- import("knockpy", convert = FALSE)
sklearn <- import("sklearn", convert = FALSE)
numpy <- import("numpy", convert = FALSE)


############################ X_CRUKPAP, linear/interaction y,  LASSO-based statistics comparison #####################################
######################################################################################################################################

# linear
list.beta <- list()
table.perf.linear.LASSO <- data.frame(rep(NA, ncol(X_real)+2))

tic()
for (i in 1:10){# Repeat the comparison for 10 simulated outcomes
  print(i)
  table.perf.part <- Knockoff_LassoBasedMethod(scale(X_real, center = TRUE, scale = TRUE), scale = FALSE, k = 10, amplitude = 10, list_target_fdr = list_target_fdr, method = "linear")
  table.perf.linear.LASSO <- cbind(table.perf.linear.LASSO, table.perf.part[[1]])
  list.beta[[i]] <- table.perf.part[[2]]
}
toc()

table.perf.linear.LASSO <- table.perf.linear.LASSO[, -1]

# non-linear

table.perf.nlinear.LASSO <- data.frame(rep(NA, ncol(X_real)+2))

tic()
for (i in 1:10){# Repeat the comparison for 10 simulated outcomes
  print(i)
  y <- list_y_interaction[[i]]
  beta <- list_beta_interaction[[i]]
  table.perf.part <- Knockoff_LassoBasedMethod(scale(X_real, center = TRUE, scale = TRUE), y = y, beta = beta, list_target_fdr = list_target_fdr)
  table.perf.nlinear.LASSO <- cbind(table.perf.nlinear.LASSO, table.perf.part)
}

toc()

table.perf.nlinear.LASSO <- table.perf.nlinear.LASSO[, -1]


save(table.perf.linear.LASSO, table.perf.nlinear.LASSO, list.beta, file = "/home/julie/Documents/Paper_codes/CUKPAP_EXPERIMENTS/Statistics_comparison/R_files/LASSOStat_trueX_fakeY_tables_perfx10.R")

######################################################################################################################################
######################################################################################################################################

############################ X_CRUKPAP, linear/interaction y,  Tree-based statistics comparison ######################################
######################################################################################################################################

# linear

list.beta <- list()
table.perf.linear.RF.synth <- data.frame(rep(NA, ncol(X_real)+2))

tic()
for (i in 1:10){# Repeat the comparison for 10 simulated outcomes
  print(i)
  table.perf.part <- Knockoff_TreeBasedMethod.synthesis(scale(X_real, center = TRUE, scale = TRUE), scale = FALSE, k = 10, amplitude = 10, list_target_fdr = list_target_fdr, method = "linear")
  table.perf.linear.RF.synth <- cbind(table.perf.linear.RF.synth, table.perf.part[[1]])
  list.beta[[i]] <- table.perf.part[[2]]
}

toc()

table.perf.linear.RF.synth <- table.perf.linear.RF.synth[, -1]

# non-linear

list.beta.nl = list()
table.perf.nlinear.RF.synth <- data.frame(rep(NA, ncol(X_real)+2))

tic()
for (i in 1:10){# Repeat the comparison for 10 simulated outcomes
  print(i)
  y <- list_y_interaction[[i]]
  beta <- list_beta_interaction[[i]]
  table.perf.part <- Knockoff_TreeBasedMethod.synthesis(scale(X_real, center = TRUE, scale = TRUE), y = y, beta = beta, list_target_fdr = list_target_fdr)
  table.perf.nlinear.RF.synth <- cbind(table.perf.nlinear.RF.synth, table.perf.part)
}

toc()

table.perf.nlinear.RF.synth <- table.perf.nlinear.RF.synth[, -1]

save(table.perf.linear.RF.synth, table.perf.nlinear.RF.synth, list.beta, file = "/home/julie/Documents/Paper_codes/CUKPAP_EXPERIMENTS/Statistics_comparison/R_files/TreeSynthStat_trueX_fakeY_tables_perfx10.R")


############################ X_CRUKPAP, linear/interaction y,  Deep-based statistics comparison ######################################
######################################################################################################################################

# linear

list.beta <- list()
table.perf.linear.Deep <- data.frame(rep(NA, ncol(X_real)+2))

tic()
for (i in 1:10){# Repeat the comparison for 10 simulated outcomes
  print(i)
  table.perf.part <- Knockoff_DeepBasedMethod(scale(X_real, center = TRUE, scale = TRUE), scale = FALSE, k = 10, amplitude = 10, list_target_fdr = list_target_fdr, method = "linear")
  table.perf.linear.Deep <- cbind(table.perf.linear.Deep, table.perf.part[[1]])
  list.beta[[i]] <- table.perf.part[[2]]
}

toc()

table.perf.linear.Deep <- table.perf.linear.Deep[, -1]

# non-linear

list.beta.nl = list()
table.perf.nlinear.Deep <- data.frame(rep(NA, ncol(X_real)+2))

tic()
for (i in 1:10){# Repeat the comparison for 10 simulated outcomes
  print(i)
  y <- list_y_interaction[[i]]
  beta <- list_beta_interaction[[i]]
  table.perf.part <- Knockoff_DeepBasedMethod(scale(X_real, center = TRUE, scale = TRUE), y = y, beta = beta, list_target_fdr = list_target_fdr)
  table.perf.nlinear.Deep <- cbind(table.perf.nlinear.Deep, table.perf.part)
}

toc()

table.perf.nlinear.Deep <- table.perf.nlinear.Deep[, -1]


save(table.perf.linear.Deep, table.perf.nlinear.Deep, list.beta, file = "/home/julie/Documents/Paper_codes/CUKPAP_EXPERIMENTS/Statistics_comparison/R_files/DeepStat_trueX_fakeY_tables_perfx10.R")


############################ X_CRUKPAP, linear/interaction y , Statistics comparison #################################################
######################################################################################################################################


# linear
list.beta <- list()
table.perf.linear.comp <- data.frame(rep(NA, ncol(X_real)+2))

tic()
for (i in 1:100){# Repeat the comparison for 100 simulated outcomes
  print(i)
  table.perf.part <- Knockoff_comp_stat_LSCIP(scale(X_real, center = TRUE, scale = TRUE), scale = FALSE, k = 10, amplitude = 10, list_target_fdr = list_target_fdr, list_lambda = list_lambda, method = "linear")
  table.perf.linear.comp <- cbind(table.perf.linear.comp, table.perf.part[[1]])
  list.beta[[i]] <- table.perf.part[[2]]
}

toc()

table.perf.linear.comp <- table.perf.linear.comp[, -1]

# non-linear

list.beta.nl = list()
table.perf.nlinear.comp <- data.frame(rep(NA, ncol(X_real)+2))

tic()

for (i in 1:10){# Repeat the comparison for 10 simulated outcomes
  print(i)
  y <- list_y_interaction[[i]]
  beta <- list_beta_interaction[[i]]
  table.perf.part <- Knockoff_comp_stat_LSCIP(scale(X_real, center = TRUE, scale = TRUE), y = y, beta = beta, list_target_fdr = list_target_fdr, list_lambda = list_lambda)
  table.perf.nlinear.comp <- cbind(table.perf.nlinear.comp, table.perf.part)
}

toc()

table.perf.nlinear.comp <- table.perf.nlinear.comp[, -1]


save(table.perf.linear.comp, table.perf.nlinear.comp, list.beta, file = "/home/julie/Documents/Paper_codes/CUKPAP_EXPERIMENTS/Statistics_comparison/R_files/compStat_trueX_fakeY_tables_perf.R")


######################################################################################################################################
######################################################################################################################################
