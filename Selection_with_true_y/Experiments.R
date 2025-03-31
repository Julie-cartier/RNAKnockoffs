## get R functions and packages ##

source("/home/julie/Documents/Paper_codes/Selection_with_true_y/Libraries_Functions.R")



## get py functions ##

# Set the py environment
Sys.setenv(RETICULATE_PYTHON = "/home/julie/anaconda3/envs/py-environment/bin/python")

# Import necessary py modules
knockpy <- import("knockpy", convert = FALSE)
sklearn <- import("sklearn", convert = FALSE)
numpy <- import("numpy", convert = FALSE)

source_python("/home/julie/Documents/Paper_codes/CUKPAP_EXPERIMENTS/Wilcoxon_vs_KO_vs_LASSO_vs_KOPI/Functions.py")


## get data for the BC cohort ##

load("/home/julie/Documents/Paper_codes/Selection_with_true_y/R_files/BC_data.R")


##################################
#### Eperiments for BC cohort ####
##################################


## Create D = 100 KO matrices ##

# This has been coomputed on a cluster
# D = 100
# list.KO.BC <- list()
# 
# tic()
# for (i in 1:D){
#   print(i)
#   list.KO.BC[[i]] <- create.LSCIP_parallel(X_BC)
# }
# toc()

load("/home/julie/Documents/Paper_codes/Data/BC/KO_BC_bank.R")

# get the table of statistics for all the KO matrices (computed once fall all classifiers)

tic()
KO_stats_BC <- create_KO_stats(X_BC, list.KO.BC, y_BC.er, num_cores = 5)
toc()


# Perform variable selection with KOPI for three different biomarkers : HER2 ER PGR

tic()
table.selected.genes.BC.er <- KOPI_selection(X_BC, y_BC.er, list.KO.BC, KO_stats = KO_stats_BC)
toc()

colnames(X_BC)[which(table.selected.genes.BC.er$`q = 0.1` != 0)]
colnames(X_BC)[which(table.selected.genes.BC.er$`q = 0.2` != 0)]
colnames(X_BC)[which(table.selected.genes.BC.er$`q = 0.3` != 0)]
colnames(X_BC)[which(table.selected.genes.BC.er$`q = 0.5` != 0)]

#save(table.selected.genes.BC.er, file = "/home/julie/Documents/Paper_codes/Selection_with_true_y/R_files/BC_er_selected_genes.R")


### hr2 ###

# compute steps separetaly to avoid python issue 
tic()
KO_stats_BC_hr2 <- create_KO_stats(X_BC, list.KO.BC, y_BC.hr2, num_cores = 5)
toc()



tic()
table.selected.genes.BC.hr2 <- KOPI_selection(X_BC, y_BC.hr2, list.KO.BC, KO_stats = KO_stats_BC_hr2)
toc()

colnames(X_BC)[which(table.selected.genes.BC.hr2$`q = 0.1` != 0)]
colnames(X_BC)[which(table.selected.genes.BC.hr2$`q = 0.2` != 0)]
colnames(X_BC)[which(table.selected.genes.BC.hr2$`q = 0.3` != 0)]
colnames(X_BC)[which(table.selected.genes.BC.hr2$`q = 0.5` != 0)]

#save(table.selected.genes.BC.hr2, file = "/home/julie/Documents/Paper_codes/Selection_with_true_y/R_files/BC_hR2_selected_genes.R")


### pgr ###

# compute steps separetaly to avoid python issue 
tic()
KO_stats_BC_pgr <- create_KO_stats(X_BC, list.KO.BC, y_BC.pgr, num_cores = 5)
toc()


tic()
table.selected.genes.BC.pgr <- KOPI_selection(X_BC, y_BC.pgr, list.KO.BC, KO_stats = KO_stats_BC_pgr)
toc()

#save(table.selected.genes.BC.pgr, file = "/home/julie/Documents/Paper_codes/Selection_with_true_y/R_files/BC_pgr_selected_genes.R")

colnames(X_BC)[which(table.selected.genes.BC.pgr$`q = 0.1` != 0)]
colnames(X_BC)[which(table.selected.genes.BC.pgr$`q = 0.2` != 0)]
colnames(X_BC)[which(table.selected.genes.BC.pgr$`q = 0.3` != 0)]
colnames(X_BC)[which(table.selected.genes.BC.pgr$`q = 0.5` != 0)]


##################################
#### Eperiments for AEGIS cohort #
##################################

## get data for the AEGIS cohort ##

load("/home/julie/Documents/Paper_codes/Selection_with_true_y/R_files/AEGIS_data.R")

## Create D = 100 KO matrices ##

D = 100
list.KO.AEGIS <- list()

tic()
for (i in 1:D){
  print(i)
  list.KO.AEGIS[[i]] <- create.LSCIP_parallel(X_AEGIS)
}
toc()


# Perform variable selection with KOPI for cancer status classification
tic()
table.selected.genes.AEGIS <- KOPI_selection(X_AEGIS, y_AEGIS, list.KO.AEGIS)
toc()


#save(table.selected.genes.AEGIS, file = "/home/julie/Documents/Paper_codes/Selection_with_true_y/R_files/AEGIS_selected_genes.R")


colnames(X_AEGIS)[which(table.selected.genes.AEGIS$`q = 0.1` != 0)]
colnames(X_AEGIS)[which(table.selected.genes.AEGIS$`q = 0.2` != 0)]
colnames(X_AEGIS)[which(table.selected.genes.AEGIS$`q = 0.3` != 0)]
colnames(X_AEGIS)[which(table.selected.genes.AEGIS$`q = 0.5` != 0)]




##################################
## Eperiments for CRUKPAP cohort #
##################################

## get data for the AEGIS cohort ##

load("/home/julie/Documents/Paper_codes/Selection_with_true_y/R_files/CRUKPAP_data.R")
load("/home/julie/Bureau/Data/CRUKPAP/KO_bank_CRUKPAP_real")

## Create D = 100 KO matrices ##

# these 100 matrices have already been created for KOPI experiments


# Perform variable selection with KOPI for cancer status classification
tic()
table.selected.genes.CRUKPAP <- KOPI_selection(X_CRUKPAP, y_CRUKPAP, list_X_k_100)
toc()


#save(table.selected.genes.CRUKPAP, file = "/home/julie/Documents/Paper_codes/Selection_with_true_y/R_files/CRUKPAP_selected_genes.R")


colnames(X_CRUKPAP)[which(table.selected.genes.CRUKPAP$`q = 0.1` != 0)]
colnames(X_CRUKPAP)[which(table.selected.genes.CRUKPAP$`q = 0.2` != 0)]
colnames(X_CRUKPAP)[which(table.selected.genes.CRUKPAP$`q = 0.3` != 0)]
colnames(X_CRUKPAP)[which(table.selected.genes.CRUKPAP$`q = 0.5` != 0)]





