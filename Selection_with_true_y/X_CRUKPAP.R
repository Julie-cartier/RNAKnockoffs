# get the data matrix 

RNA.data <- read.table(file = "/home/julie/Bureau/Data/CRUKPAP/dat_vst_NS_no_Inel_selected.txt", head = T) 

X_CRUKPAP <- as.matrix(RNA.data)

# get the real cancer status

Clinic_data_CRUKPAP_tot <- read.table(file = "/home/julie/Bureau/Data/CRUKPAP/clinical_no_Inel.tsv", head = TRUE, sep = "\t")

y_CRUKPAP <- Clinic_data_CRUKPAP_tot$cancernum

save(X_CRUKPAP, y_CRUKPAP, file = "/home/julie/Documents/Paper_codes/Selection_with_true_y/R_files/CRUKPAP_data.R")

# load 100 KO matrices from the bank created with /home/julie/Documents/Paper_codes/CUKPAP_EXPERIMENTS/Wilcoxon_vs_KO_vs_LASSO_vs_KOPI/KO_banks_cluster_script.R script

load(str_glue("/media/julie/T5 EVO/Data_cluster/KO_bank_under_sets/X_k_",10,"_100.R"))


save(list_X_k_100, file = "/home/julie/Documents/Paper_codes/Data/CRUKPAP/KO_partial_bank_100_CRUKPAP_LSCIP/KO_bank_CRUKPAP_real")
