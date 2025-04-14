library(stringr)

list.cv.ind <- list()

for (n_cv in (1:10)){
  load(str_glue('/media/julie/T5 EVO/Data_cluster/KO_bank_CV_sets/CV_indices/CV_',n_cv,'.R')) # get indices saved in the KO_banks_with_CV_cluster_script.R 
  list.cv.ind <-append(list.cv.ind, split_vec)
}

save(list.cv.ind, file = "/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Stability/R_files/list_cv_indices.R") # Save the indices of sample used in the cross validation for comparison

