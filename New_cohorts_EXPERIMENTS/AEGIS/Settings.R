rm(list = ls()) 


library(stringr)
library(MASS)
library(tictoc)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(randomForest)
library(knockoff)
library(glmnet)
library(mboost)
library(ggpubr) 
library(mvtnorm)
library(doParallel)
library(xgboost)
library(reticulate)



# load CRUKPAP data

RNA.data <- read.table(file = "/home/julie/Bureau/Data/CRUKPAP/dat_vst_NS_no_Inel_selected.txt", head = T) 

X_CRUKPAP <- as.matrix(RNA.data)


# create the AEGIS data from the set of the 749 risk genes 

RNA_data_perez_tot <- read.table(file = "/home/julie/Documents/Paper_codes/Data/AEGIS/perez2017_exprmat_geneID_filtered_and_ENS_converted.tsv",head = TRUE, sep = "\t") 

X_AEGIS <- t(as.matrix(RNA_data_perez_tot))

X_AEGIS_exp <- X_AEGIS[,intersect(colnames(X_CRUKPAP), colnames(X_AEGIS))]


# List of target FDR for the knockoffs
list_target_fdr <- seq(0,1,0.01)

#list of lambda for the lasso penalized logistic regression
list_lambda <- seq(0, 0.40, 0.004)


