rm(list = ls())


# loading libraries

library(stringr)
library(MASS)
library(tictoc)
library(ggplot2)
library(reshape2)
library(tidyverse)


library(knockoff)
library(glmnet)
library(ggpubr) 
library(mvtnorm)
library(reticulate)
library(patchwork)
library(doParallel)



# get the real X

#RNA.data <- read.table(file = "/home/julie/Documents/Paper_codes/Data/CRUKPAP/dat_vst_NS_no_Inel_selected.txt", head = T) 
RNA.data <- read.table(file = "/home/julie/Bureau/Data/CRUKPAP/dat_vst_NS_no_Inel_selected.txt", head = T)
X_real <- as.matrix(RNA.data)
n <- nrow(X_real) # Number of samples
p <- ncol(X_real) # Number of variables


# Target FDR for KO 
list_target_fdr <- seq(0,1,0.01)

# Penalization parameters for the LASSO (for the sake of graphical representation list_lambda < 0.4 which corresponds to a limit case where no features tend to be selected)
list_lambda <- seq(0, 0.40, 0.004)

# Target FDR for the BH procedure (post wilcoxon test)
exponents <- -10:0
vec <- c(10**exponents, 5*10**exponents)
list_target_fdr_wilcoxon <- vec[-length(vec)]
list_target_fdr_wilcoxon <- list_target_fdr_wilcoxon[order(list_target_fdr_wilcoxon)]