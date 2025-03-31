rm(list = ls())


# loading libraries

library(stringr)
library(MASS)
library(tictoc)
library(ggplot2)
library(reshape2)
library(tidyverse)


#library(Metrics)
library(randomForest)
library(knockoff)
library(glmnet)
library(mboost)
library(ggpubr) 
library(mvtnorm)

library(doParallel)

library(xgboost)
library(KOBT)

# get the real X

RNA.data <- read.table(file = "/home/julie/Documents/Paper_codes/Data/CRUKPAP/dat_vst_NS_no_Inel_selected.txt", head = T) 
X_real <- as.matrix(RNA.data)

#settings

n <- nrow(X_real) # Number of samples
p <- ncol(X_real) # Number of variables

# For KO selection

list_target_fdr <- seq(0,1,0.01)

# For lasso selection

list_lambda <- seq(0, 0.40, 0.004)
