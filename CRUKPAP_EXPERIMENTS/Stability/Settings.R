rm(list = ls())


# loading libraries

library(stringr)
library(MASS)
library(tictoc)
library(ggplot2)
library(reshape2)
library(tidyverse)

library(PRROC)
library(knockoff)
library(glmnet)
library(ggpubr) 
library(mvtnorm)

library(doParallel)

# get the real X

RNA.data <- read.table(file = "/home/julie/Documents/Paper_codes/Data/CRUKPAP/dat_vst_NS_no_Inel_selected.txt", head = T) 
X_real <- as.matrix(RNA.data)


n <- nrow(X_real) # Number of samples
p <- ncol(X_real) # Number of variables

list_target_FDP_fdr = seq(0,1,0.01) 


