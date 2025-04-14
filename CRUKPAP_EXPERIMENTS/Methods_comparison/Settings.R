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

library(doParallel)

# get the real X

RNA.data <- read.table(file = "/home/julie/Documents/Paper_codes/Data/CRUKPAP/dat_vst_NS_no_Inel_selected.txt", head = T) 
X_real <- as.matrix(RNA.data)

n <- nrow(X_real)
p <- ncol(X_real)

# get the sequence of target FDR

list_target_fdr <- seq(0, 1, 0.01)


# Simulate gaussian data with known covariance structure

# Create gaussian matrix with bloc diagonal covariance structure
X_BC <- function(rho, p){# BC : bloc diagonal
  Sigma_mat <- matrix(data = rho, nrow = p, ncol = p)
  Sigma_mat[(round(p*0.55)+1):p, (round(p*0.55)+1):p] = rho - 0.1
  Sigma_mat[1:round(p*0.55), (round(p*0.55)+1):p] = -rho + 0.1
  Sigma_mat[(round(p*0.55)+1):p, 1:round(p*0.55)] = -rho + 0.1

  diag(Sigma_mat) <- 1
  X <- mvrnorm(369, mu = rep(0,p), Sigma_mat)
  return(list(X,Sigma_mat))
}

list_X_Sigma <- X_BC(0.7,749)
X.sim2 <- list_X_Sigma[[1]]
S.sim2 <- list_X_Sigma[[2]]



