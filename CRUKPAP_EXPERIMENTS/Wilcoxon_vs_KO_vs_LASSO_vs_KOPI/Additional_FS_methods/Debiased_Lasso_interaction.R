# This code reuses example and functions from https://github.com/luxia-bios/DebiasedLassoGLMs/blob/main/DBL_GLMs_functions.R

library(mvtnorm)
library(glmnet)
library(tictoc)
library(doParallel)
library(stringr)

registerDoParallel(10)
set.seed(811)
### Settings ###

load("/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Statistics_comparison/R_files/list_beta_interaction.R")
RNA.data <- read.table(file = "/home/julie/Bureau/Data/CRUKPAP/dat_vst_NS_no_Inel_selected.txt", head = T)
X_real <- as.matrix(RNA.data)

### Functions ###



########## https://github.com/luxia-bios/DebiasedLassoGLMs/blob/main/DBL_GLMs_functions.R ###########
### The repository documents the R code accompanying the article "De-biased Lasso forGeneralized ####
##### Linear Models with A Diverging Number of Covariates" by Lu Xia, Bin Nan and Yi Li (2021+). ####

ORIG_DS_inf <- function(x, y, family, lasso_est, nfold=5, n_lambda=100, lambda_ratio=0.005){
  nn <- length(y)
  pp <- ncol(x)
  X <- cbind(rep(1, nrow(x)), x)
  if(ncol(X) != length(lasso_est)) {
    stop("The length of lasso_est is incompatible with the covariate matrix.")
  }
  if(family == "binomial") {
    mu <- as.vector(exp(X%*%lasso_est)/(1+exp(X%*%lasso_est)))
    neg_dloglik_glmnet <- 0 - as.vector(t(X)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- t(X)%*%diag(mu*(1-mu))%*%X/nn
    C_glmnet <- sqrt(diag(mu*(1-mu))/nn)%*%X
  } else if(family == "poisson") {
    mu <- as.vector(exp(X%*%lasso_est))
    neg_dloglik_glmnet <- 0 - as.vector(t(X)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- t(X)%*%diag(mu)%*%X/nn
    C_glmnet <- sqrt(diag(mu)/nn)%*%X
  } else {
    stop("Input family is not supported.")
  }
  
  theta_glmnet <- diag(pp+1)
  tau_glmnet <- rep(NA, pp+1)
  for(j in 1:(pp+1)) { # for: nodewise lasso
    current_x <- sqrt(nn)*C_glmnet[,-j]
    current_y <- sqrt(nn)*as.vector(C_glmnet[,j])
    lam_max <- max(abs(t(current_x)%*%current_y)/nn)
    lam_min <- lam_max*lambda_ratio
    lam_seq <- exp(seq(from=log(lam_max), to=log(lam_min), length.out=n_lambda))
    gamma_j_glmnet <- cv.glmnet(x=current_x, y=current_y,
                                family="gaussian", alpha=1, standardize=F, intercept=F,
                                nfolds=nfold, lambda=lam_seq, parallel = TRUE)
    gamma_j_glmnet <- as.vector(glmnet(x=sqrt(nn)*C_glmnet[,-j], y=sqrt(nn)*as.vector(C_glmnet[,j]),
                                       family="gaussian", alpha=1, standardize=F, intercept=F,
                                       lambda=gamma_j_glmnet$lambda.min)$beta)
    theta_glmnet[j,-j] <- (-1)*t(gamma_j_glmnet)
    tau_glmnet[j] <- as.numeric(neg_ddloglik_glmnet[j,j]-neg_ddloglik_glmnet[j,-j]%*%gamma_j_glmnet)
  } # end for: nodewise lasso
  theta_glmnet <- diag(1/tau_glmnet)%*%theta_glmnet 
  
  b_hat_nw <- as.vector(beta_glmnet - theta_glmnet%*%neg_dloglik_glmnet)
  se_nw <- sqrt(diag(theta_glmnet%*%neg_ddloglik_glmnet%*%t(theta_glmnet)))/sqrt(nn)
  pval_nw <- 2*pnorm(abs(b_hat_nw/se_nw), lower.tail=F)
  
  return(list(est=b_hat_nw, se=se_nw, pvalue=pval_nw, theta=theta_glmnet))
}


### Experiments ###
pvalues.list <- list()
tic()
for (i in 1:10){
print(i)
  y <- list_y_interaction[[i]]
  beta <- list_beta_interaction[[i]]
  X <- scale(X_real)
  cvobj_glmnet <- cv.glmnet(x=X, y=y, family="binomial", type.measure = "deviance", standardize=F, intercept=T)
  beta_glmnet <- as.vector(coef(glmnet(x=X, y=y, family="binomial", standardize=F, intercept=T), s=cvobj_glmnet$lambda.min))
  
  obj_orig_ds <- ORIG_DS_inf(x=X, y=y, family="binomial", lasso_est=beta_glmnet)
  
  pvalues.list[[i]] <- obj_orig_ds[[3]]
}
toc()

#save(pvalues.list, file = "/cluster/CBIO/data1/jcartier/DebiasedLasso/R_files/pvalues.RData")
save(pvalues.list, file = "/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Wilcoxon_vs_KO_vs_LASSO_vs_KOPI/Additional_FS_methods/R_files/pvalues_interaction.RData")


