source('/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Wilcoxon_vs_KO_vs_LASSO_vs_KOPI/Settings.R')
source('/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Wilcoxon_vs_KO_vs_LASSO_vs_KOPI/Functions.R')
load("/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Wilcoxon_vs_KO_vs_LASSO_vs_KOPI/R_files/LPLR_KO_BH_linear_table_perfx10.R")


library(caret)
library(e1071)
library(tictoc)
# clbk("mlr3fselect.svm_rfe")
# #> <CallbackBatchFSelect:mlr3fselect.svm_rfe>: SVM-RFE Callback
# #> * Active Stages: on_optimization_begin
# 
# library(mlr3learners)
# 
# # Create instance with classification svm with linear kernel
# instance = fsi(
#   task = tsk("sonar"),
#   learner = lrn("classif.svm", type = "C-classification", kernel = "linear"),
#   resampling = rsmp("cv", folds = 3),
#   measures = msr("classif.ce"),
#   terminator = trm("none"),
#   callbacks = clbk("mlr3fselect.svm_rfe"),
#   store_models = TRUE
# )
# 
# fselector = fs("rfe", feature_number = 5, n_features = 10)
# 
# # Run recursive feature elimination on the Sonar data set
# fselector$optimize(instance)



#### Test Bayesian variable selection conditional Logistic regression ###
i = 1
beta <- list.beta[[i]]
y <- y.sample(scale(X_real), beta, method = "linear", scale = FALSE, threshold = 0.5)

library(mRMRe)
df = data.frame(y = y, scale(X_real))
data <- new("mRMRe.Data",(as.data.frame(cbind(y,scale(X_real)))))
data$
filter <- mRMR.classic("mRMRe.Filter", data = data, target_indices = 1, feature_count = 749)


dd <- mRMR.data(data = df)
test <- mRMR.classic(data = dd, target_indices = c(1), feature_count = 10)

test@scores



get_power_FDP <- function(feature_importance, beta){
  
  FDP <- c()
  power <- c()
  
  rank.value.ind <- order(feature_importance, decreasing = TRUE)
  selection.state <- beta[rank.value.ind]
  k <- length(which(beta!=0))
  
  for (n in 1:length(beta)){
    selected.set <- selection.state[1:n]
    power <- c(power, length(which(selected.set != 0))/k)
    FDP <- c(FDP, length(which(selected.set == 0))/n)
  }
  
  return(rbind(FDP, power))
}
MRMR_importance <- unlist(test@scores)
test2 <- get_power_FDP(MRMR_importance, list.beta[[i]])

test3 <- MRMR(scale(X_real), as.factor(y), positive = FALSE, threads = 4)
