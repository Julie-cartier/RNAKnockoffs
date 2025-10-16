source('/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Wilcoxon_vs_KO_vs_LASSO_vs_KOPI/Settings.R')
source('/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Wilcoxon_vs_KO_vs_LASSO_vs_KOPI/Functions.R')
load("/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Wilcoxon_vs_KO_vs_LASSO_vs_KOPI/R_files/LPLR_KO_BH_linear_table_perfx10.R")


library(caret)
library(e1071)

clbk("mlr3fselect.svm_rfe")
#> <CallbackBatchFSelect:mlr3fselect.svm_rfe>: SVM-RFE Callback
#> * Active Stages: on_optimization_begin

library(mlr3learners)

# Create instance with classification svm with linear kernel
instance = fsi(
  task = tsk("sonar"),
  learner = lrn("classif.svm", type = "C-classification", kernel = "linear"),
  resampling = rsmp("cv", folds = 3),
  measures = msr("classif.ce"),
  terminator = trm("none"),
  callbacks = clbk("mlr3fselect.svm_rfe"),
  store_models = TRUE
)

fselector = fs("rfe", feature_number = 5, n_features = 10)

# Run recursive feature elimination on the Sonar data set
fselector$optimize(instance)



#### Test Bayesian variable selection conditional Logistic regression ###

library(varbvs)
fit.varbvs <- varbvs(X = scale(X_real), y = y, Z = NULL, family = "binomial")
