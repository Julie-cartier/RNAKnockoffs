LPLR_min_oracle <- function(X, y, list_lambda){
  
  # X: data matrix of size nxp
  # y: outcomes (vector of size n)
  # list_lambda : vector of penalization values (lambda) for the selection with LASSO
  
  n <- nrow(X)
  p <- ncol(X)
  
  mat.features <- data.frame(matrix(0, ncol = 1*length(list_target_fdr), nrow = p+2 ))
  
  # Select features
  
  for (i in 1:length(list_target_fdr)){
    
    lambda <- list_lambda[i]
    
    #LASSO
    
    fitlog <- glmnet(X, y, type.measure = "deviance", family = "binomial", lambda = sort(list_lambda, decreasing = TRUE))
    coef_lasso <- coef(fitlog, s = lambda)
    selected_lasso <- which(coef_lasso[-1]!=0)
    
    predicted_y  <- predict(fitlog, s = lambda, X, type='response')
    predicted_y[predicted_y>0.5]  <- 1
    predicted_y[predicted_y<=0.5] <- 0
    
    mat.features[1, i] <- fdp(selected_lasso, beta)
    mat.features[2, i] <- power(selected_lasso, beta)
    mat.features[selected_lasso + 2, i] = 1
    
    colnames(mat.features)[i] <- str_glue("lasso_", lambda)
  }
  
  rownames(mat.features) = c("FDP", "power", c(1:p))
  
  # get the FDP and power for lambda.min : 
  
  cv.fitlog <- cv.glmnet(X, y, type.measure = "deviance", family = "binomial")
  beta.lambda.min <- coef(cv.fitlog, s = "lambda.min")
  selected_lambda.min <- which(beta.lambda.min[-1]!=0)
  
  vec.lambda.min <- c(cv.fitlog$lambda.min, fdp(selected_lambda.min, beta), power(selected_lambda.min, beta))
  
  ## get lambda oracle 
  
  new_lambda_sequence <- seq(1, 0, -0.001)
  fitlog_oracle <- glmnet(X, y, type.measure = "deviance", family = "binomial", lambda = seq(1, 0, -0.001))
  
  vec.ratio.TD.FD = c()
  
  # Find the best lambda
  
  for (lambda in new_lambda_sequence){
    coef_lasso <- coef(fitlog_oracle, s = lambda)
    selected_lasso <- which(coef_lasso[-1]!=0)
    vec.ratio.TD.FD <- c(vec.ratio.TD.FD, sum(beta[selected_lasso] != 0)/sum(beta[selected_lasso] == 0))
  }
  
  # get the power and FDP corresponding to this lambda 
  
  lambda.oracle <- new_lambda_sequence[which.max(vec.ratio.TD.FD)]
  beta.lambda.oracle <- coef(fitlog_oracle, s = lambda.oracle)
  selected_lambda.oracle <- which(beta.lambda.oracle[-1]!=0)
  vec.lambda.oracle <- c(lambda.oracle, fdp(selected_lambda.oracle, beta), power(selected_lambda.oracle, beta))
  
  return(list(mat.features, vec.lambda.oracle, vec.lambda.min))
}

