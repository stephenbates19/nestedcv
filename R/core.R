##############################################
# standard CV
##############################################
naive_cv <- function(X, Y, funcs, n_folds = 10, alpha = .1) {
  fold_id <- (1:nrow(X)) %% n_folds + 1
  fold_id <- sample(fold_id)

  errors <- c()
  for(k in 1:n_folds) {
    fit <- funcs$fitter(X[fold_id !=k, ], Y[fold_id != k])
    y_hat <- funcs$predictor(fit, X[fold_id == k, ])
    error_k <- funcs$loss(y_hat, Y[fold_id == k])
    errors <- c(errors, error_k)
  }

  return(c(mean(errors), sd(errors)))
}
##############################################


##############################################
# nested CV
##############################################
nested_cv <- function(X, Y, funcs, reps = 50, n_folds = 10,  alpha = .1, trans = list(identity), running = FALSE) {
  #compute out-of-fold errors on SE scale
  var_pivots <- c()
  ho_errs <- c()
  for(j in 1:reps) {
    temp <- nested_cv_helper(X, Y, funcs, n_folds, alpha, trans)
    var_pivots <- rbind(var_pivots, temp$pivots)
    ho_errs <- c(ho_errs, temp$errs)
  }

  #look at the estimate after each rep to evaluate the variance
  if(running == TRUE) {
    infl_est <- sapply(1:reps * n_folds, function(i) {
      (var(as.vector(var_pivots[1:i])) / var(ho_errs[1:i]) * length(Y) / n_folds - 1) * (n_folds)
    })
    return(infl_est)
  }

  infl_est <- (var(as.vector(var_pivots)) / var(ho_errs) * length(Y) / n_folds - 1) * (n_folds)
  max(1, min(infl_est, n_folds))
}


nested_cv_helper <- function(X, Y, funcs, n_folds = 10,  alpha = .1, percentile = FALSE, trans = list(identity)) {
  fold_id <- 1:nrow(X) %% n_folds + 1
  fold_id <- sample(fold_id)

  #nested CV model fitting
  ho_errors <- array(0, dim = c(n_folds, n_folds, nrow(X) / n_folds))
  for(f1 in 1:(n_folds - 1)) {
    for(f2 in (f1+1):n_folds) {
      test_idx <- c(which(fold_id == f1), which(fold_id == f2))
      fit <- funcs$fitter(X[-test_idx, ], Y[-test_idx])
      preds <- funcs$predictor(fit, X)
      ho_errors[f1, f2, ] <- funcs$loss(preds[fold_id == f1], Y[fold_id == f1])
      ho_errors[f2, f1, ] <- funcs$loss(preds[fold_id == f2], Y[fold_id == f2])
    }
  }

  out_vec <- rep(0, n_folds)
  for(f1 in 1:(n_folds)) {
    i <- sample(1:(n_folds-1), 1)
    if(i >= f1) {i <- i + 1} #sample random other fold
    out_vec[f1] <- sum(ho_errors[, f1, ]) / (nrow(X) * (n_folds - 1) / n_folds) -  mean(ho_errors[f1, i, ])
  }

  all_ho_errs <- c()
  for(f1 in 1:(n_folds - 1)) {
    for(f2 in (f1+1):n_folds) {
      all_ho_errs <- c(all_ho_errs, ho_errors[f1, f2, ], ho_errors[f2, f1, ])
    }
  }

  return(list("pivots" = out_vec, "errs" = all_ho_errs))
}
##############################################
