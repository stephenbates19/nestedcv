##############################################
# standard CV
##############################################

#' Cross-validation
#'
#' Performs standard k-fold cross-validation to estimate prediction
#' error and provide a confidence interval for it.
#'
#' @param X A nxp data matrix.
#' @param Y A n-vector of responses
#' @param n_folds Number of folds, an integer.
#' @param alpha Nominal type-I error level. Must be in (0,.5)
#' @param funcs A list of containing three functions used as subroutines:
#'   \describe{
#'     \item{\code{fitter}}{A function taking X and Y and returning a model fit object.}
#'     \item{\code{predictor}}{A function taking a model fit object from above
#'     and X' and returning a vector of predictions.}
#'     \item{\code{loss}}{A function taking a vector of predictions and true values returning
#'     a vector of losses.}
#'   }
#'
#' @return
#' \describe{
#'   \item{\code{mean}}{Point estimate of prediction error: mean of CV prediction errors.}
#'   \item{\code{ci_lo}}{Lower endpoint of naive CV interval for prediction error.}
#'   \item{\code{ci_hi}}{Lower endpoint of naive CV interval for prediction error.}
#'   \item{\code{sd}}{Standard deviation of CV prediction errors.}
#'   \item{\code{raw_errors}}{Errors for each point, a vector of same length as input Y.}
#'   \item{\code{fold_id}}{The fold assignmnet for each point, a vector of same length as input Y.}
#' }
#'
#' @export
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

  return(list("err_hat" = mean(errors),
              "ci_lo" = mean(errors) - qnorm(1-alpha/2) * sd(errors) / sqrt(length(Y)),
              "ci_hi" = mean(errors) + qnorm(1-alpha/2) * sd(errors) / sqrt(length(Y)),
              "raw_mean" = mean(errors),
              "sd" = sd(errors),
              "raw_errors" = errors,
              "fold_id" = fold_id))
}
##############################################


##############################################
# nested CV
##############################################

#' nested Cross-validation
#'
#' Performs nested k-fold cross-validation, aggregating over many random splits of the data.
#' Provides a confidence interval for prediction interval that is more accurate than that of
#' standard cross-validation.
#'
#' @param X A nxp data matrix.
#' @param Y A n-vector of responses
#' @param funcs A list of containing three functions used as subroutines:
#'   \describe{
#'     \item{\code{fitter}}{A function taking X and Y and returning a model fit object.}
#'     \item{\code{predictor}}{A function taking a model fit object from above
#'     and X' and returning a vector of predictions.}
#'     \item{\code{loss}}{A function taking a vector of predictions and true values returning
#'     a vector of losses.}
#'   }
#' @param n_folds Number of folds, an integer.
#' @param reps Number of repitions of nested CV to combine. Many iterations are needed for stability.
#' @param alpha Nominal type-I error level. Must be in (0,.5)
#'
#' @return
#' \describe{
#'   \item{\code{mean}}{Point estimate of prediction error: mean of CV prediction errors.}
#'   \item{\code{ci_lo}}{Lower endpoint of naive CV interval for prediction error.d}
#'   \item{\code{ci_hi}}{Lower endpoint of naive CV interval for prediction error.d}
#'   \item{\code{sd_infl}}{The multiplicative factor needed to correct the interval width.}
#'   \item{\code{sd}}{Standard deviation of CV prediction errors.}
#'   \item{\code{raw_errors}}{Errors for each point, a vector of same length as input Y.}
#'   \item{\code{fold_id}}{The fold assignmnet for each point, a vector of same length as input Y.}
#'   \item{\code{running_sd_infl}}{The value of "sd_infl" after each repitition. This is used
#'     to evaluate how stable the estimate is and if more replicates are needed.}
#' }
#'
#' @export
nested_cv <- function(X, Y, funcs, reps = 50, n_folds = 10,  alpha = .1, bias_reps = NA) {
  #compute out-of-fold errors on SE scale
  var_pivots <- c()
  ho_errs <- c()
  for(i in 1:reps) {
    temp <- nestedcv:::nested_cv_helper(X, Y, funcs, n_folds)
    var_pivots <- rbind(var_pivots, temp$pivots)
    ho_errs <- c(ho_errs, temp$errs)
  }

  # inflation estimate on a variance scale
  infl_est <- (var(as.vector(var_pivots)) / var(ho_errs) * length(Y) / n_folds - 1) * (n_folds - 1)
  infl_est <- max(1, min(infl_est, n_folds))

  # look at the estimate of inflation after each repetition
  infl_est2 <- sapply(1:reps, function(i) {
      temp <- (var(as.vector(var_pivots[1:(i*n_folds)])) / var(ho_errs[1:(i*n*(n_folds-1))]) * length(Y) / n_folds - 1) * (n_folds - 1)
      max(1, min(temp, n_folds))
  })

  #bias correction
  cv_means <- c() #estimated pred error from normal CV
  bias_est <- 0
  if(is.na(bias_reps)) {
    bias_reps <- ceiling(reps / 10) #fewer reps for bias estimation
  }
  if(bias_reps == 0) {
    bias_est <- 0
  }
  else {
    for(i in 1:bias_reps) {
      temp <- nestedcv:::naive_cv(X, Y, funcs, n_folds)
      cv_means <- c(cv_means, temp$err_hat)
    }

    bias_est <- ( mean(ho_errs) - mean(cv_means)) * (1 + ((n_folds - 2) / (n_folds - 1))^(1.5))
  }
  pred_est <- mean(ho_errs) - bias_est #debiased estimate

  list("sd_infl" = sqrt(infl_est),
      "err_hat" = pred_est,
      "ci_lo" = pred_est - qnorm(1-alpha/2) * sd(ho_errs) / sqrt(length(Y)) * sqrt(infl_est),
      "ci_hi" = pred_est + qnorm(1-alpha/2) * sd(ho_errs) / sqrt(length(Y)) * sqrt(infl_est),
      "raw_mean" = mean(ho_errs),
      "bias_est" = bias_est,
      "sd" = sd(ho_errs),
      "running_sd_infl" = sqrt(infl_est2))
}

#' Internal helper for nested_cv
#'
#' Randomly assigns folds and does one iteration of nested cross-validation.
#'
#' arguments as in the "nested_cv" function
#'
#' @return
#' \describe{
#'   \item{\code{pivots}}{A vector of same length as Y of difference statistics.}
#'   \item{\code{errors}}{A vector of all errors of observations not used in model training.}
#' }
nested_cv_helper <- function(X, Y, funcs, n_folds = 10) {
  fold_id <- 1:nrow(X) %% n_folds + 1
  fold_id <- sample(fold_id)

  #nested CV model fitting
  ho_errors <- array(0, dim = c(n_folds, n_folds, nrow(X) / n_folds))
    #^ entry i, j is error on fold i,
    # when folds i & j are not used for model fitting.
  for(f1 in 1:(n_folds - 1)) {
    for(f2 in (f1+1):n_folds) {
      test_idx <- c(which(fold_id == f1), which(fold_id == f2))
      fit <- funcs$fitter(X[-test_idx, ], Y[-test_idx])
      preds <- funcs$predictor(fit, X)
      ho_errors[f1, f2, ] <- funcs$loss(preds[fold_id == f1], Y[fold_id == f1])
      ho_errors[f2, f1, ] <- funcs$loss(preds[fold_id == f2], Y[fold_id == f2])
    }
  }

  #e_bar - f_bar in the notation of the paper
  out_vec <- rep(0, n_folds)
  for(f1 in 1:(n_folds)) {
    i <- sample(1:(n_folds-1), 1)
    if(i >= f1) {i <- i + 1} #sample random other fold
    out_vec[f1] <- sum(ho_errors[, f1, ]) / (nrow(X) * (n_folds - 1) / n_folds) -  mean(ho_errors[f1, i, ])
  }

  #errors on points not used for fitting, combined across all runs
  all_ho_errs <- c()
  for(f1 in 1:(n_folds - 1)) {
    for(f2 in (f1+1):n_folds) {
      all_ho_errs <- c(all_ho_errs, ho_errors[f1, f2, ], ho_errors[f2, f1, ])
    }
  }

  return(list("pivots" = out_vec,
              "errs" = all_ho_errs))
}
##############################################
