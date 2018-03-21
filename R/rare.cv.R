#' Perform K-fold cross validation
#'
#' The function does K-fold cross validaton (CV) to choose an optimal pair of (\code{lambda}, \code{alpha})
#' on which the model performs best according to the chosen error metric: mean squared error
#' or mean absolute error.
#'
#' @param fitObj Output of \code{rarefit}
#' @param y Response variable.
#' @param X \code{nobs}-by-\code{nvars} input matrix:
#' each row is an observation vector and each column stores
#' a count covariate.
#' @param errtype Type of error metric used in cross validation.
#' Available choices are \emph{mean-squared-error} (default)
#' and \emph{mean-absolute-error}.
#' @param nfolds Number of folds (default is 5)
#' @param ... Other arguments that can be passed to \code{rarefit}
#'
#' @return
#' \item{folds}{A length-\code{nfolds} list with the kth element being elements in the \code{k}th fold.}
#' \item{errs}{A \code{nlam}-by-\code{nalpha}-by-\code{nfolds} 3-dimensional array of errors.
#' \code{errs[i,j,k]} is error incurred in using \code{lambda[i]} and \code{alpha[j]} on the \code{k}th fold.}
#' \item{m}{A \code{nlam}-by-\code{nalpha} matrix for storing CV error (i.e., mean error across folds).
#' \code{m[i,j]} is CV error incurred in using \code{lambda[i]} and \code{alpha[j]}.}
#' \item{se}{A \code{nlam}-by-\code{nalpha} matrix for storing standard error across folds.
#' \code{se[i,j]} is standard error incurred in using \code{lambda[i]} and \code{alpha[j]}.}
#' \item{ibest}{Indices of pair of (\code{lambda}, \code{alpha}) minimizing CV error.}
#' \item{lambda.best}{Value of \code{lambda} minimizing CV error.}
#' \item{alpha.best}{Value of \code{alpha} minimizing CV error.}
#'
#' @examples
#' \dontrun{
#' # See vignette for more details.
#' set.seed(100)
#' ts <- sample(1:length(data.rating), 400) # Train set indices
#' # Fit the model on train set
#' ourfit <- rarefit(y = data.rating[ts], X = data.dtm[ts, ], hc = data.hc, lam.min.ratio = 1e-6,
#'                   nlam = 20, nalpha = 10, rho = 0.01, eps1 = 1e-5, eps2 = 1e-5, maxite = 1e4)
#' # Cross validation
#' ourfit.cv <- rarefit.cv(ourfit, y = data.rating[ts], X = data.dtm[ts, ],
#'                         rho = 0.01, eps1 = 1e-5, eps2 = 1e-5, maxite = 1e4)
#' }
#'
#' @seealso \code{\link{rarefit}}, \code{\link{rarefit.predict}}
#'
#' @export
rarefit.cv <- function(fitObj, y, X, errtype = "mean-squared-error", nfolds = 5, ...) {
  n <- length(y)
  nlam <- length(fitObj$lambda)
  nalpha <- length(fitObj$alpha)
  errs <- array(NA, dim=c(nlam, nalpha, nfolds))

  # define error function
  errfun <- function(est, truth) colMeans((est - truth)^2)
  if (errtype == "mean-absolute-error") {
    errfun <- function(est, truth) colMeans(abs(est - truth))
  } else if (errtype != "mean-squared-error") {
    stop("The error function needs to be either mean squared error or mean absolute error.")
  }

  # make folds
  nn <- round(n / nfolds)
  sizes <- rep(nn, nfolds)
  sizes[nfolds] <- sizes[nfolds] + n - nn * nfolds
  b <- c(0, cumsum(sizes))
  set.seed(100) # set.seed for random number generator
  ii <- sample(n)
  folds <- list()
  for (i in seq(nfolds))
    folds[[i]] <- ii[seq(b[i] + 1, b[i + 1])]
  folds

  # Fit based on folds and compute error metric
  for (i in seq(nfolds)) {
    # fit model on all but the ith fold
    fit_cv <- rarefit(y = y[-folds[[i]]], X = X[-folds[[i]], ], A = fitObj$A, Q = fitObj$Q,
                      intercept = fitObj$intercept, lambda = fitObj$lambda, alpha = fitObj$alpha, ...)
    pred_te <- lapply(seq(nalpha), function(k) {
                                                  if (fitObj$intercept) {
                                                    X[folds[[i]], ] %*% fit_cv$beta[[k]] + rep(fit_cv$beta0[[k]], each = length(folds[[i]]))
                                                  } else {
                                                    X[folds[[i]], ] %*% fit_cv$beta[[k]]
                                                  }
                                                })
    for (k in seq(nalpha)) errs[, k, i] <- errfun(pred_te[[k]], y[folds[[i]]])
    cat("##########################\n")
    cat(sprintf("Finished model fits for fold[%s].\n", i))
    cat("##########################\n")
  }
  m <- apply(errs, c(1, 2), mean)
  se <- apply(errs, c(1, 2), stats::sd) / sqrt(nfolds)
  ibest <- which(m == min(m), arr.ind = TRUE)[1, , drop = FALSE]

  list (folds = folds, errs = errs, m = m, se = se, ibest = ibest,
        lambda.best = fitObj$lambda[ibest[1]], alpha.best = fitObj$alpha[ibest[2]])
}


#' Make predictions from a rarefit object and a rarefit.cv object
#'
#' The function makes predictions using a \code{rarefit} object at optimal
#' (\code{lambda}, \code{alpha}) chosen by \code{rarefit.cv}.
#'
#' @param fitObj Output of \code{rarefit}.
#' @param cvObj Output of \code{rarefit.cv}.
#' @param newx Matrix of new values for x at which predictions are made.
#'
#' @return Returns a sequence of predictions.
#'
#' @examples
#' \dontrun{
#' # See vignette for more details.
#' set.seed(100)
#' ts <- sample(1:length(data.rating), 400) # Train set indices
#' # Fit the model on train set
#' ourfit <- rarefit(y = data.rating[ts], X = data.dtm[ts, ], hc = data.hc, lam.min.ratio = 1e-6,
#'                   nlam = 20, nalpha = 10, rho = 0.01, eps1 = 1e-5, eps2 = 1e-5, maxite = 1e4)
#' # Cross validation
#' ourfit.cv <- rarefit.cv(ourfit, y = data.rating[ts], X = data.dtm[ts, ],
#'                         rho = 0.01, eps1 = 1e-5, eps2 = 1e-5, maxite = 1e4)
#' # Prediction on test set
#' pred <- rarefit.predict(ourfit, ourfit.cv, data.dtm[-ts, ])
#' pred.error <- mean((pred - data.rating[-ts])^2)
#' }
#'
#' @seealso \code{\link{rarefit}}, \code{\link{rarefit.cv}}
#'
#' @export
rarefit.predict <- function(fitObj, cvObj, newx) {
  ibest.lambda <- cvObj$ibest[1]
  ibest.alpha <- cvObj$ibest[2]
  if (fitObj$intercept) {
    as.vector(newx %*% fitObj$beta[[ibest.alpha]][, ibest.lambda] + fitObj$beta0[[ibest.alpha]][ibest.lambda])
  } else {
    as.vector(newx %*% fitObj$beta[[ibest.alpha]][, ibest.lambda])
  }
}
