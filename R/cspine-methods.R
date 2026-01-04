#' Extract estimated coefficients from a 'cspine' object
#'
#' @param object An object of class [cspine].
#' @param ... Additional arguments (not used).
#' @return A list with components 'beta', an array of precision matrix
#'  components, and 'gamma', the estimate mean component.
#' @method coef cspine
#' @export
coef.cspine <- function(object, ...) {
  result <- list()
  result$beta <- object$beta
  result$gamma <- object$gamma
  return(result)
}

#' Predict from a 'cspine' object
#'
#' This function predicts the mean and precision matrix for a new vector
#'  of covariates.
#' @param object An object of class [cspine].
#' @param u A vector of covariates.
#' @param ... Additional arguments (not used).
#' @return A list with components 'precision', the predicted precision matrix,
#'  and 'mean', the predicted mean given the new covariates.
#' @method predict cspine
#' @export
predict.cspine <- function(object, u, ...) {
  q <- dim(object$gamma)[2]
  u <- as.numeric(u)
  if (length(u) != q) {
    stop("Expected covariate vector of length ", q, ".")
  }
  omega <- apply(object$beta, c(1, 2), \(b) b %*% c(1, u))
  diag(omega) <- 1 / object$sigma2
  omega <- Matrix::Diagonal(x = object$sigma2) %*% omega
  eps <- 1e-10 # ensure stable inversion for prediction
  mu <- NULL
  while (is.null(mu)) {
    mu <- tryCatch(
      solve(omega + eps * diag(nrow(omega)), object$gamma %*% u),
      error = function(e) {
        eps <- eps * 10
        NULL
      }
    )
  }
  return(list(precision = omega, mean = mu))
}
