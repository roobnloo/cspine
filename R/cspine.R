#' Run Covariate-adjusted Sparse Precision with Natural Estimation
#'
#' @param responses \eqn{n \times p} matrix of responses
#' @param covariates \eqn{n \times q} matrix of covariates
#' @param sglmixpath A path of sparse-group lasso mixing parameter with \eqn{0\leq \alpha \leq 1}. Default is \eqn{0.1, 0.2, \dotsc, 0.9}.
#' @param nlambda The number of lambda values to use for cross-validation - default is 100.
#' @param lam_max The maximum lambda considered. Automatically calculated if NULL.
#' @param lambda_factor The smallest value of lambda as a fraction of the maximum lambda.
#' @param maxit The maximum number of iterations. Default is \eqn{3\times 10^6}.
#' @param tol The convergence threshhold for optimization. Default is \eqn{10^{-10}}.
#' @param nfolds Number of folds for cross-validation. Default is 5.
#' @param verbose If TRUE, prints progress messages. Default is TRUE.
#' @param ncores Runs the nodewise regressions in parallel using that many cores. Default is 1.
#' @param adaptive Use adaptive weights when fitting nodewise regressions. Default is FALSE.
#' @importFrom Matrix colMeans colSums
#' @importFrom stats sd
#' @importFrom sparsegl sparsegl
#' @import parallel
#' @export
cspine <- function(responses, covariates, sglmixpath = seq(0.1, 0.9, 0.1), nlambda = 100,
                  lam_max = NULL, lambda_factor = 1e-4, maxit = 3e6, tol = 1e-8, nfolds = 5,
                  verbose = TRUE, ncores = 1, adaptive = FALSE) {
  stopifnot(
    is.matrix(responses), is.matrix(covariates),
    nrow(responses) == nrow(covariates),
    all(sglmixpath >= 0), all(sglmixpath <= 1)
  )

  p <- ncol(responses)
  q <- ncol(covariates)
  bveclength <- (p - 1) * (q + 1)

  nsglmix <- length(sglmixpath)
  lambda <- matrix(nrow = nlambda, ncol = p)
  beta <- matrix(nrow = p, ncol = bveclength)
  gamma <- matrix(nrow = p, ncol = q)
  cvm <- array(dim = c(nlambda, nsglmix, p))
  sigma2 <- numeric(p)
  mse <- numeric(p)
  cv_lambda_idx <- numeric(p)
  cv_alpha_idx <- numeric(p)

  cov_scale <- scale(covariates, scale = FALSE)
  intx <- intxmx(responses, covariates)
  intx_scale <- scale(intx, scale = FALSE)

  nodewise <- function(node) {
    y <- responses[, node] - mean(responses[, node])
    intx_scale_node <- intx_scale[, -(seq(0, q) * p + node)]
    nodereg <- cv_cspine_node(
      y, cbind(cov_scale, intx_scale_node), p, q, nlambda, lam_max, lambda_factor, sglmixpath,
      maxit, tol, nfolds
    )
    message(node, " ", appendLF = FALSE)

    return(list(
      gamma = nodereg$gamma,
      beta = nodereg$beta,
      sigma2 = nodereg$sigma2,
      lambda = nodereg$lambda,
      mse = nodereg$mse,
      cvm = nodereg$cvm,
      cv_lambda_idx = nodereg$cv_lambda_idx,
      cv_alpha_idx = nodereg$cv_alpha_idx
    ))
  }

  message("Running nodewise regressions...")

  if (ncores > 1) {
    reg_result <- parallel::mclapply(seq_len(p), nodewise, mc.cores = ncores)
  } else {
    reg_result <- lapply(seq_len(p), nodewise)
  }

  for (node in seq_len(p)) {
    gamma[node, ] <- reg_result[[node]]$gamma
    beta[node, ] <- reg_result[[node]]$beta
    sigma2[node] <- reg_result[[node]]$sigma2
    lambda[, node] <- reg_result[[node]]$lambda
    mse[node] <- reg_result[[node]]$mse
    cvm[, , node] <- reg_result[[node]]$cvm
    cv_lambda_idx[node] <- reg_result[[node]]$cv_lambda_idx
    cv_alpha_idx[node] <- reg_result[[node]]$cv_alpha_idx
  }

  message("\nFinished regressions.")

  bhat_tens <- array(0, dim = c(p, p, q + 1))
  bhat_symm <- array(0, dim = c(p, p, q + 1))
  for (i in seq_len(p)) {
    bhat_tens[i, -i, ] <- beta[i, ]
  }

  for (h in seq_len(q + 1)) {
    bhat_symm[, , h] <- symmetrize(-diag(1 / sigma2) %*% bhat_tens[, , h])
  }

  outlist <- list(
    gamma = gamma,
    beta = bhat_symm,
    beta_raw = bhat_tens,
    sigma2 = sigma2,
    mse = mse,
    lambda = lambda,
    alpha = sglmixpath,
    cvm = cvm,
    cv_lambda_idx = cv_lambda_idx,
    cv_gmix_idx = cv_alpha_idx
  )
  class(outlist) <- "cspine"

  return(outlist)
}
