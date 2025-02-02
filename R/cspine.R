#' Run Covariate-adjusted Sparse Precision with Natural Estimation
#'
#' @param responses \eqn{n \times p} matrix of responses
#' @param covariates \eqn{n \times q} matrix of covariates
#' @param sglmixpath A path of sparse-group lasso mixing parameter with \eqn{0 < \alpha \leq 1}.
#' @param nlambda The number of lambda values to use for cross-validation.
#' @param lam_max The maximum lambda considered. Automatically calculated if NULL.
#' @param lambda_factor The smallest value of lambda as a fraction of the maximum lambda.
#' @param symmetrize_rule Which rule to use to symmetrize the precision matrix components.
#' @param maxit The maximum number of iterations.
#' @param tol The convergence threshhold for optimization.
#' @param nfolds Number of folds for cross-validation.
#' @param ncores Runs the nodewise regressions in parallel using specified number cores. Defaults to no parallelization.
#' @param adaptive Use adaptive weights when fitting nodewise regressions.
#' @importFrom Matrix colMeans colSums
#' @importFrom stats sd
#' @importFrom sglssnal sglssnal
#' @import parallel
#' @export
cspine <- function(responses, covariates, sglmixpath = seq(0.1, 1, 0.1), nlambda = 100,
                   lam_max = NULL, lambda_factor = 1e-2, symmetrize_rule = c("and", "or"),
                   maxit = 3e6, tol = 1e-8, nfolds = 5,
                   ncores = 1, adaptive = FALSE) {
  stopifnot(
    is.matrix(responses), is.matrix(covariates),
    nrow(responses) == nrow(covariates),
    all(sglmixpath > 0), all(sglmixpath <= 1)
  )
  symmetrize_rule <- match.arg(symmetrize_rule)

  p <- ncol(responses)
  q <- ncol(covariates)
  n <- nrow(responses)
  bveclength <- (p - 1) * (q + 1)

  nsglmix <- length(sglmixpath)
  lambda <- matrix(nrow = nlambda, ncol = p)
  beta <- matrix(nrow = p, ncol = bveclength)
  gamma <- matrix(nrow = p, ncol = q)
  cvm <- array(dim = c(nlambda, nsglmix, p))
  cv_idx <- matrix(nrow = p, ncol = 2)
  sigma2 <- numeric(p)
  mse <- numeric(p)

  sdu <- sqrt(Matrix::colSums(covariates^2) / n)
  sdx <- sqrt(Matrix::colSums(responses^2) / n)
  muu <- Matrix::colMeans(covariates)
  mux <- Matrix::colMeans(responses)

  uc <- sweep(covariates, 2, muu, "-")
  uc <- sweep(uc, 2, sdu, "/")
  ux <- sweep(responses, 2, mux, "-")
  ux <- sweep(ux, 2, sdx, "/")

  nodewise <- function(node) {
    y <- responses[, node]
    intx_node <- intxmx(ux[, -node], uc)
    nodereg <- cv_cspine_node(
      y, cbind(uc, intx_node), p, q, nlambda, lam_max, lambda_factor, sglmixpath,
      maxit, tol, nfolds
    )
    muxj <- mux[-node]
    sdxj <- sdx[-node]
    gamma <- nodereg$gamma
    beta <- nodereg$beta
    temp <- sapply(
      1:q, \(h) sum(beta[seq((p - 1) * h + 1, (p - 1) * (h + 1))] * muxj / sdxj)
    )
    gamma_tilde <- (gamma - temp) / sdu

    # for fixed k, index over blocks h = 1 to q
    temp <- sapply(
      1:(p - 1), \(k) sum(beta[p - 1 + k + seq(0, q - 1) * (p - 1)] * muu / sdu)
    )
    beta_tilde <- beta
    beta_tilde[1:(p - 1)] <- (beta_tilde[1:(p - 1)] - temp) / sdxj
    for (h in seq_len(q)) {
      beta_tilde[(p - 1) * h + 1:(p - 1)] <- beta_tilde[(p - 1) * h + 1:(p - 1)] / (sdu[h] * sdxj)
    }

    gamma0_tilde <- mux[node] + nodereg$gamma0 - sum(gamma * muu / sdu) - sum(beta[1:(p - 1)] * muxj / sdxj)
    gamma0_tilde <- gamma0_tilde +
      sum(Reduce(c, lapply(1:q, \(h) muu[h] * muxj / (sdu[h] * sdxj))) * beta[-(1:(p - 1))])

    # hard-threshold for numerical stability
    gamma_tilde[abs(gamma_tilde) < 1e-9] <- 0
    beta_tilde[abs(beta_tilde) < 1e-9] <- 0
    if (abs(gamma0_tilde) < 1e-9) {
      gamma0_tilde <- 0
    }

    message(node, " ", appendLF = FALSE)

    return(list(
      gamma0 = gamma0_tilde,
      gamma = gamma_tilde,
      beta = beta_tilde,
      sigma2 = nodereg$sigma2,
      lambda = nodereg$lambda,
      mse = nodereg$mse,
      cvm = nodereg$cvm,
      cv_idx = nodereg$cv_idx
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
    cv_idx[node, ] <- reg_result[[node]]$cv_idx
  }

  message("\nFinished regressions.")

  bhat_tens <- array(0, dim = c(p, p, q + 1))
  bhat_symm <- array(0, dim = c(p, p, q + 1))
  for (i in seq_len(p)) {
    bhat_tens[i, -i, ] <- beta[i, ]
  }

  bhat_symm[, , 1] <- symmetrize(-diag(1 / sigma2) %*% bhat_tens[, , 1], "and")
  for (h in seq(2, q + 1)) {
    bhat_symm[, , h] <- symmetrize(-diag(1 / sigma2) %*% bhat_tens[, , h], symmetrize_rule)
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
    cv_idx = cv_idx
  )
  class(outlist) <- "cspine"

  return(outlist)
}
