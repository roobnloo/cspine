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
#' @importFrom Matrix colMeans colSums Diagonal norm t
#' @importFrom stats sd
#' @importFrom sparsegl sparsegl
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

  nodewise <- function(node) {
    # Extract raw data for this node
    y_raw <- responses[, node]
    x_raw <- responses[, -node, drop = FALSE]
    u_raw <- covariates

    # Compute means
    muy <- mean(y_raw)
    mux_j <- Matrix::colMeans(x_raw)
    muu <- Matrix::colMeans(u_raw)

    # Center variables
    y_centered <- y_raw - muy
    x_centered <- sweep(x_raw, 2, mux_j, "-")
    u_centered <- sweep(u_raw, 2, muu, "-")

    # Compute standard deviations
    sdy <- sqrt(sum(y_centered^2) / n)
    sdx_j <- sqrt(Matrix::colSums(x_centered^2) / n)
    sdu <- sqrt(Matrix::colSums(u_centered^2) / n)

    # Prevent division by zero
    if (sdy < 1e-9) sdy <- 1
    sdx_j[sdx_j < 1e-9] <- 1
    sdu[sdu < 1e-9] <- 1

    # Standardize main effects: tilde_u and tilde_x
    u_std <- u_centered %*% Matrix::Diagonal(x = 1 / sdu)
    x_std <- x_centered %*% Matrix::Diagonal(x = 1 / sdx_j)

    # Form interaction columns w_kh = (u_h - mean(u_h)) * (x_k - mean(x_k))
    # Store means and sds of w_kh for back-transformation
    w_means <- matrix(0, nrow = p - 1, ncol = q)
    w_sds <- matrix(0, nrow = p - 1, ncol = q)
    w_std_list <- list()

    for (h in seq_len(q)) {
      for (k in seq_len(p - 1)) {
        # Form w_kh (already centered since u and x are centered)
        w_kh <- u_centered[, h] * x_centered[, k]

        # Compute mean and sd of w_kh
        w_means[k, h] <- mean(w_kh)
        w_centered <- w_kh - w_means[k, h]
        w_sds[k, h] <- sqrt(sum(w_centered^2) / n)

        # Prevent division by zero
        if (w_sds[k, h] < 1e-9) w_sds[k, h] <- 1

        # Standardize w_kh
        w_std_list[[(h - 1) * (p - 1) + k]] <- w_centered / w_sds[k, h]
      }
    }

    # Combine all standardized interaction columns
    w_std <- do.call(cbind, w_std_list)

    # Form full design matrix: [tilde_u, tilde_x, tilde_w]
    design_std <- cbind(u_std, x_std, w_std)

    # Fit regression on standardized data
    nodereg <- cv_cspine_node(
      y_centered, design_std, p, q, nlambda, lam_max, lambda_factor, sglmixpath,
      maxit, tol, nfolds
    )

    # Extract standardized coefficients
    gamma_std <- nodereg$gamma
    beta_std <- nodereg$beta
    gamma0_std <- nodereg$gamma0

    # Back-transform to original scale
    # beta_kh = beta_kh^std / sd(w_kh)
    beta_unstd <- numeric(length(beta_std))

    # Main effect betas: beta_k0^std (first p-1 elements)
    # Interaction betas: beta_kh^std (remaining elements, organized as blocks of p-1)

    # Transform interaction coefficients first (we need these for main effects)
    beta_interact_unstd <- matrix(0, nrow = p - 1, ncol = q)
    for (h in seq_len(q)) {
      for (k in seq_len(p - 1)) {
        idx <- (p - 1) + (h - 1) * (p - 1) + k
        beta_interact_unstd[k, h] <- beta_std[idx] / w_sds[k, h]
      }
    }

    # Transform gamma_h: gamma_h = gamma_h^std / sd(u_h) - sum_k mean(x_k) * beta_kh^std / sd(w_kh)
    gamma_unstd <- numeric(q)
    for (h in seq_len(q)) {
      gamma_unstd[h] <- gamma_std[h] / sdu[h] - sum(mux_j * beta_interact_unstd[, h])
    }

    # Transform beta_k0: beta_k0 = beta_k0^std / sd(x_k) - sum_h mean(u_h) * beta_kh^std / sd(w_kh)
    beta_main_unstd <- numeric(p - 1)
    for (k in seq_len(p - 1)) {
      beta_main_unstd[k] <- beta_std[k] / sdx_j[k] - sum(muu * beta_interact_unstd[k, ])
    }

    # Transform intercept
    gamma0_unstd <- muy + gamma0_std
    gamma0_unstd <- gamma0_unstd - sum(muu * gamma_unstd)
    gamma0_unstd <- gamma0_unstd - sum(mux_j * beta_main_unstd)

    # Add the interaction correction term: sum_{h,k} (mean(u_h)*mean(x_k) - mean(w_kh)) * beta_kh^std / sd(w_kh)
    for (h in seq_len(q)) {
      for (k in seq_len(p - 1)) {
        gamma0_unstd <- gamma0_unstd + (muu[h] * mux_j[k] - w_means[k, h]) * beta_interact_unstd[k, h]
      }
    }

    # Reassemble beta vector: [beta_main, beta_interactions]
    beta_unstd[1:(p - 1)] <- beta_main_unstd
    for (h in seq_len(q)) {
      beta_unstd[(p - 1) + (h - 1) * (p - 1) + 1:(p - 1)] <- beta_interact_unstd[, h]
    }

    # Hard-threshold for numerical stability
    gamma_unstd[abs(gamma_unstd) < 1e-9] <- 0
    beta_unstd[abs(beta_unstd) < 1e-9] <- 0
    if (abs(gamma0_unstd) < 1e-9) {
      gamma0_unstd <- 0
    }

    message(node, " ", appendLF = FALSE)

    return(list(
      gamma0 = gamma0_unstd,
      gamma = gamma_unstd,
      beta = beta_unstd,
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

  sigma2_diag <- Matrix::Diagonal(x = 1 / sigma2)
  bhat_symm[, , 1] <- symmetrize(-sigma2_diag %*% bhat_tens[, , 1], "and")
  for (h in seq(2, q + 1)) {
    bhat_symm[, , h] <- symmetrize(-sigma2_diag %*% bhat_tens[, , h], symmetrize_rule)
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
