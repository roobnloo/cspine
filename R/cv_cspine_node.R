cv_cspine_node <- function(y, uw, p, q, nlambda, lam_max, lambda_factor, alpha,
                           maxit, tol, nfolds) {
  n <- length(y)
  nvars <- q + (p - 1) * (q + 1)

  foldid <- cut(sample(seq_len(n)), nfolds, labels = FALSE)
  start_ids <- c(1, q + seq(0, q) * (p - 1) + 1)
  end_ids <- q + seq(0, q + 1) * (p - 1)
  grp_vec <- seq(1, nvars)
  grp_idx <- rbind(start_ids, end_ids)
  lambda <- numeric(nlambda)
  if (is.null(lam_max)) {
    lam_max <- norm(crossprod(uw, y), type = "I")
  }
  lambda <- lam_max * exp(seq(log(1), lambda_factor, length = nlambda))
  pf_group <- c(0, 0, rep(1, q))
  sgl1 <- sglssnal::cv.sglssnal(
    uw, y,
    grp_vec, grp_idx,
    lambdas = lambda,
    alphas = alpha,
    foldid = foldid,
    pfgroup = pf_group,
    quietall = TRUE,
    stoptol = tol
  )

  # for (asid in seq_along(alpha)) {
  #   asparse <- alpha[asid]
  #   cvm_mx[, asid] <- sgl1$cvm
  #   lambda_min_ind <- which.min(sgl1$cvm)
  #   fit <- sgl1$sparsegl.fit
  #   coefs[, asid] <- as.numeric(fit$beta[, lambda_min_ind])
  #   mse[asid] <- fit$mse[lambda_min_ind]
  # }

  # cv_ind <- arrayInd(which.min(cvm_mx), dim(cvm_mx))
  # alpha_min_ind <- cv_ind[2]
  gamma <- sgl1$x[1:q]
  beta <- sgl1$x[(q + 1):nvars]

  nnz <- sgl1$info$nnz
  sigma2 <- 1
  if (n > nnz) {
    sigma2 <- sgl1$info$mse * n / abs(n - nnz)
  }

  return(list(
    gamma = gamma,
    beta = beta,
    sigma2 = sigma2,
    lambda = lambda,
    mse = sgl1$info$mse,
    cvm = sgl1$cv_info$cvm,
    cv_lambda_idx = sgl1$cv_info[1],
    cv_alpha_idx = sgl1$cv_info[2]
  ))
}
