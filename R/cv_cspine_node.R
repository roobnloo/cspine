cv_cspine_node <- function(y, uw, p, q, nlambda, lam_max, lambda_factor, alpha,
                           maxit, tol, nfolds) {
  n <- length(y)
  nvars <- q + (p - 1) * (q + 1)
  muy <- mean(y)
  muuw <- Matrix::colMeans(uw)

  # This centering allows us to not include an intercept in the regression.
  yc <- y - muy
  uwc <- sweep(uw, 2, muuw, "-")

  foldid <- cut(sample(seq_len(n)), nfolds, labels = FALSE)
  start_ids <- c(1, q + seq(0, q) * (p - 1) + 1)
  end_ids <- q + seq(0, q + 1) * (p - 1)
  grp_vec <- seq(1, nvars)
  grp_idx <- rbind(start_ids, end_ids)
  lambda <- numeric(nlambda)
  if (is.null(lam_max)) {
    lam_max <- norm(crossprod(uwc, yc), type = "I")
  }
  lambda <- lam_max * exp(seq(log(1), log(lambda_factor), length = nlambda))
  pf_group <- c(0, 0, rep(1, q))
  sgl1 <- sglssnal::cv.sglssnal(
    uwc, yc,
    grp_vec, grp_idx,
    lambdas = lambda,
    alphas = alpha,
    foldid = foldid,
    pfgroup = pf_group,
    quietall = TRUE,
    stoptol = tol
  )

  gamma <- sgl1$x[1:q]
  gamma0 <- muy - muuw %*% sgl1$x # intercept
  beta <- sgl1$x[(q + 1):nvars]

  nnz <- sgl1$info$nnz
  sigma2 <- 1
  if (n > nnz) {
    sigma2 <- sgl1$info$mse * n / abs(n - nnz)
  }

  return(list(
    gamma0 = gamma0,
    gamma = gamma,
    beta = beta,
    sigma2 = sigma2,
    lambda = lambda,
    mse = sgl1$info$mse,
    cvm = sgl1$cv_info$cvm,
    cv_idx = sgl1$cv_info$cv_idx
  ))
}
