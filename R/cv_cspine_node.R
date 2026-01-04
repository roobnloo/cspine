cv_cspine_node <- function(y, uw, p, q, nlambda, lam_max, lambda_factor, alpha,
                           maxit, tol, nfolds) {
  n <- length(y)
  nvars <- q + (p - 1) * (q + 1)
  nalpha <- length(alpha)

  muy <- mean(y)
  muuw <- Matrix::colMeans(uw)

  # This centering allows us to not include an intercept in the regression.
  yc <- y - muy
  uwc <- sweep(uw, 2, muuw, "-")

  # Convert to format accepted by sparsegl (base matrix or sparseMatrix)
  if (!inherits(uwc, "sparseMatrix") && !is.matrix(uwc)) {
    uwc <- as.matrix(uwc)
  }

  foldid <- cut(sample(seq_len(n)), nfolds, labels = FALSE)
  groupid <- c(rep(0, q), rep(1:(q + 1), each = p - 1)) + 1
  cvm_mx <- matrix(0, nrow = nlambda, ncol = nalpha)
  coefs <- matrix(nrow = nvars, ncol = nalpha)
  lambda <- numeric(nlambda)
  mse <- numeric(nalpha)
  if (is.null(lam_max)) {
    amin <- min(alpha)
    if (amin == 0) {
      amin <- 1
    }
    lam_max <- Matrix::norm(crossprod(uwc, yc), type = "I") / (n * amin)
  }
  lambda <- lam_max * exp(seq(log(1), log(lambda_factor), length = nlambda))

  for (asid in seq_along(alpha)) {
    asparse <- alpha[asid]
    pf_group <- c(0, 0, rep(sqrt(p - 1), q))
    sgl1 <- sparsegl::cv.sparsegl(
      uwc, yc,
      group = groupid,
      foldid = foldid,
      lambda = lambda,
      pf_group = pf_group,
      asparse = asparse,
      eps = tol, maxit = maxit,
      intercept = FALSE,
      standardize = FALSE
    )
    cvm_mx[, asid] <- sgl1$cvm
    lambda_min_ind <- which.min(sgl1$cvm)
    fit <- sgl1$sparsegl.fit
    coefs[, asid] <- as.numeric(fit$beta[, lambda_min_ind])
    mse[asid] <- fit$mse[lambda_min_ind]
  }

  cv_ind <- arrayInd(which.min(cvm_mx), dim(cvm_mx))
  alpha_min_ind <- cv_ind[2]
  gamma <- coefs[1:q, alpha_min_ind]
  beta <- coefs[(q + 1):nvars, alpha_min_ind]

  # Compute nnz based on largest magnitude coefficients
  bcs <- cumsum(sort(abs(beta), decreasing = TRUE))
  nnz <- which(bcs >= 0.999 * sum(abs(beta)))[1] + 1
  if (nnz >= n) {
    sigma2 <- 1
  } else {
    sigma2 <- mse[alpha_min_ind] * n / abs(n - nnz)
  }

  return(list(
    gamma0 = muy - muuw %*% c(gamma, beta),
    gamma = gamma,
    beta = beta,
    sigma2 = sigma2,
    lambda = lambda,
    mse = mse[alpha_min_ind],
    cvm = cvm_mx,
    cv_idx = cv_ind
  ))
}
