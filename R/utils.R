#' @return symmetrized version of matrix mx
#' result_ij = result_ji is nonzero iff both mx_ij and mx_ji are nonzero,
#' in which case we choose the smaller value in magnitude.
#' @noRd
symmetrize <- function(mx, rule) {
  if (rule == "and") {
    result <- mx * (abs(mx) < Matrix::t(abs(mx))) + Matrix::t(mx) * (Matrix::t(abs(mx)) < abs(mx))
  } else {
    result <- mx * (abs(mx) >= Matrix::t(abs(mx))) + Matrix::t(mx) * (Matrix::t(abs(mx)) >= abs(mx))
  }
  return(as.matrix(result))
}
