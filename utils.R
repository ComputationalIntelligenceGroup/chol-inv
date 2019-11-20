#' Calculates the inverse of the Cholesky factor via the
#' regression coefficients formula
#' 
#' @param coeffs Lower-triangular matrix of regression coefficients
#' 
#' @return Matrix inverse of the Cholesky factor
chol_inv <- function(coeffs) {
	o <- diag(nrow(coeffs))
	dimnames(o) <- list(rownames(coeffs), colnames(coeffs))
	p <- ncol(o)

	for (i in 2:p) {
		o[i, i - 1] <- coeffs[i, i - 1]
	}

	for (j in 1:(p - 2)) {
		for (i in (j + 2):p) {
			o[i, j] <- coeffs[i, j] +
				coeffs[i, (j + 1):(i - 1)] %*%
					o[(j + 1):(i - 1), j]
		}
	}

	return(o)
}
