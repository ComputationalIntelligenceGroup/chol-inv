#' Calculates the inverse of the Cholesky factor via the
#' regression coefficients formula
#' 
#' @param coeffs Lower-triangular matrix of regression coefficients
#' 
#' @return Matrix inverse of the Cholesky factor
chol_inv <- function(coeffs) {
	o <- diag(nrow(coeffs))
	dimnames(o) <- list(rownames(coeffs), colnames(coeffs))

	for (i in 2:ncol(o)) {
		o[i, i - 1] <- - coeffs[i, i - 1]
	}

	for (i in 3:nrow(o)) {
		for (j in 1:(i - 2)) {
			o[i, j] <- - coeffs[i, j] -
				coeffs[i, (j + 1):(i - 1)] %*%
					o[(j + 1):(i - 1), j]
		}
	}

	return(o)
}

get_norms <- function(U, D, dag) {
	L <- chol_inv(U)
	tryCatch(L_solve <- solve(U), error = function(e) { print(e) })
	tryCatch(L_true <- solve(diag(p) - pcalg::wgtMatrix(dag)),
					 error = function(e) { print(e) })
	
	S <- L %*% D %*% t(L)
	S_solve <- L_solve %*% D %*% t(L_solve)
	tryCatch(S_true <- pcalg::trueCov(dag), error = function(e) { print(e) })
	
	norms <- eval(formals(base::norm)$type)
	
	norm_est_mle <- sapply(norms, norm, x = L - L_solve)
	norm_est_true <- sapply(norms, norm, x = L - L_true)
	norm_mle_true <- sapply(norms, norm, x = L_solve - L_true)
	
	L_norms <- data.frame(norm_est_mle, norm_est_true, norm_mle_true)
	
	norm_est_mle <- sapply(norms, norm, x = S - S_solve)
	norm_est_true <- sapply(norms, norm, x = S - S_true)
	norm_mle_true <- sapply(norms, norm, x = S_solve - S_true)
	tryCatch(
		norm_sample_true <- sapply(norms, norm, x = cov(data) - S_true),
		error = function(e) {	print(e) }
	)
	
	S_norms <- data.frame(norm_est_mle, norm_est_true, norm_mle_true,
														norm_sample_true)
	
	return(list(L = L_norms, S = S_norms))
}
