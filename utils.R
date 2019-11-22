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

fit_o <- function(dag, data) {

	p <- graph::numNodes(dag)

	B <- diag(p)
	colnames(B) <- rownames(B) <- 1:p
	for (i in 2:p) {
		dag_rev <- graph::reverseEdgeDirections(dag)
		edge_list <- dag_rev@edgeL[[i]]$edges
		parents <- edge_list[edge_list < i] #sobra
		if (length(parents) > 0) {
			model <- lm(data[, i] ~ data[, parents])
			B[i, parents] <- - model$coefficients[-1]
		}
	}
	
	return(B)
}

fit_oinv <- function(dag, data) {
	p <- graph::numNodes(dag)

	B <- diag(p)
	colnames(B) <- rownames(B) <- 1:p
	for (i in 2:p) {
		dag_rev <- graph::reverseEdgeDirections(dag)
		parents <- dag_rev@edgeL[[i]]$edges
		for (j in 1:(i - 1)) {
			if (length(parents) == 0) { 
				B[i, j] <- 0
			} else if (parents[length(parents)] < j) {
				B[i, j] <- 0
			} else if (parents[length(parents)] == j) {
				model <- lm(data[, i] ~ data[, parents] - 1)
				B[i, j] <- model$coefficients[length(model$coefficients)]
			} else {
				model <- lm(data[, i] ~ data[, 1:j] - 1)
				B[i, j] <- model$coefficients[j]
			}
		}
	}
	
	return(B)
}

