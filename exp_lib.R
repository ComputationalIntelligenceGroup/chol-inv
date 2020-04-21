source("nestedlasso.r")
source("lasso.r")

# Returns an estimate of the Cholesky factor
gradient_est <- function(ntrain, X) {
	
	n <- nrow(X)

	#Covtest <- cov(X[(ntrain + 1):n,])
	datatest <-  scale(X[(ntrain+1):n,], center = TRUE, scale = FALSE) 
	llpath <- covchol::cholpath(X = X[1:ntrain,])
	frobs <- lapply(
		X = llpath,
		FUN = function(res) {
		  h <- forwardsolve(res$L, t(datatest))
			#norm(res$L %*% t(res$L) - Covtest, type = "F")
		  2*sum(log(diag(res$L))) + sum(crossprod(h))
		}
	)
	
	Lest <- llpath[[which.min(frobs)]]$L
	
	return(Lest)
}

# Returns an estimate of the Cholesky factor
gradient_est_f <- function(ntrain, X) {
	
	n <- nrow(X)

	Covtest <- cov(X[(ntrain + 1):n,])
	
	path <- covchol::cholpathf(X = X[1:ntrain,])
	frobs <- lapply(
		X = path,
		FUN = function(res) {
			norm(res$L %*% t(res$L) - Covtest, type = "F")
		}
	)
	
	Lest <- path[[which.min(frobs)]]$L
	
	return(Lest)
}

# Returns an estimate of the covariance matrix
band_est <- function(ntrain, X) {

	Sigmaest <- PDSCE::band.chol.cv(x = X, n.tr = ntrain, nsplits = 1)$sigma
	
	return(Sigmaest)
}

# Returns an estimate of the Cholesky factor
nestedlasso_est <- function(ntrain, X) {
	
	n <- nrow(X)
	
	Covtest <- cov(X[(ntrain + 1):n,])
	
	llpath <- nested.lasso.path(X = X[1:ntrain,])
	frobs <- lapply(
		X = llpath,
		FUN = function(res) {
			norm(res$sigma - Covtest, type = "F")
		}
	)
	
	selected <- llpath[[which.min(frobs)]]
	Lest <- selected$cholesky %*% diag(sqrt(selected$sigma2))
	
	return(Lest)
}

# Returns an estimate of the Cholesky factor
lasso_est <- function(ntrain, X) {
	n <- nrow(X)
	
	Covtest <- cov(X[(ntrain + 1):n,])
	
	llpath <- lasso.path(X = X[1:ntrain,])
	frobs <- lapply(
		X = llpath,
		FUN = function(res) {
			norm(res$sigma - Covtest, type = "F")
		}
	)
	
	selected <- llpath[[which.min(frobs)]]
	Lest <- selected$cholesky %*% diag(sqrt(selected$sigma2))
	
	return(Lest)
}
