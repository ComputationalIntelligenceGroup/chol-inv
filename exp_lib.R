source("lasso.r")

gradient_est <- function(ntrain, X) {
	
	n <- nrow(X)

	datatest <-  scale(X[(ntrain+1):n,], center = TRUE, scale = FALSE) 
	
	llpath <- covchol::cholpath(X = X[1:ntrain,])
	lls <- lapply(
		X = llpath,
		FUN = function(res) {
		  h <- forwardsolve(res$L, t(datatest))
		  sum(log(diag(res$L))) + sum(crossprod(h))
		}
	)
	
	Lest <- llpath[[which.min(lls)]]$L
	
	return(Lest)
}

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

band_est <- function(ntrain, X) {

	Sigmaest <- PDSCE::band.chol.cv(x = X, n.tr = ntrain, nsplits = 1)$sigma
	
	return(chol(Sigmaest))
}

lasso_est <- function(ntrain, X) {
	n <- nrow(X)
	
	datatest <-  scale(X[(ntrain+1):n,], center = TRUE, scale = FALSE) 
	
	llpath <- lasso.path(X = X[1:ntrain,])
	lls <- lapply(
		X = llpath,
		FUN = function(res) {
		  L <- res$cholesky %*% diag(sqrt(res$sigma2))
		  h <- forwardsolve(L, t(datatest))
		  sum(log(diag(L))) + sum(crossprod(h))
		}
	)
	
	selected <- llpath[[which.min(lls)]]
	Lest <- selected$cholesky %*% diag(sqrt(selected$sigma2))
	
	return(Lest)
}

f_chol <- c("grad_lik" = gradient_est,
			"grad_frob" = gradient_est_f,
			"band" = band_est,
			"lasso" = lasso_est)

