###### Adapted from the code from Rothman and Liu for nested lasso
# Authors: Irene Cordoba and Gherardo Varando
lasso.path <- function(X){
	X_centered <- scale(X, center = TRUE, scale = FALSE)
	
	lambdas = seq(0, 1, length = 10)
	
	results <- list()
	for (i in 1:length(lambdas)){
		results[[i]] <- lasso.cov(x = X_centered, lam = lambdas[i])
	}
	return(results)
}


## The main function is lasso.cov
## Arguments:
##     x is an n by p data matrix that is column-centered.
##     lam is the non-negative tuning parmater in the penalty

## the function returns a list with elements
##    cholesky is the estimate of L
##    sigma2 is the vector for which diag(sigma2) is D
##    sigma is the covarinace estimate LDL'


lasso.cov <- function(x, lam) 
{
	### x is a column-centered n by p matrix
	p <- dim(x)[2]
	cholesky <- diag(p)
	resid = x[,1,drop=FALSE]
	sigma2 = mean(resid^2)
	
	# First iteration
	newx <- cbind(resid, 0)
	newy <- x[, 2]
	tmp <- glmnet::glmnet(x = newx, y = newy, lambda = lam, standardize = FALSE,
						intercept = FALSE)  
	cholesky[2, 1] <- tmp$beta[1, 1]
	tmp_resid <- newy - resid%*%tmp$beta[1, 1]
	resid <- cbind(resid, tmp_resid)
	sigma2 <- c(sigma2, mean(tmp_resid^2))
	for (j in 3:p) 
	{
		newy <- x[, j]
		newx <- resid
		tmp <- glmnet::glmnet(x = newx, y = newy, lambda = lam, standardize = FALSE,
													intercept = FALSE)  
		cholesky[j, 1:(j-1)] <- tmp$beta[,1]
		
		tmp_resid <- newy - newx%*%tmp$beta[,1]
		resid <- cbind(resid, tmp_resid)
		sigma2 <- c(sigma2, mean(tmp_resid^2))
	}
	sigma <- tcrossprod(cholesky %*% diag(sigma2),cholesky)
	return(list(cholesky=cholesky, sigma2=sigma2, sigma=sigma) )
	
}
