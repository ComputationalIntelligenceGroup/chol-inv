# Authors: Irene Cordoba and Gherardo Varando
nested.lasso.path <- function(X){
	X_centered <- scale(X, center = TRUE, scale = FALSE)

	lambdas = seq(0, 1, length = 10)
	
	results <- list()
	for (i in 1:length(lambdas)){
		results[[i]] <- nested.lasso.cov(x = X_centered, lam = lambdas[i])
	}
	return(results)
}


## The main function is nest.lasso.cov
## Arguments:
##     x is an n by p data matrix that is column-centered.
##     lam is the non-negative tuning parmater in the penalty

## the function returns a list with elements
##    cholesky is the estimate of L
##    sigma2 is the vector for which diag(sigma2) is D
##    sigma is the covarinace estimate LDL'


nested.lasso.cov <- function(x, lam) 
{
  ### x is a column-centered n by p matrix
  p <- dim(x)[2]
  cholesky <- diag(p)
  resid = x[,1,drop=FALSE]
  sigma2 = mean(resid^2)
  for (j in 2:p) 
  {
    newy <- x[, j]
    newx <- resid
    tmp <- autoReg(x=newx, y=newy, lambda1=lam)  
    cholesky[j, 1:(j-1)] <- tmp$phi
    resid = cbind(resid, newy - newx%*%tmp$phi)
    sigma2 <- c(sigma2, tmp$sigma2)
  }
  sigma <- tcrossprod(cholesky %*% diag(sigma2),cholesky)
  return(list(cholesky=cholesky, sigma2=sigma2, sigma=sigma) )

}

autoReg <- function(x, y, lambda1, lambda2, eps1=10^(-10), eps2=10^(-4)) 
{
  ## Written by Ji Zhu and Adam Rothman
  ### The penalty is \lambda (|\phi_{j,j-1}|/|\phi_{j,j-1}^0| + |\phi_{j,j-2}|/|\phi_{j,j-1}| + ... + |\phi_{j,1}|/|\phi_{j,2}|)
 
  xtx <- crossprod(x)
  xty <- crossprod(x,y)
  ##Initial estimate
  phi.old <- xty / apply(x^2, 2, sum)
  phi.j0 <- phi.old[length(phi.old)] ##

  diff.p <- 1
  kp <- 1
  while (diff.p > eps2 && kp <= 200) 
  {
    Dphi <- diag( length(phi.old) )
    phi.old.s <- c(phi.old[-1], phi.j0)
    diag(Dphi) <- lambda1 * as.vector(1/abs(phi.old * phi.old.s))
    phi.new <- qr.solve(xtx + Dphi, xty, tol=1e-200)
    if ( any(abs(phi.new) <= eps1) ) 
    {
      tmpind <- rev(seq(phi.new)[abs(phi.new) <= eps1])[1]
      phi.new[1:tmpind] <- eps1
    }
    diff.p <- sum(abs(phi.new - phi.old)) / sum(abs(phi.old))
    phi.old <- phi.new
    kp <- kp + 1
  }
  phi.new[ abs(phi.new) <= eps1 ] <- 0
  sigma2 <- mean((y - x %*% phi.new)^2)
  return( list(phi=phi.new, sigma2=sigma2) )
}
