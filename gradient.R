#devtools::install_github("irenecrsn/covchol")

p <- c(30, 100, 200, 500, 1000)
r <- 200

# Returns an estimate of the covariance matrix
band_est <- function(n, ntrain, Sigmatrue) {
	p <- ncol(Sigmatrue)
	X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
	
	Sigmaest <- PDSCE::band.chol.cv(x = X, n.tr = ntrain)$sigma
	
	return(Sigmaest)
}

# Returns an estimate of the Cholesky factor
gradient_est <- function(n, ntrain, Sigmatrue) {
	p <- ncol(Sigmatrue)
	X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
	Covtest <- cov(X[(ntrain + 1):n,])
	
	llpath <- covchol::cholpath(X = X[1:ntrain,])
	frobs <- lapply(
		X = llpath,
		FUN = function(res) {
			norm(res$L %*% t(res$L) - Covtest, type = "F")
		}
	)
	
	Lest <- llpath[[which.min(frobs)]]$L
	
	return(Lest %*% t(Lest))
}

# for comparison of sigmas instead of l factors
sigma_exp <- function(p, d, n, ntrain) {
	Ltrue <- covchol::rlower(p = p, d = d)
	diag(Ltrue) <- runif(p, 0.1, 1)
	Sigmatrue <- Ltrue %*% t(Ltrue)
	
	Sigmasparse <- gradient_est(n, ntrain, Sigmatrue = Sigmatrue)
	Sigmaband <- band_est(n, ntrain, Sigmatrue = Sigmatrue)
	
	return(list("sigmatrue" = Sigmatrue, 
							"sigmasparse" = Sigmasparse,
							"sigmaband" = Sigmaband))
}

# d is ignored in this experiment
rothman_exp <- function(p, d, n, ntrain) {
	## Sigma1 = AR1
	sigma1 <- matrix(ncol = p,
									 nrow = p,
									 data = 1)
	for (i in 2:p) {
		for (j in 1:(i - 1)) {
			sigma1[i, j] <- sigma1[j, i] <- 0.7 ^ (i - j)
		}
	}
	
	Lest1band <- band_est(n, ntrain, Sigmatrue = sigma1)
	Lest1sparse <- gradient_est(n, ntrain, Sigmatrue = sigma1)
	
	sigma2 <- matrix(ncol = p,
									 nrow = p,
									 data = 0)
	for (i in 1:p) {
		if (i <= p - 1) {
			sigma2[i, i + 1] <- sigma2[i + 1, i] <- 0.4
			if (i <= p - 2) {
				sigma2[i, i + 2] <- sigma2[i + 2, i] <- 0.2
				if (i <= p - 3) {
					sigma2[i, i + 3] <- sigma2[i + 3, i] <- 0.2
					if (i <= p - 4) {
						sigma2[i, i + 4] <- sigma2[i + 4, i] <- 0.1
					}
				}
			}
		}
	}
	diag(sigma2) <- 1
	
	Lest2band <- band_est(n, ntrain, Sigmatrue = sigma2)
	Lest2sparse <- gradient_est(n, ntrain, Sigmatrue = sigma2)
	
	sigma3 <- matrix(ncol = p,
									 nrow = p,
									 data = 0.5)
	diag(sigma3) <- 1
	
	Lest3band <- band_est(n, ntrain, Sigmatrue = sigma3)
	Lest3sparse <- gradient_est(n, ntrain, Sigmatrue = sigma3)
	
	return(
		list(
			"Lest1band" = Lest1band,
			"Lest2band" = Lest2band,
			"Lest3band" = Lest3band,
			"Lest1sparse" = Lest1sparse,
			"Lest2sparse" = Lest2sparse,
			"Lest3sparse" = Lest3sparse
		)
	)
}

execute_experiment <- function(p, r, ename, emethod, ...) {
	n_cores <- min(r, parallel::detectCores() - 2)
	cl <- parallel::makeCluster(n_cores, outfile = "")
	doParallel::registerDoParallel(cl)
	
	iter <- foreach::foreach(repetition = seq(r), .combine = rbind,
													 .export = c("band_est", "gradient_est"))
	foreach::"%dopar%"(iter, {
		dir.create(ename, showWarnings = FALSE)
		
		for (nnode in p) {
			densities <- c(1 / nnode, 2 / nnode, 3 / nnode)
			for (d in densities) {
				result <- emethod(p = nnode, d = d, ...)
				saveRDS(result,
								file = paste0(ename, "/", nnode, "_", d, "_r", repetition, ".rds"))
			}
		}
	})
	
	parallel::stopCluster(cl)
}

execute_experiment(p = p, r = r, ename = "sigma_exp", emethod = sigma_exp, n = 200, ntrain = 100)
#execute_experiment(p = p, r = r, ename = "rothman_exp", emethod = rothman_exp, n = 200, ntrain = 100)
