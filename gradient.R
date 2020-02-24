devtools::install_github("irenecrsn/covchol")

# Returns an estimate of the covariance matrix
band_est <- function(n, ntrain, X) {

	Sigmaest <- PDSCE::band.chol.cv(x = X, n.tr = ntrain)$sigma
	
	return(Sigmaest)
}

# Returns an estimate of the Cholesky factor
gradient_est <- function(n, ntrain, X) {

	Covtest <- cov(X[(ntrain + 1):n,])
	
	llpath <- covchol::cholpath(X = X[1:ntrain,])
	frobs <- lapply(
		X = llpath,
		FUN = function(res) {
			norm(res$L %*% t(res$L) - Covtest, type = "F")
		}
	)
	
	Lest <- llpath[[which.min(frobs)]]$L
	
	return(Lest)
}

# for comparison of sigmas instead of l factors
sigma_exp <- function(repetition, nodes, n, ntrain) {
	for (p in nodes) {
		densities <- c(1/p, 2/p, 3/p)
		for (d in densities) {
			Ltrue <- covchol::rlower(p = p, d = d)
			diag(Ltrue) <- runif(p, 0.1, 1)
			Sigmatrue <- Ltrue %*% t(Ltrue)
			
			X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
			
			Lsparse <- gradient_est(n, ntrain, X = X)
			Sigmaband <- band_est(n, ntrain, X = X)
			Sigmasparse <- Lsparse %*% t(Lsp)
			
			result <- list("sigmatrue" = Sigmatrue, 
									"sigmasparse" = Sigmasparse,
									"sigmaband" = Sigmaband)		
			saveRDS(result,
							file = paste0("sigma_exp/", p, "_", d, "_r", repetition, ".rds"))
		}
	}
}

# 
rothman_exp <- function(repetition, nodes, n, ntrain) {
	
	for (p in nodes) {
		## Sigma1 = AR1
		sigma1 <- matrix(ncol = p,
										 nrow = p,
										 data = 1)
		for (i in 2:p) {
			for (j in 1:(i - 1)) {
				sigma1[i, j] <- sigma1[j, i] <- 0.7 ^ (i - j)
			}
		}
		
		X <- MASS::mvrnorm(n, rep(0, p), Sigma = sigma1)
		Sigma1band <- band_est(n, ntrain, X = X)
		Lest1sparse <- gradient_est(n, ntrain, X = X)
		
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
		
		X <- MASS::mvrnorm(n, rep(0, p), Sigma = sigma2)
		Sigma2band <- band_est(n, ntrain, X = X)
		Lest2sparse <- gradient_est(n, ntrain, X = X)
		
		sigma3 <- matrix(ncol = p,
										 nrow = p,
										 data = 0.5)
		diag(sigma3) <- 1
		
		
		X <- MASS::mvrnorm(n, rep(0, p), Sigma = sigma3)
		Sigma3band <- band_est(n, ntrain, X = X)
		Lest3sparse <- gradient_est(n, ntrain, X = X)
		
		result1 <-	list(
				"sigmaband" = Sigma1band,
				"sigmasparse" = Lest1sparse %*% t(Lest1sparse),
				"sigmatrue" = sigma1
			)
		result2 <-	list(
			"sigmaband" = Sigma2band,
			"sigmasparse" = Lest2sparse %*% t(Lest2sparse),
			"sigmatrue" = sigma2
		)
		result3 <-	list(
			"sigmaband" = Sigma3band,
			"sigmasparse" = Lest3sparse %*% t(Lest3sparse),
			"sigmatrue" = sigma3
		)
		
		saveRDS(result1,
						file = paste0("rothman_exp/sigma1_", p, "_r", repetition, ".rds"))
		saveRDS(result2,
						file = paste0("rothman_exp/sigma2_", p, "_r", repetition, ".rds"))
		saveRDS(result3,
						file = paste0("rothman_exp/sigma3_", p, "_r", repetition, ".rds"))
	}
}

l_exp <- function(repetition, nodes, n, ntrain) {
	for (p in nodes) {
		densities <- c(1/p, 2/p, 3/p)
		for (d in densities) {
			Ltrue <- covchol::rlower(p = p, d = d)
			diag(Ltrue) <- runif(p, 0.1, 1)
			Sigmatrue <- Ltrue %*% t(Ltrue)
			
			X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
			
			Sigmasparse <- gradient_est(n, ntrain, X = X)
			Sigmaband <- band_est(n, ntrain, X = X)
			
			result <- list("sigmatrue" = Sigmatrue, 
										 "sigmasparse" = Sigmasparse,
										 "sigmaband" = Sigmaband)		
			saveRDS(result,
							file = paste0("sigma_exp/", p, "_", d, "_r", repetition, ".rds"))
		}
	}
}

execute_experiment <- function(r, ename, emethod, ...) {
	n_cores <- min(r, parallel::detectCores() - 2)
	cl <- parallel::makeCluster(n_cores, outfile = "")
	doParallel::registerDoParallel(cl)
	
	iter <- foreach::foreach(repetition = seq(r), .combine = rbind,
													 .export = c("band_est", "gradient_est"))
	foreach::"%dopar%"(iter, {
		dir.create(ename, showWarnings = FALSE)
		
		emethod(repetition = repetition, ...)

	})
	
	parallel::stopCluster(cl)
}

nodes <- c(30, 100, 200, 500, 1000)
#execute_experiment(r = 200, ename = "sigma_exp", emethod = sigma_exp, nodes = nodes, n = 200, ntrain = 100)
execute_experiment(r = 200, ename = "rothman_exp", emethod = rothman_exp, nodes = nodes, n = 200, ntrain = 100)
