source("nestedlasso.r")
source("lasso.r")
devtools::install_github("irenecrsn/covchol")

# Returns an estimate of the covariance matrix
band_est <- function(n, ntrain, X) {

	Sigmaest <- PDSCE::band.chol.cv(x = X, n.tr = ntrain)$sigma
	
	return(Sigmaest)
}

# Returns an estimate of the Cholesky factor
nestedlasso_est <- function(n, ntrain, X) {
	
	Covtest <- cov(X[(ntrain + 1):n,])
	
	llpath <- nested.lasso.path(X = X[1:ntrain,])
	frobs <- lapply(
		X = llpath,
		FUN = function(res) {
			norm(res$sigma - Covtest, type = "F")
		}
	)
	
	selected <- llpath[[which.min(frobs)]]
	Lest <- covchol::cholfromldl(L = selected$cholesky, D = selected$sigma2)
	
	return(Lest)
}

# Returns an estimate of the Cholesky factor
lasso_est <- function(n, ntrain, X) {
	
	Covtest <- cov(X[(ntrain + 1):n,])
	
	llpath <- lasso.path(X = X[1:ntrain,])
	frobs <- lapply(
		X = llpath,
		FUN = function(res) {
			norm(res$sigma - Covtest, type = "F")
		}
	)
	
	selected <- llpath[[which.min(frobs)]]
	Lest <- covchol::cholfromldl(L = selected$cholesky, D = selected$sigma2)
	
	return(Lest)
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

####### Sigma exp: for comparison of sigmas instead of l factors

# Generate true matrices
sigma_exp_gen <- function(repetition, nodes) {
	for (p in nodes) {
		densities <- c(1/p, 2/p, 3/p)
		for (d in densities) {
			Ltrue <- 
			Sigmatrue <- Ltrue %*% t(Ltrue)
			
			saveRDS(Sigmatrue,
							file = paste0("sigma_exp/sigmatrue_", p, "_", d, "_r", repetition, ".rds"))
		}
	}
}

# Banding estimate
sigma_exp_band <- function(repetition, nodes, n, ntrain) {
	for (p in nodes) {
		densities <- c(1/p, 2/p, 3/p)
		for (d in densities) {
			Sigmatrue <- readRDS(file = paste0("sigma_exp/sigmatrue_", p, "_", d, "_r", repetition, ".rds"))
			
			X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
			Sigmaband <- band_est(n, ntrain, X = X)

			saveRDS(Sigmaband,
							file = paste0("sigma_exp/sigmaband_", p, "_", d, "_r", repetition, ".rds"))
		}
	}
}

# Likelihood estimate
sigma_exp_sparse <- function(repetition, nodes, n, ntrain) {
	for (p in nodes) {
		densities <- c(1/p, 2/p, 3/p)
		for (d in densities) {
			Sigmatrue <- readRDS(file = paste0("sigma_exp/sigmatrue_", p, "_", d, "_r", repetition, ".rds"))
			
			X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
			
			Lsparse <- gradient_est(n, ntrain, X = X)
			Sigmasparse <- Lsparse %*% t(Lsparse)
			
			saveRDS(Sigmasparse,
							file = paste0("sigma_exp/sigmasparse_", p, "_", d, "_r", repetition, ".rds"))
		}
	}
}


####### Rothman exp: for comparison of fixed covariance matrices
rothman_exp_gen <- function(nodes) {
	
	dir.create("rothman_exp", showWarnings = FALSE)
	
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
		
		sigma3 <- matrix(ncol = p,
										 nrow = p,
										 data = 0.5)
		diag(sigma3) <- 1
		
		saveRDS(sigma1,
						file = paste0("rothman_exp/sigma1true_", p, ".rds"))
		saveRDS(sigma2,
						file = paste0("rothman_exp/sigma2true_", p, ".rds"))
		saveRDS(sigma3,
						file = paste0("rothman_exp/sigma3true_", p, ".rds"))
	}
}

rothman_exp_band <- function(repetition, nodes, n, ntrain) {
	
	for (p in nodes) {
		Sigmatrue <- readRDS(file = paste0("rothman_exp/sigma1true_", p, ".rds"))
		X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
		Sigma1band <- band_est(n, ntrain, X = X)

		Sigmatrue <- readRDS(file = paste0("rothman_exp/sigma2true_", p, ".rds"))
		X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
		Sigma2band <- band_est(n, ntrain, X = X)

		Sigmatrue <- readRDS(file = paste0("rothman_exp/sigma3true_", p, ".rds"))
		X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
		Sigma3band <- band_est(n, ntrain, X = X)

		saveRDS(Sigma1band,
						file = paste0("rothman_exp/sigma1band_", p, "_r", repetition, ".rds"))
		saveRDS(Sigma2band,
						file = paste0("rothman_exp/sigma2band_", p, "_r", repetition, ".rds"))
		saveRDS(Sigma3band,
						file = paste0("rothman_exp/sigma3band_", p, "_r", repetition, ".rds"))
	}
}

rothman_exp_sparse <- function(repetition, nodes, n, ntrain) {
	
	for (p in nodes) {
		Sigmatrue <- readRDS(file = paste0("rothman_exp/sigma1true_", p, ".rds"))
		X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
		L1sparse <- gradient_est(n, ntrain, X = X)
		
		Sigmatrue <- readRDS(file = paste0("rothman_exp/sigma2true_", p, ".rds"))
		X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
		L2sparse <- gradient_est(n, ntrain, X = X)
		
		Sigmatrue <- readRDS(file = paste0("rothman_exp/sigma3true_", p, ".rds"))
		X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
		L3sparse <- gradient_est(n, ntrain, X = X)
		
		saveRDS(L1sparse %*% t(L1sparse),
						file = paste0("rothman_exp/sigma1sparse_", p, "_r", repetition, ".rds"))
		saveRDS(L2sparse %*% t(L2sparse),
						file = paste0("rothman_exp/sigma2sparse_", p, "_r", repetition, ".rds"))
		saveRDS(L3sparse %*% t(L3sparse),
						file = paste0("rothman_exp/sigma3sparse_", p, "_r", repetition, ".rds"))
	}
}

####### L exp: for comparison of l factors
l_exp_gen <- function(repetition, nodes) {
	for (p in nodes) {
		densities <- c(1/p, 2/p, 3/p)
		for (d in densities) {
			Ltrue <- covchol::rlower(p = p, d = d)
			diag(Ltrue) <- runif(p, 0.1, 1)

			saveRDS(Ltrue,
							file = paste0("l_exp/ltrue_", p, "_", d, "_r", repetition, ".rds"))
		}
	}
}

l_exp_sparse <- function(repetition, nodes, n, ntrain) {
	for (p in nodes) {
		densities <- c(1/p, 2/p, 3/p)
		for (d in densities) {
			Ltrue <- readRDS(file = paste0("l_exp/ltrue_", p, "_", d, "_r", repetition, ".rds"))
			
			X <- MASS::mvrnorm(n, rep(0, p), Sigma = Ltrue %*% t(Ltrue))
			
			Lest <- gradient_est(n, ntrain, X = X)
			
			saveRDS(Lest,
							file = paste0("l_exp/sparse_", p, "_", d, "_r", repetition, ".rds"))
		}
	}
}

l_exp_lasso <- function(repetition, nodes, n, ntrain) {
	for (p in nodes) {
		densities <- c(1/p, 2/p, 3/p)
		for (d in densities) {
			Ltrue <- readRDS(file = paste0("l_exp/ltrue_", p, "_", d, "_r", repetition, ".rds"))
			
			X <- MASS::mvrnorm(n, rep(0, p), Sigma = Ltrue %*% t(Ltrue))
			
			Lest <- lasso_est(n, ntrain, X = X)
			
			saveRDS(Lest,
							file = paste0("l_exp/lasso_", p, "_", d, "_r", repetition, ".rds"))
		}
	}
}

l_exp_nestedlasso <- function(repetition, nodes, n, ntrain) {
	for (p in nodes) {
		densities <- c(1/p, 2/p, 3/p)
		for (d in densities) {
			Ltrue <- readRDS(file = paste0("l_exp/ltrue_", p, "_", d, "_r", repetition, ".rds"))
			
			X <- MASS::mvrnorm(n, rep(0, p), Sigma = Ltrue %*% t(Ltrue))
			
			Lest <- nestedlasso_est(n, ntrain, X = X)
			
			saveRDS(Lest,
							file = paste0("l_exp/nestedlasso_", p, "_", d, "_r", repetition, ".rds"))
		}
	}
}

execute_parallel <- function(r, ename, emethod, ...) {
	n_cores <- min(r, parallel::detectCores() - 2)
	cl <- parallel::makeCluster(n_cores, outfile = "")
	doParallel::registerDoParallel(cl)
	
	iter <- foreach::foreach(repetition = seq(r), .combine = rbind,
													 .export = c("band_est", "gradient_est", 
													 						"nestedlasso_est", "lasso_est",
													 						"nested.lasso.path", "autoReg", "nested.lasso.cov",
													 						"lasso.path", "lasso.cov"))
	foreach::"%dopar%"(iter, {
		dir.create(ename, showWarnings = FALSE)
		
		emethod(repetition = repetition, ...)

	})
	
	parallel::stopCluster(cl)
}

nodes <- c(30, 100, 200, 500, 1000)

#### Experiment over random covariance matrices
#execute_parallel(r = 200, ename = "sigma_exp", emethod = sigma_exp_gen, nodes = nodes)
#execute_parallel(r = 200, ename = "sigma_exp", emethod = sigma_exp_sparse, nodes = nodes, n = 200, ntrain = 100)
#execute_parallel(r = 200, ename = "sigma_exp", emethod = sigma_exp_band, nodes = nodes, n = 200, ntrain = 100)

#### Experiment over fixed covariance matrices
#rothman_exp_gen(nodes = nodes)
#execute_parallel(r = 200, ename = "rothman_exp", emethod = rothman_exp_sparse, nodes = nodes, n = 200, ntrain = 100)
#execute_parallel(r = 200, ename = "rothman_exp", emethod = rothman_exp_band, nodes = nodes, n = 200, ntrain = 100)

#### Experiment over L factors
#execute_parallel(r = 200, ename = "l_exp", emethod = l_exp_gen, nodes = nodes)
#execute_parallel(r = 200, ename = "l_exp", emethod = l_exp_sparse, nodes = nodes, n = 200, ntrain = 100)
#nodes <- c(30, 100, 200)
#execute_parallel(r = 200, ename = "l_exp", emethod = l_exp_lasso, nodes = nodes, n = 200, ntrain = 100)
#execute_parallel(r = 200, ename = "l_exp", emethod = l_exp_nestedlasso, nodes = nodes, n = 200, ntrain = 100)

