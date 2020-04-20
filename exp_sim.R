source("exp_lib.R")

devtools::install_github("irenecrsn/covchol")

####### Sigma exp: for comparison of sigmas instead of l factors

# Generate true matrices
sigma_exp_gen <- function(repetition, nodes) {
	for (p in nodes) {
		densities <- c(1/p, 2/p, 3/p)
		for (d in densities) {
			Sigmatrue <- gmat::chol_mh(N = 1, p = p, d = d)[, , 1]
			
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
			Sigmaband <- band_est(ntrain, X = X)

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
			
			Lsparse <- gradient_est(ntrain, X = X)
			Sigmasparse <- Lsparse %*% t(Lsparse)
			
			saveRDS(Sigmasparse,
							file = paste0("sigma_exp/sigmasparse_", p, "_", d, "_r", repetition, ".rds"))
		}
	}
}

# Likelihood estimate
sigma_exp_sparse_f <- function(repetition, nodes, n, ntrain) {
	for (p in nodes) {
		densities <- c(1/p, 2/p, 3/p)
		for (d in densities) {
			Sigmatrue <- readRDS(file = paste0("sigma_exp/sigmatrue_", p, "_", d, "_r", repetition, ".rds"))
			
			X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
			
			Lsparse <- gradient_est_f(ntrain, X = X)
			Sigmasparse <- Lsparse %*% t(Lsparse)
			message("d:", d, " | p:",p)
			saveRDS(Sigmasparse,
							file = paste0("sigma_exp/sigmasparse_f_", p, "_", d, "_r", repetition, ".rds"))
		}
	}
}

sigma_exp_sample <- function(repetition, nodes, n) {
	for (p in nodes) {
		densities <- c(1/p, 2/p, 3/p)
		for (d in densities) {
			Sigmatrue <- readRDS(file = paste0("sigma_exp/sigmatrue_", p, "_", d, "_r", repetition, ".rds"))
			
			X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
			
			saveRDS(cov(X),
							file = paste0("sigma_exp/sigmasample_", p, "_", d, "_r", repetition, ".rds"))
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

rothman_exp_sample <- function(repetition, nodes, n) {
	
	for (p in nodes) {
		Sigmatrue <- readRDS(file = paste0("rothman_exp/sigma1true_", p, ".rds"))
		X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
		Sigma1sample <- cov(X)
		
		Sigmatrue <- readRDS(file = paste0("rothman_exp/sigma2true_", p, ".rds"))
		X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
		Sigma2sample <- cov(X)
		
		Sigmatrue <- readRDS(file = paste0("rothman_exp/sigma3true_", p, ".rds"))
		X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
		Sigma3sample <- cov(X)
		
		saveRDS(Sigma1sample,
						file = paste0("rothman_exp/sigma1sample_", p, "_r", repetition, ".rds"))
		saveRDS(Sigma2sample,
						file = paste0("rothman_exp/sigma2sample_", p, "_r", repetition, ".rds"))
		saveRDS(Sigma3sample,
						file = paste0("rothman_exp/sigma3sample_", p, "_r", repetition, ".rds"))
	}
}

rothman_exp_band <- function(repetition, nodes, n, ntrain) {
	
	for (p in nodes) {
		Sigmatrue <- readRDS(file = paste0("rothman_exp/sigma1true_", p, ".rds"))
		X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
		Sigma1band <- band_est(ntrain, X = X)

		Sigmatrue <- readRDS(file = paste0("rothman_exp/sigma2true_", p, ".rds"))
		X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
		Sigma2band <- band_est(ntrain, X = X)

		Sigmatrue <- readRDS(file = paste0("rothman_exp/sigma3true_", p, ".rds"))
		X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
		Sigma3band <- band_est(ntrain, X = X)

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
		L1sparse <- gradient_est(ntrain, X = X)
		
		Sigmatrue <- readRDS(file = paste0("rothman_exp/sigma2true_", p, ".rds"))
		X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
		L2sparse <- gradient_est(ntrain, X = X)
		
		Sigmatrue <- readRDS(file = paste0("rothman_exp/sigma3true_", p, ".rds"))
		X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
		L3sparse <- gradient_est(ntrain, X = X)
		
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
			Ltrue <- gmat::mh_u(N = 1, p = p, dag = gmat::rgraph(p = p, d = d, dag = TRUE))[, , 1][p:1, p:1]

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
			
			Lest <- gradient_est(ntrain, X = X)
		        
                        message("d: ",d, " | p: ",p) 	
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
			
			Lest <- lasso_est(ntrain, X = X)
			
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
			
			Lest <- nestedlasso_est(ntrain, X = X)
			
			saveRDS(Lest,
							file = paste0("l_exp/nestedlasso_", p, "_", d, "_r", repetition, ".rds"))
		}
	}
}

execute <- function(r, ename, emethod, ...) {
	dir.create(ename, showWarnings = FALSE)

	emethod(repetition = r, ...)
}

