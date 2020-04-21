source("exp_lib.R")

devtools::install_github("irenecrsn/covchol")

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

rothman_exp <- function(m, repetition, nodes, n, ntrain) {
	
	for (p in nodes) {
		Sigmatrue <- readRDS(file = paste0("rothman_exp/sigma1true_", p, ".rds"))
		X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
		L1 <- f_chol[[m]](ntrain, X = X)
		
		Sigmatrue <- readRDS(file = paste0("rothman_exp/sigma2true_", p, ".rds"))
		X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
		L2 <- f_chol[[m]](ntrain, X = X)
		
		Sigmatrue <- readRDS(file = paste0("rothman_exp/sigma3true_", p, ".rds"))
		X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
		L3 <- f_chol[[m]](ntrain, X = X)
		
		saveRDS(L1 %*% t(L1),
				file = paste0("rothman_exp/sigma1", m, "_", p, "_r", repetition, ".rds"))
		saveRDS(L2 %*% t(L2),
				file = paste0("rothman_exp/sigma2", m, "_", p, "_r", repetition, ".rds"))
		saveRDS(L3 %*% t(L3),
				file = paste0("rothman_exp/sigma3", m, "_", p, "_r", repetition, ".rds"))
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

l_exp <- function(m, repetition, nodes, n, ntrain) {
	for (p in nodes) {
		densities <- c(1/p, 2/p, 3/p)
		for (d in densities) {
			Ltrue <- readRDS(file = paste0("l_exp/ltrue_", p, "_", d, "_r", repetition, ".rds"))
			
			X <- MASS::mvrnorm(n, rep(0, p), Sigma = Ltrue %*% t(Ltrue))
			
			Lest <- f_chol[[m]](ntrain, X = X)
			
			saveRDS(Lest,
							file = paste0("l_exp/m_", p, "_", d, "_r", repetition, ".rds"))
		}
	}
}

execute <- function(r, ename, emethod, ...) {
	dir.create(ename, showWarnings = FALSE)

	emethod(repetition = r, ...)
}

