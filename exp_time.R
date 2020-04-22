source("exp_lib.R")

nodes <- c(30, 100, 200, 500)
repetitions <- 10
n <- 200
ntrain <- 100
dirname <- "time_exp/"

dir.create(dirname, showWarnings = FALSE)

for (p in nodes) {
	densities <- c(1/p, 2/p, 3/p)
	for (d in densities) {
		for (r in 1:repetitions) {
			Sigmatrue <- gmat::chol_mh(N = 1, p = p, d = d)[, , 1]
			X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
			for (m in names(f_chol)) {
				t_begin <- Sys.time()
				f_chol[[m]](ntrain = ntrain, X = X)
				t_end <- Sys.time()
	
				time <- as.double(difftime(t_end, t_begin, unit = "secs"), unit = "secs")
				saveRDS(time, file = paste0(dirname, m, "_", p, "_", d,
											"_r", r, ".rds"))
			}
		}
	}
}

