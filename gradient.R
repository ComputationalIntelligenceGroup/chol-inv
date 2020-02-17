devtools::install_github("irenecrsn/covchol")

execute_experiment <- function(p, r, ename, emethod, ...) {
	
	n_cores <- min(r, parallel::detectCores() - 2)
	cl <- parallel::makeCluster(n_cores, outfile = "")
	doParallel::registerDoParallel(cl)
	
	iter <- foreach::foreach(repetition = seq(r), .combine = rbind)
	foreach::"%dopar%"(iter, {
		dir.create(ename, showWarnings = FALSE)
		
		for (nnode in p) {
			densities <- c(1/nnode, 2/nnode, 3/nnode)
			for (d in densities) {
				result <- emethod(p = nnode, d = d, ...)
				saveRDS(result, file = paste0(ename, "/", nnode, "_", d, "_r", repetition, ".rds")) 
			}
		}
	})
	
	parallel::stopCluster(cl)
}

p <- c(30, 100, 200, 500, 1000)
r <- 200

##### Gradient experiment
sparse_chol <- function(p, d, n, ntrain) {

	Ltrue <- covchol::rlower(p = p, d = d)
	diag(Ltrue) <- runif(p, 0.1, 1)

	X <- MASS::mvrnorm(n, rep(0, p), Sigma = Ltrue %*% t(Ltrue))
	Covtest <- cov(X[(ntrain + 1): n, ])
	Ltest <- t(zapsmall(chol(Covtest, pivot = TRUE)))

	llpath <- covchol::cholpath(X = X[1:ntrain, ])
	frobs <- lapply(X = llpath, FUN = function(res) {
		norm(res$L - Ltest, type = "F")
	})

	Lest <- llpath[[which.min(frobs)]]$L

	return(list("ltrue" = Ltrue, "lest" = Lest))
}

band_chol <- function(p, d, n, ntrain) {
	Ltrue <- covchol::rlower(p = p, d = d)
	diag(Ltrue) <- runif(p, 0.1, 1)
	Sigmatrue <- Ltrue %*% t(Ltrue)
	
	X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
	
	Sigmaest <- PDSCE::band.chol.cv(x = X, n.tr = ntrain)$sigma
	Lest <- zapsmall(t(chol(Sigmaest, pivot = TRUE)))
		
	return(list("ltrue" = Ltrue, "lest" = Lest))
}

#execute_experiment(p = p, r = r, ename = "sparse_chol", emethod = sparse_chol, n = 200, ntrain = 100)
execute_experiment(p = p, r = r, ename = "band_chol", emethod = band_chol, n = 200, ntrain = 100)