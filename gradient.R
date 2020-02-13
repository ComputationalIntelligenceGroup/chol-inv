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

p <- c(10, 30, 50, 100)
r <- 200

##### Gradient experiment
sparse_chol <- function(p, d) {
	n <- 20*p
	res <- list()
	
	Ltrue <- covchol::rlower(p = p, d = d)
	diag(Ltrue) <- runif(p, 0.1, 1)
	Sigmatrue <- Ltrue %*% t(Ltrue)

	X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
	train <- floor(n/2)
	Cortrue <- cov2cor(Sigmatrue)
	
	Cortrain <- cor(X[1:train, ])
	Cortest <- cor(X[(train + 1):n, ])

	llpath <- covchol::cholpath(Sigma = Cortrain)
	frobs <- lapply(X = llpath, FUN = function(res) {
		norm(res$L %*% t(res$L) - Cortest, type = "F")
	})

	Lest <- llpath[[which.min(frobs)]]$L
	Corest <- Lest %*% t(Lest)
	
	return(list("ltrue" = Ltrue, "lest" = Lest))
}

band_chol <- function(p, d) {
	n <- 2*p

	Ltrue <- covchol::rlower(p = p, d = d)
	diag(Ltrue) <- runif(p, 0.1, 1)
	Sigmatrue <- Ltrue %*% t(Ltrue)
	
	X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
	
	Sigmaest <- PDSCE::band.chol.cv(x = X)$sigma
	Lest <- zapsmall(t(chol(Sigmaest)))
		
	return(list("ltrue" = Ltrue, "lest" = Lest))
}

execute_experiment(p = p, r = r, ename = "sparse_chol", emethod = sparse_chol)
execute_experiment(p = p, r = r, ename = "band_chol", emethod = band_chol)