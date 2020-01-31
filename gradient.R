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
r <- 50

#p <- 5; d <- 0.5; n <- 100

##### Gradient experiment
sparse_chol <- function(p, d) {
	n <- 20*p
	res <- list()
	
	Ltrue <- covchol::rlower(p = p, d = d)
	diag(Ltrue) <- 1
	Sigmatrue <- Ltrue %*% t(Ltrue)

	X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
	train <- floor(n/2)
	Cortrue <- cov2cor(Sigmatrue)
	
	Cortrain <- cor(X[1:train, ])

	llpath <- covchol::llpathL(Sigma = Cortrain)
	frobs <- lapply(X = llpath, FUN = function(res) {
		ldebias <- covchol::proxgradL(Sigma = Cortrain, L = res$L, lambda = 0, eps = 1e-10)
		norm(ldebias$L %*% t(ldebias$L) - cor(X[(train + 1):n, ]), type = "F")
	})

	Lest <- llpath[[which.min(frobs)]]$L
	Corest <- Lest %*% t(Lest)
	
	return(list("ltrue" = Ltrue, "lest" = Lest))
}

band_chol <- function(p, d) {
	n <- 2*p
	res <- list()
	
	Ltrue <- covchol::rlower(p = p, d = d)
	diag(Ltrue) <- 1
	Sigmatrue <- Ltrue %*% t(Ltrue)
	
	X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigmatrue)
	
	Sigmaest <- PDSCE::band.chol.cv(x = X)$sigma
	Lest <- zapsmall(t(chol(Sigmaest)))
		
	return(list("ltrue" = Ltrue, "lest" = Lest))
}

execute_experiment(p = p, r = r, ename = "sparse_chol", emethod = sparse_chol)
execute_experiment(p = p, r = r, ename = "band_chol", emethod = band_chol)