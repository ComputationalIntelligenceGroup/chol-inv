devtools::install_github("irenecrsn/ggmexp")
devtools::install_github("irenecrsn/covchol")

r <- 10
p <- c(10, 30, 50, 100)
d <- c(0.0025, 0.005, 0.025, 0.05, 0.25, 0.5)

##### Gradient experiment
gradient <- function(p, d, N) {
	
	Ltrue <- covchol::rlower(p = p, d = d)
	diag(Ltrue) <- 1
	D <-  diag(p)
	Sigmatrue <- Ltrue %*% D %*% t(Ltrue)
	
	X <- MASS::mvrnorm(N, rep(0, p), Sigma = Sigmatrue)
	Cortrue <- cov2cor(Sigmatrue)

	Lstart <- t(chol(cor(X)))
	diag(Lstart) <- 1
	
	Lest <- covchol::proxgradL(Sigma = cor(X), L = Lstart, D = diag(p),
										lambda = 0.02, h = TRUE, trace = 0, alpha = 1, beta = 0.5, 
										maxIter = 100, eps = 0.0001)
	Corest <- Lest %*% t(Lest)
	res_mll <- covchol::mll(solve(Corest), Cortrue)
	min_mll <- covchol::mll(solve(Cortrue), Cortrue)
	
	
	return(list("res_mll" = res_mll, "min_mll" = min_mll,
							"ltrue" = Ltrue, "lest" = Lest))
}

ggmexp::execute_experiment(p = p, d = d, r = r, N = 10*p, ename = "gradient", emethod = gradient)