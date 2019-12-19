source("utils.R")
library("ggplot2")

devtools::install_github("irenecrsn/ggmexp")

r <- 10
p <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 250, 300, 400, 500, 750, 1000)
d <- c(0.0025, 0.005, 0.025, 0.05, 0.25, 0.5)

##### First experiment -> True regression coefficients available
true <- function(p, d) {
	dag <- pcalg::randomDAG(n = p, prob = d)
	B <- pcalg::wgtMatrix(dag)
	L <- diag(p) - B
	res <- tryCatch(
	{
		O_inv_solve <- solve(L)
		O_inv_chol <- chol_inv(B)	

		return(norm(O_inv_chol - O_inv_solve))
	}, 
		error = function(e) {
		print(paste0("Error at p = ", p, " and d = ", d, ":"))
		print(e)
		return(NA)
	})
	
	return(res)
}

##### Second experiment -> Estimated regression coefficients with standard MLE
mle <- function(p, d, N) {
	dag <- pcalg::randomDAG(n = p, prob = d)
	data <- pcalg::rmvDAG(n = N, dag = dag)

	# Maximum likelihood estimation via regression
	dag_am <- ggm::grMAT(dag)
	fit <- ggm::fitDag(dag_am, cov(data), N)
	L <- fit$Ahat
	D <- solve(diag(fit$Dhat))
	
	# Calculating the inverse with the two methods
	O_solve <- solve(L)
	O_inv <- chol_inv(diag(p) - L)
	
	S_true <- pcalg::trueCov(dag)
	S_inv <- O_inv %*% D %*% t(O_inv)
	S_solve <- O_solve %*% D %*% t(O_solve)
	
	norms <- eval(formals(base::norm)$type)
	
	norm_inv_solve <- sapply(norms, norm, x = S_inv - S_solve)
	norm_inv_true <- sapply(norms, norm, x = S_inv - S_true)
	norm_solve_true <- sapply(norms, norm, x = S_solve - S_true)
	
	S_norms <- data.frame(norm_inv_solve, norm_inv_true, norm_solve_true)
	
	return(S_norms)
}

##### Third experiment -> Estimating both structure and coefficients 
pcmle <- function(p, d, N) {
	dag <- pcalg::randomDAG(n = p, prob = d)
	data <- pcalg::rmvDAG(n = N, dag = dag)
	pdag_learned <- pcalg::pc(suffStat = list(C = cov(data), n = p),
							indepTest = pcalg::gaussCItest,
							alpha = 0.01,
							p = p)
	dag_learned <- pcalg::pdag2dag(pdag_learned@graph)$graph

	# Maximum likelihood estimation via regression
	dag_am <- ggm::grMAT(dag_learned)
	colnames(data) <- colnames(data, do.NULL = FALSE, prefix = "")
	fit <- ggm::fitDag(dag_am, cov(data), N)
	L <- fit$Ahat
	D <- solve(diag(fit$Dhat))

	# Calculating the inverse with the two methods
	O_solve <- solve(L)
	O_inv <- chol_inv(diag(p) - L)

	S_true <- pcalg::trueCov(dag)
	S_inv <- O_inv %*% D %*% t(O_inv)
	S_solve <- O_solve %*% D %*% t(O_solve)

	norms <- eval(formals(base::norm)$type)

	norm_inv_solve <- sapply(norms, norm, x = S_inv - S_solve)
	norm_inv_true <- sapply(norms, norm, x = S_inv - S_true)
	norm_solve_true <- sapply(norms, norm, x = S_solve - S_true)

	S_norms <- data.frame(norm_inv_solve, norm_inv_true, norm_solve_true)

	print(S_norms)
	print(norm(O_inv - O_solve))
}

wd <- getwd()
plot_dname <- "plot_true"
dir.create(paste0(wd, "/", plot_dname), showWarnings = FALSE)
ggmexp::execute(p = p, d = d, r = r, experiment = true)
pl <- ggmexp::plot_map_reduce(p = p, d = d, r = r, N = 1, reduce = mean, exp_name = "true") 
pl <- pl + ylab("") + ggtitle("")
ggsave(filename = "true.pdf", plot = pl, device = "pdf", path = plot_dname)

#ggmexp::execute(p = p, d = d, r = r, experiment = mle, N = 2*p)
#ggmexp::execute(p = p, d = d, r = r, experiment = pcmle, N = 2*p)
