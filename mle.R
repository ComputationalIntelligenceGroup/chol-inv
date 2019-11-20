source("utils.R")

##### Second experiment -> Estimated regression coefficients with standard MLE
p <- 1000
d <- 0.025
N <- 2000
rep <- 10
error_solve <- 0
error_chol <- 0
for (i in 1:rep) {
	dag <- pcalg::randomDAG(n = p, prob = d)
	data <- pcalg::rmvDAG(n = N, dag = dag)

	# Maximum likelihood estimation via regression
	dag_am <- ggm::grMAT(dag)
	colnames(data) <- colnames(data, do.NULL = FALSE, prefix = "")
	fit <- ggm::fitDag(dag_am, cov(data), N)
	L <- fit$Ahat
	D <- solve(diag(fit$Dhat))

	# Calculating the inverse with the two methods
	O_solve <- tryCatch(solve(L), error = function(e) {
			print(e)
			error_solve <- error_solve + 1
		})
	O_chol <- tryCatch(chol_inv(diag(p) - L), error = function(e)
		{
			print(e)
			error_chol <- error_chol + 1
		})
	S_true <- pcalg::trueCov(dag)
	S_chol <- O_chol %*% D %*% t(O_chol)
	S_solve <- O_solve %*% D %*% t(O_solve)
	
	norms <- eval(formals(base::norm)$type)
	
	norm_chol_solve <- sapply(norms, norm, x = S_chol - S_solve)
	norm_chol_true <- sapply(norms, norm, x = S_chol - S_true)
	norm_solve_true <- sapply(norms, norm, x = S_solve - S_true)
	
	S_norms <- data.frame(norm_chol_solve, norm_chol_true, norm_solve_true)
	
	print(S_norms)
	print(norm(O_chol - O_solve))
}
