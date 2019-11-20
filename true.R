source("utils.R")

p <- 1000
d <- 0.0025

##### First experiment -> True regression coefficients available
rep <- 10
error_solve <- 0
error_chol <- 0
for (i in 1:rep) {
	dag <- pcalg::randomDAG(n = p, prob = d)
	B <- pcalg::wgtMatrix(dag)
	L <- diag(p) - B
	O_solve <- tryCatch(solve(L), error = function(e) {
			print(e)
			error_solve <- error_solve + 1
		})
	O_inv <- tryCatch(chol_inv(B), error = function(e)
		{
			print(e)
			error_chol <- error_chol + 1
		})
	
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
