source("utils.R")

p <- 1000
d <- 0.0025

##### First experiment -> True regression coefficients available
rep <- 10
error_solve <- 0
error_chol <- 0
for (i in 1:rep) {
	dag <- pcalg::randomDAG(n = p, prob = d)
	b_true <- pcalg::wgtMatrix(dag)
	o_solve <- tryCatch(solve(diag(p) - b_true), error = function(e) {
			print(e)
			error_solve <- error_solve + 1
		})
	o_chol <- tryCatch(chol_inv(b_true), error = function(e)
		{
			print(e)
			error_chol <- error_chol + 1
		})
	print(norm(o_chol - o_solve))
}

#tryCatch(S_true <- pcalg::trueCov(dag), error = function(e) { print(e) })
