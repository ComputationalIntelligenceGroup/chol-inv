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
}

##### Second experiment -> Estimated regression coefficients

#N <- 50
#dag <- pcalg::randomDAG(n = p, prob = d)
#data <- pcalg::rmvDAG(n = N, dag = dag)

# Maximum likelihood estimation via regression
#dag_am <- ggm::grMAT(dag)

#colnames(data) <- colnames(data, do.NULL = FALSE, prefix = "")
#fit_ggm <- ggm::fitDag(dag_am, cov(data), N)
#norms_ggm <- get_norms(U = fit_ggm$Ahat, D = diag(fit_ggm$Dhat), dag = dag)
