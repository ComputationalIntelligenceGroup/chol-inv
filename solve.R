library("ggplot2")

devtools::install_github("irenecrsn/ggmexp")
devtools::install_github("irenecrsn/covchol")

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
		O_inv_chol <- covchol::chol_inv(B)	

		norm(O_inv_chol - O_inv_solve)
	}, 
		error = function(e) {
		print(paste0("Error at p = ", p, " and d = ", d, ":"))
		print(e)
		return(NA)
	})
	
	return(res)
}

wd <- getwd()
plot_dname <- "plot_true"
dir.create(paste0(wd, "/", plot_dname), showWarnings = FALSE)
ggmexp::execute(p = p, d = d, r = r, experiment = true)

pl <- ggmexp::plot_map_reduce(p = p, d = d, r = 1, N = 1, reduce = mean, exp_name = "true") 
pl <- pl + ylab("") + ggtitle("")
ggsave(filename = "true.pdf", plot = pl, device = "pdf", path = plot_dname)
