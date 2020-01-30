library("ggplot2")

devtools::install_github("irenecrsn/ggmexp")
devtools::install_github("irenecrsn/covchol")

r <- 10
p <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 250, 300, 400, 500, 750, 1000)
d <- c(0.0025, 0.005, 0.025, 0.05, 0.25, 0.5)

##### First experiment --> comparison with ggm::fitDag
learn_chol_conc <- function(p, d, N) {
	dag <- pcalg::randomDAG(n = p, prob = d)
	data <- pcalg::rmvDAG(n = N, dag = dag)
	amat <- t(as(dag, "matrix"))
	chol_conc_fit <- covchol::fit_chol_conc(amat = amat, data = data)
	chol_conc_true <- diag(p) - pcalg::wgtMatrix(dag)
	
	return(norm(chol_conc_fit - chol_conc_true))
}

ggmexp::execute(p = p, d = d, r = r, experiment = learn_chol_conc, N = 2*p)

wd <- getwd()
plot_dname <- "plot_comp_true"
dir.create(paste0(wd, "/", plot_dname), showWarnings = FALSE)
pl <- ggmexp::plot_map_reduce(p = p, d = d, r = 1, N = 1, reduce = mean, exp_name = "comp_true") 
pl <- pl + ylab("") + ggtitle("")
ggsave(filename = "comp_true.pdf", plot = pl, device = "pdf", path = plot_dname)
