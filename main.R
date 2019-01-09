source("utils.R")

p <- 10
d <- 0.25
N <- 50

dag <- pcalg::randomDAG(n = p, prob = d)
data <- pcalg::rmvDAG(n = N, dag = dag)

# Maximum likelihood estimation via regression
dag_am <- ggm::grMAT(dag)

colnames(data) <- colnames(data, do.NULL = FALSE, prefix = "")
fit_ggm <- ggm::fitDag(dag_am, cov(data), N)

norms_ggm <- get_norms(U = fit_ggm$Ahat, D = diag(fit_ggm$Dhat), dag = dag)
