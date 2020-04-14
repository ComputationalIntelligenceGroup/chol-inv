library("ggplot2")
library("dplyr")

type <- c("R", "M")
est <- c("sample", "sparse", "sparse_f", "band", "lasso")
covs <- readRDS(file = "data/covs.rds")

### Plot covariance heatmaps
df <- covs %>% as.tbl_cube(met_name = "covs") %>% as_tibble()
df$type <- as.factor(df$type)
df$est <- as.factor(df$est)

pl <- ggplot(df, aes(x = rows, y = cols, z = covs, fill = covs)) +
	facet_grid(cols = vars(type), rows = vars(est)) +
	geom_tile() + coord_equal() +
	geom_contour(color = "white", alpha = 0.75) +
	scale_fill_distiller(palette = "Spectral", na.value = "white") +
	theme_bw() +
	xlab("") +
	ylab("") +
	ylim(60, 0)

ggsave(filename = paste0("sonar_covs.pdf"), plot = pl, device = "pdf",
	path = "../sparsecholeskycovariance/img/")

### Plot eigenvalue freqpol (scree plot)
eigens <- array(dim = c(length(type), length(est), 60),
			dimnames = list(type = type, est = est, value = 1:60))
			
for (t in type) {
	for (e in est) {
		eigens[t, e, ] <- eigen(covs[t, e, , ])$values
	}
}

df <- eigens %>% as.tbl_cube(met_name = "eigens") %>% as_tibble()
df$type <- as.factor(df$type)
df$est <- as.factor(df$est)

pl <- ggplot(df, aes(eigens, group = est, color = est)) +
	facet_grid(cols = vars(type)) +
	geom_freqpoly(bins = nclass.FD(df$eigens)) +
	theme_bw() +
	xlab("Eigenvalue") +
	ylab("Absolute frequency") 

ggsave(filename = paste0("sonar_eigens.pdf"), plot = pl, device = "pdf",
	path = "../sparsecholeskycovariance/img/")

