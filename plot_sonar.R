library("ggplot2")
library("dplyr")

method <- c("grad_lik",
			"grad_frob",
			"band",
			"lasso")
type <- c("M" = 1, "R" = 2)

get_choleskys <- function(dirname) {

	chols <- array(dim = c(length(type), length(method), 60, 60),
				dimnames = list(type = names(type), method = method,
								rows = 1:60, cols = 1:60))
	for (m in method) {
		for (t in names(type)) {
			chols[t, m, , ] <- 
				readRDS(file = paste0(dirname, m, "_", t, ".rds"))
		}
	}

	return(chols)
}

plot_heatmaps <- function(chols, fname) {
	
	df <- chols %>% as.tbl_cube(met_name = "covs") %>% as_tibble()
	df$type <- as.factor(df$type)
	df$method <- as.factor(df$method)

	pl <- ggplot(df, aes(x = cols, y = rows, z = covs, fill = covs)) +
		facet_grid(rows = vars(type), cols = vars(method)) +
		geom_tile() + coord_equal() +
		geom_raster(alpha = 0.75) +
		scale_fill_distiller(palette = "Spectral") +
		theme_bw() +
		theme(legend.position = "bottom", text = element_text(size = 20)) +
		xlab("") +
		ylab("") +
		ylim(60, 0)

	ggsave(filename = paste0(fname, ".pdf"), plot = pl, device = "pdf",
		path = "../sparsecholeskycovariance/img/", width = 11, height = 6)
}
plot_scree <- function(chols) {
	eigens <- array(dim = c(length(type), length(method), 60),
			dimnames = list(type = type, method = method, value = 1:60))
			
	for (t in type) {
		for (m in method) {
			sigma <- chols[t, m, , ] %*% t(chols[t, m, , ])
			eigens[t, m, ] <- eigen(sigma, symmetric = TRUE, 
									only.values = TRUE)$values
		}
	}

	df <- eigens %>% as.tbl_cube(met_name = "eigens") %>% as_tibble()
	df$type <- as.factor(df$type)
	df$method <- as.factor(df$method)

	pl <- ggplot(df, aes(x = value, y = eigens, group = method, color = method)) +
		facet_grid(cols = vars(type)) +
		geom_line() +
		theme_bw() +
		theme(legend.position = "bottom", text = element_text(size = 20)) +
		xlab("") +
		ylab("Eigenvalue") 

	ggsave(filename = paste0("sonar_eigens.pdf"), plot = pl, device = "pdf",
		path = "../sparsecholeskycovariance/img/", width = 11, height = 5)
}

### Plot covariance heatmaps and eigenvalues
dirname <- "sonar_heat_exp/"
chols <- get_choleskys(dirname)
plot_heatmaps(chols, fname = "sonar_chols")
covs <- array(dim = c(length(type), length(method), 60, 60),
				dimnames = list(type = names(type), method = method,
								rows = 1:60, cols = 1:60))
for (t in type) {
	for (m in method) {
		covs[t, m, , ] <- chols[t, m, , ] %*% t(chols[t, m, , ])
	}
}
plot_heatmaps(covs, fname = "sonar_covs")
plot_scree(chols)

### Generate predictions
preds <- readRDS(file = "data/preds.rds")
stat_tpr <- function(pred, true) {
	positive <- 1 # Mines
	if (sum(true == positive) == 0) {
		return(1)
	}
	tp <- sum(true == positive & pred == positive)
	return(tp/(sum(true == positive)))
}
stat_tnr <- function(pred, true) {
	negative <- 2 # Rocks
	if (sum(true == negative) == 0) {
		return(1)
	}
	tn <- sum(true == negative & pred == negative)
	return(tn/sum(true == negative))
}
stat_f1 <- function(pred, true, val) {
	positive <- 1 # Mines
	negative <- 2 # Rocks
	tp <- sum(true == positive & pred == positive)
	fn <- sum(pred == negative & true == positive)
	fp <- sum(true == negative & pred == positive)
	return(2*tp/(2*tp + fp + fn))
}

est <- c("grad_frob", "band", "grad_lik", "lasso")
fstat <- c("tpr" = stat_tpr,
			"tnr" = stat_tnr,
			 "f1" = stat_f1)
stat <- matrix(
	nrow = length(fstat), ncol = length(est),
	dimnames = list(fstat = names(fstat), est = est)
)

for (e in est) {
	for (f in names(fstat)) {
		stat[f, e] <- fstat[[f]](pred = preds[, e], 
								true = preds[, "true"])
	}
}

saveRDS(stat, file = "data/stat_preds.rds")

