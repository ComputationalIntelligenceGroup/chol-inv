library("ggplot2")
library("dplyr")

source("real_lib.R")

method <- c("band",
			"grad_frob",
			"grad_lik",
			"lasso")

get_choleskys <- function(dirname) {

	chols <- array(dim = c(length(type), length(method), class - 1, class - 1),
				dimnames = list(type = type, method = method,
								rows = 1:(class - 1), cols = 1:(class - 1)))
	for (m in method) {
		for (t in type) {
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
		ylim(class, 0)

	ggsave(filename = paste0(fname, ".pdf"), plot = pl, device = "pdf",
		path = "../sparsecholeskycovariance/img/", width = 11, height =
		3*length(type))
}
plot_scree <- function(chols) {
	eigens <- array(dim = c(length(type), length(method), class - 1),
			dimnames = list(type = type, method = method, value = 1:(class - 1)))
			
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

	ggsave(filename = paste0(exp, "_eigen.pdf"), plot = pl, device = "pdf",
		path = "../sparsecholeskycovariance/img/", width = 11, height = 5)
}

### Plot covariance heatmaps and eigenvalues
dirname <- paste0(exp, "_heat_exp/")
chols <- get_choleskys(dirname)
plot_heatmaps(chols, fname = paste0(exp, "_chols"))
covs <- array(dim = c(length(type), length(method), class - 1, class - 1),
				dimnames = list(type = type, method = method,
								rows = 1:(class - 1), cols = 1:(class - 1)))
for (t in type) {
	for (m in method) {
		covs[t, m, , ] <- chols[t, m, , ] %*% t(chols[t, m, , ])
	}
}
plot_heatmaps(covs, fname = paste0(exp, "_covs"))
plot_scree(chols)

### Generate predictions
stat_tpr <- function(pred, true, t) {
	if (sum(true == t) == 0) {
		return(1)
	}
	tp <- sum(true == t & pred == t)
	return(tp/(sum(true == t)))
}
stat_tnr <- function(pred, true, t) {
	if (sum(true != t) == 0) {
		return(1)
	}
	tn <- sum(true != t & pred != t)
	return(tn/sum(true != t))
}
stat_f1 <- function(pred, true, t) {
	tp <- sum(true == t & pred == t)
	fn <- sum(pred != t & true == t)
	fp <- sum(true != t & pred == t)
	return(2*tp/(2*tp + fp + fn))
}

fstat <- c("tnr" = stat_tnr,
		   "f1" = stat_f1)
stat <- array(
	dim = c(length(fstat), length(method), length(type)),
	dimnames = list(fstat = names(fstat), method = method, type = type)
)
preds <- readRDS(file = paste0("data/", exp, "_preds.rds"))

for (i in 1:length(type)) {
	for (m in method) {
		for (f in names(fstat)) {
			stat[f, m, type[i]] <- fstat[[f]](pred = preds[, m], 
										true = preds[, "true"],
										t = i)
		}
	}
}
print(stat)
saveRDS(stat, file = paste0("data/", exp, "_stat_preds.rds"))

