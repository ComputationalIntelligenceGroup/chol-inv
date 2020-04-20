library("ggplot2")
library("dplyr")

covs <- readRDS(file = "data/covs.rds")

### Plot covariance heatmaps
df <- covs %>% as.tbl_cube(met_name = "covs") %>% as_tibble()
df$type <- as.factor(df$type)
df$est <- as.factor(df$est)

pl <- ggplot(df, aes(x = rows, y = cols, z = covs, fill = covs)) +
	facet_grid(rows = vars(type), cols = vars(est)) +
	geom_tile() + coord_equal() +
	geom_raster(alpha = 0.75) +
	scale_fill_distiller(palette = "Spectral") +
	theme_bw() +
	theme(legend.position = "bottom", text = element_text(size = 20)) +
	xlab("") +
	ylab("") +
	ylim(60, 0)

ggsave(filename = paste0("sonar_covs.pdf"), plot = pl, device = "pdf",
	path = "../sparsecholeskycovariance/img/", width = 11, height = 6)

### Plot eigenvalue freqpol (scree plot)
type <- c("R", "M")
est <- c("sample", "grad_frob", "grad_lik", "band", "lasso")
eigens <- array(dim = c(length(type), length(est), 60),
			dimnames = list(type = type, est = est, value = 1:60))
			
for (t in type) {
	for (e in est) {
		eigens[t, e, ] <- eigen(covs[t, e, , ], symmetric = TRUE, 
			only.values = TRUE)$values
	}
}

df <- eigens %>% as.tbl_cube(met_name = "eigens") %>% as_tibble()
df$type <- as.factor(df$type)
df$est <- as.factor(df$est)

pl <- ggplot(df, aes(x = value, y = eigens, group = est, color = est)) +
	facet_grid(cols = vars(type)) +
	geom_line() +
	theme_bw() +
	theme(legend.position = "bottom", text = element_text(size = 20)) +
	xlab("") +
	ylab("Eigenvalue") 

ggsave(filename = paste0("sonar_eigens.pdf"), plot = pl, device = "pdf",
	path = "../sparsecholeskycovariance/img/", width = 11, height = 5)


preds <- readRDS(file = "data/preds.rds")
stat_tpr <- function(pred, true, val) {
	if (sum(true == val) == 0) {
		return(1)
	}
	tp <- sum(true == val & pred == val)
	return(tp/(sum(true == val)))
}
stat_tnr <- function(pred, true, val) {
	if (sum(true != val) == 0) {
		return(1)
	}
	tn <- sum(true != val & pred != val)
	return(tn/sum(true != val))
}
stat_f1 <- function(pred, true, val) {
	tp <- sum(true == val & pred == val)
	fn <- sum(pred != val & true == val)
	fp <- sum(true != val & pred == val)
	return(2*tp/(2*tp + fp + fn))
}

est <- c("grad_frob", "band")
fstat <- c("tpr" = stat_tpr,
			"tnr" = stat_tnr,
			 "f1" = stat_f1)
stat <- array(
	dim = c(length(type), length(fstat), length(est)),
	dimnames = list(type = type, fstat = names(fstat), est = est)
)

for (t in type) {
	for (e in est) {
		for (f in names(fstat)) {
			stat[t, f, e] <- fstat[[f]](pred = preds[, e], 
										true = preds[, "true"], 
										val = t)
		}
	}
}

saveRDS(stat, file = "data/stat_preds.rds")

