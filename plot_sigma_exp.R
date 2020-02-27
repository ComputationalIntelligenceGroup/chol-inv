library("ggplot2")
library("dplyr")

get_norms_sigma_exp <- function(p, r) {

	method <- c("sparse", "band")
	norms <- array(
		dim = c(length(p), 3, length(method)),
		dimnames = list(p = p, d = 1:3, method = method)
	)
	norms_sd <- array(
		dim = c(length(p), 3, length(method)),
		dimnames = list(p = p, d = 1:3, method = method)
	)
	norms_res <- matrix(nrow = r, ncol = length(method), 
											dimnames = list(r = 1:r, method = method))
	
	for (i in 1:length(p)) {
		d <- c(1/p[i], 2/p[i], 3/p[i])
		for (j in seq_along(d)) {
			for (k in 1:r) {
				res_sparse <- readRDS(file = paste0("sigma_exp/sigmasparse_", p[i], "_", d[j], "_r", k, ".rds"))
				res_band <- readRDS(file = paste0("sigma_exp/sigmaband_", p[i], "_", d[j], "_r", k, ".rds"))
				res_true <- readRDS(file = paste0("sigma_exp/sigmatrue_", p[i], "_", d[j], "_r", k, ".rds"))
				norms_res[k, "sparse"] <- norm(res_true - res_sparse, type = "2") 
				norms_res[k, "band"] <- norm(res_true - res_band, type = "2")
			}
			norms[i, j, "sparse"] <- mean(norms_res[, "sparse"])
			norms[i, j, "band"] <- mean(norms_res[, "band"])
			
			norms_sd[i, j, "sparse"] <- stats::sd(norms_res[, "sparse"])
			norms_sd[i, j, "band"] <- stats::sd(norms_res[, "band"])
		}
	}
	
	df <- norms %>% as.tbl_cube(met_name = "norms") %>% as_tibble()
	df$method <- as.factor(df$method)
	df_sd <- norms_sd %>% as.tbl_cube(met_name = "norms_sd") %>% as_tibble()
	df$norms_sd <- df_sd$norms_sd
	
	return(df)
}

plot_sigma_exp <- function(df, plot_title = "", plot_ylab = "") {
	
	lab_densities <- function(str) {
		return(paste0("Density = ", str, "/p"))
	}
	
	pl <- ggplot(df, aes(x = p, y = norms, group = method)) +
		facet_grid(cols = vars(d), labeller = labeller(d = lab_densities),
							 scales = "free") +
		geom_line(aes(color = method)) +
		geom_point(aes(color = method)) +
		theme_bw() +
		theme(text = element_text(size = 20), legend.position = "bottom") +
		xlab("Number of nodes (p)") +
		ylab("") 
	
		pl <- pl +
			geom_ribbon(aes(ymin = norms - norms_sd, ymax = norms + norms_sd, fill = method),
									alpha = .2) +
			labs(fill = "Method", color = "Method") +
			scale_fill_discrete(labels = c("Banding", "Likelihood")) +
			scale_color_discrete(labels = c("Banding", "Likelihood"))
	
	return(pl)
}

r <- 200
p <- c(30, 100, 200, 500, 1000)
df <- get_norms_sigma_exp(p = p, r = r)
pl <- plot_sigma_exp(df)
ggplot2::ggsave(filename = "sigma_exp.pdf", plot = pl, device = "pdf", width = 12, height = 4,
								path = "../sparsecholeskycovariance/img/")
	