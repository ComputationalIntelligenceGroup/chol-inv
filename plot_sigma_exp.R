library("ggplot2")
library("dplyr")


stat_tpr <- function(ltrue, lest) {
	p <- ncol(ltrue)
	return((sum(ltrue != 0 & lest != 0) - p)/(sum(ltrue != 0) - p))
}
stat_tnr <- function(ltrue, lest) {
	p <- ncol(ltrue)
	tn <- sum(ltrue == 0 & lest == 0) - p*(p - 1)/2
	return(tn/(sum(ltrue == 0) - (p * (p - 1) / 2)))
}
stat_frob <- function(ltrue, lest) {
	return(norm(lest - ltrue, type = "F"))
}
stat_f1 <- function(ltrue, lest) {
	p <- ncol(ltrue)
	tp <- sum(ltrue != 0 & lest != 0) - p
	fn <- sum(ltrue != 0 & lest == 0)
	fp <- sum(ltrue == 0 & lest != 0)
	return(2*tp/(2*tp + fp + fn))
}

get_frobs_sigma_exp <- function(p, r) {

	methods <- c("sparse", "band")
	frobs <- array(
		dim = c(length(p), 3, length(methods)),
		dimnames = list(p = p, d = 1:3)
	)
	frobs_sd <- array(
		dim = c(length(p), 3, length(methods)),
		dimnames = list(p = p, d = 1:3)
	)
	frobs_res <- matrix(nrow = r, ncol = 2, 
											dimnames = list(r = 1:r, method = c("sparse", "band")))
	
	for (i in 1:length(p)) {
		d <- c(1/p[i], 2/p[i], 3/p[i])
		for (j in seq_along(d)) {
			for (k in 1:r) {
				res <- readRDS(file = paste0("sigma_exp/", p[i], "_", d[j], "_r", k, ".rds"))
				frobs_res[k, "sparse"] <- stat_frob(res$sigmatrue, res$sigmasparse) 
				frobs_res[k, "band"] <- stat_frob(res$sigmatrue, res$sigmaband) 

				frobs[i, j, "sparse"] <- stats::mean(frobs_res[, "sparse"])
				frobs[i, j, "band"] <- stats::mean(frobs_res[, "band"])

				frobs_sd[i, j, "sparse"] <- stats::sd(frobs_res[, "sparse"])
				frobs_sd[i, j, "band"] <- stats::sd(frobs_res[, "band"])
			}
		}
	}
	
	df <- frobs %>% as.tbl_cube(met_name = "frobs") %>% as_tibble()
	df$method <- as.factor(df$method)
	df_sd <- frobs_sd %>% as.tbl_cube(met_name = "frobs_sd") %>% as_tibble()
	df$frobs_sd <- df_sd$frobs_sd
	
	return(df)
}

plot_sigma_exp <- function(df, plot_title = "", plot_ylab = "") {
	
	lab_densities <- function(str) {
		return(paste0(str, "/p"))
	}
	
	pl <- ggplot(df, aes(x = p, y = frobs, group = method)) +
		facet_grid(cols = vars(d), labeller = labeller(d = lab_densities),
							 scales = "free") +
		geom_line(aes(color = method)) +
		geom_point(aes(color = method)) +
		theme_bw() +
		theme(text = element_text(size = 20), legend.position = "bottom") +
		xlab("Number of nodes (p)") +
		ylab("") 
	
		pl <- pl +
			geom_ribbon(aes(ymin = frobs - frobs_sd, ymax = frobs + frobs_sd, fill = method),
									alpha = .2) +
			labs(fill = "Method", color = "Method") +
			scale_fill_discrete(labels = c("Banding", "Likelihood")) +
			scale_color_discrete(labels = c("Banding", "Likelihood"))
	
	return(pl)
}

r <- 200
p <- c(10, 30, 50, 100, 200, 500, 1000)
df <- get_frobs_sigma_exp(p = p, r = r)
pl <- plot_sigma_exp(df)
ggplot2::ggsave(filename = "sigma_exp.pdf", plot = pl, device = "pdf", width = 11, height = 9,
								path = "../sparsecholeskycovariance/img/")

	