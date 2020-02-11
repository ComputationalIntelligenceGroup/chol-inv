library("ggplot2")
library("dplyr")

r <- 200
p <- c(10, 30, 50, 100)

stat_tpr <- function(ltrue, lest) {
	p <- ncol(ltrue)
	return((sum(ltrue != 0 & lest != 0) - p)/(sum(ltrue != 0) - p))
}
stat_tnr <- function(ltrue, lest) {
	p <- ncol(ltrue)
	tn <- sum(ltrue == 0 & lest == 0) - p*(p - 1)/2
	return(tn/(sum(ltrue == 0) - (p * (p - 1) / 2)))
}
stat_acc <- function(ltrue, lest) {
	p <- ncol(ltrue)
	tp <- sum(ltrue != 0 & lest != 0) - p
	tn <- sum(ltrue == 0 & lest == 0) - p*(p - 1)/2
	return((tp + tn)/(p*(p - 1)/2))
}
stat_f1 <- function(ltrue, lest) {
	p <- ncol(ltrue)
	tp <- sum(ltrue != 0 & lest != 0) - p
	fn <- sum(ltrue != 0 & lest == 0)
	fp <- sum(ltrue == 0 & lest != 0)
	return(2*tp/(2*tp + fp + fn))
}

get_statistics <- function(p, r, ename) {
	fstat <- c("tpr" = stat_tpr,
						 "tnr" = stat_tnr,
						 "acc" = stat_acc,
						 "f1" = stat_f1)
	data <- array(
		dim = c(length(p), 3, length(ename), length(fstat)),
		dimnames = list(p = p, d = 1:3, ename = ename, fstat = names(fstat))
	)
	data_sd <- array(
		dim = c(length(p), 3, length(ename), length(fstat)),
		dimnames = list(p = p, d = 1:3, ename = ename, fstat = names(fstat))
	)
	stat_res <- matrix(nrow = r, ncol = length(fstat))
	
	for (i in 1:length(p)) {
		d <- c(1/p[i], 2/p[i], 3/p[i])
		for (j in seq_along(d)) {
			for (m in ename) {
				for (k in 1:r) {
					atomic_res <- readRDS(file = paste0(m, "/", p[i], "_", d[j], "_r", k, ".rds"))
					for (l in seq(length(fstat))) {
						stat_res[k, l] <- fstat[[l]](atomic_res$ltrue, atomic_res$lest) 
					}
				}
				for (l in seq(length(fstat))) {
					data[i, j, m, l] <- mean(stat_res[, l])
					data_sd[i, j, m, l] <- stats::sd(stat_res[, l])
				}
			}
		}
	}
	
	df <- data %>% as.tbl_cube(met_name = "data") %>% as_tibble()
	df$ename <- as.factor(df$ename)
	df_sd <- data_sd %>% as.tbl_cube(met_name = "data_sd") %>% as_tibble()
	df$data_sd <- df_sd$data_sd
	
	return(df)
}

plot_comparison <- function(df, plot_title = "", plot_ylab = "", ename) {
	
	lab_densities <- function(str) {
		return(paste0(str, "/p"))
	}
	
	pl <- ggplot(df, aes(x = p, y = data, group = ename)) +
		facet_grid(cols = vars(d), rows = vars(fstat), 
							 labeller = labeller(fstat = toupper, d = lab_densities)) +
		geom_line(aes(color = ename)) +
		geom_point(aes(color = ename)) +
		theme(text = element_text(size = 20), legend.position = "bottom") +
		xlab("Number of nodes (p)") +
		ylab("")

		pl <- pl +
			geom_ribbon(aes(ymin = data - data_sd, ymax = data + data_sd, fill = ename),
									alpha = .2) +
			labs(fill = "Method", color = "Method") +
			scale_fill_discrete(labels = c("Banding", "Likelihood")) +
			scale_color_discrete(labels = c("Banding", "Likelihood"))
	
	return(pl)
}

df <- get_statistics(p = p, r = r, ename = c("sparse_chol", "band_chol"))
pl <- plot_comparison(df, ename = c("sparse_chol", "band_chol"))
ggplot2::ggsave(filename = "stats.pdf", plot = pl, device = "pdf", width = 11, height = 7,
								path = "../sparsecholeskycovariance/img/")


################################################

	