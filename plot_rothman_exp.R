library("ggplot2")
library("dplyr")

r <- 200
p <- c(10, 30, 50, 100, 200, 500, 1000)

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

get_statistics <- function(p, r, ename) {
	fstat <- c("tpr" = stat_tpr,
						 "tnr" = stat_tnr,
						 "frob" = stat_frob,
						 "f1" = stat_f1)
	methods <- c("sparse", "band")
	data <- array(
		dim = c(length(p), 3, length(methods), length(fstat)),
		dimnames = list(p = p, sigma = 1:3, methods = methods, fstat = names(fstat))
	)
	data_sd <- array(
		dim = c(length(p), 3, length(methods), length(fstat)),
		dimnames = list(p = p, sigma = 1:3, methods = methods, fstat = names(fstat))
	)
	stat_res <- array(
		dim = c(r, length(fstat), length(methods)),
		dimnames = list(r = 1:r, fstat = names(fstat), methods = methods)
	)
		
	for (i in 1:length(p)) {
		for (sigma in 1:3) {
			for (k in 1:r) {
				res <- readRDS(file = paste0("rothman_exp/sigma", sigma, "_", p[i], "_r", k, ".rds"))
				for (l in seq(length(fstat))) {
					stat_res[k, l, "sparse"] <- fstat[[l]](res$sigmatrue, res$sigmasparse)
					stat_res[k, l, "band"] <- fstat[[l]](res$sigmatrue, res$sigmaband)
				}
			}
			for (m in methods) {
				for (l in seq(length(fstat))) {
					data[i, sigma, m, l] <- stats::mean(stat_res[, l, m])
					data_sd[i, sigma, m, l] <- stats::sd(stat_res[, l, m])
				}
			}
		}
	}
	
	df <- data %>% as.tbl_cube(met_name = "data") %>% as_tibble()
	df$m <- as.factor(df$m)
	df$sigma <- as.factor(df$sigma)
	df_sd <- data_sd %>% as.tbl_cube(met_name = "data_sd") %>% as_tibble()
	df$data_sd <- df_sd$data_sd
	
	return(df)
}

plot_comparison <- function(df, plot_title = "", plot_ylab = "") {
	
	pl <- ggplot(df, aes(x = p, y = data, group = m)) +
		facet_grid(cols = vars(sigma), rows = vars(fstat), 
							 labeller = labeller(fstat = toupper),
							 scales = "free") +
		geom_line(aes(color = m)) +
		geom_point(aes(color = m)) +
		theme_bw() +
		theme(text = element_text(size = 20), legend.position = "bottom") +
		xlab("Number of nodes (p)") +
		ylab("") 
	
	pl <- pl +
		geom_ribbon(aes(ymin = data - data_sd, ymax = data + data_sd, fill = m),
								alpha = .2) +
		labs(fill = "Method", color = "Method") #+
		#scale_fill_discrete(labels = c("Likelihood")) +
		#scale_color_discrete(labels = c("Likelihood"))
	
	return(pl)
}

df <- get_statistics(p = p, r = r, ename = "rothman_exp")
pl <- plot_comparison(df)
ggplot2::ggsave(filename = "rothman_exp.pdf", plot = pl, device = "pdf", width = 11, height = 9,
								path = "../sparsecholeskycovariance/img/")

