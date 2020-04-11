library("ggplot2")
library("dplyr")

stat_tpr <- function(sigmatrue, sigmaest) {
	p <- ncol(sigmatrue)
	tp <- sum(sigmatrue != 0 & sigmaest != 0) - p
	return(tp/(sum(sigmatrue != 0) - p))
}
stat_tnr <- function(sigmatrue, sigmaest) {
	p <- ncol(sigmatrue)
	tn <- sum(sigmatrue == 0 & sigmaest == 0)
	if (sum(sigmatrue == 0) == 0) {
		return(1)
	}
	return(tn/sum(sigmatrue == 0))
}
stat_opnorm <- function(sigmatrue, sigmaest) {
	return(norm(sigmaest - sigmatrue, type = "2"))
}
stat_f1 <- function(sigmatrue, sigmaest) {
	p <- ncol(sigmatrue)
	tp <- sum(sigmatrue != 0 & sigmaest != 0) - p
	fn <- sum(sigmatrue != 0 & sigmaest == 0)
	fp <- sum(sigmatrue == 0 & sigmaest != 0)
	return(2*tp/(2*tp + fp + fn))
}

get_statistics <- function(p, r) {

	fstat <- c("tpr" = stat_tpr,
						 "tnr" = stat_tnr,
						 "opnorm" = stat_opnorm,
						 "f1" = stat_f1)
	method <- c("sample", "band", "sparse", "sparse_f")
	data <- array(
		dim = c(length(p), 3, length(method), length(fstat)),
		dimnames = list(p = p, d = 1:3, method = method, fstat = names(fstat))
	)
	data_se <- array(
		dim = c(length(p), 3, length(method), length(fstat)),
		dimnames = list(p = p, d = 1:3, method = method, fstat = names(fstat))
	)
	stat_res <- array(
		dim = c(r, length(fstat), length(method)),
		dimnames = list(r = 1:r, fstat = names(fstat), method = method)
	)
	
	for (i in 1:length(p)) {
		d <- c(1/p[i], 2/p[i], 3/p[i])
		for (j in seq_along(d)) {
			for (k in 1:r) {
				sigmasparse <- readRDS(file = paste0("sigma_exp/sigmasparse_", p[i], "_", d[j], "_r", k, ".rds"))
				sigmasparse_f <- readRDS(file = paste0("sigma_exp/sigmasparse_f_", p[i], "_", d[j], "_r", k, ".rds"))
				sigmaband <- readRDS(file = paste0("sigma_exp/sigmaband_", p[i], "_", d[j], "_r", k, ".rds"))
				sigmatrue <- readRDS(file = paste0("sigma_exp/sigmatrue_", p[i], "_", d[j], "_r", k, ".rds"))
				sigmasample <- readRDS(file = paste0("sigma_exp/sigmasample_", p[i], "_", d[j], "_r", k, ".rds"))
				for (l in seq(length(fstat))) {
					stat_res[k, l, "sparse_f"] <- fstat[[l]](sigmatrue, sigmasparse_f) 
					stat_res[k, l, "sparse"] <- fstat[[l]](sigmatrue, sigmasparse) 
					stat_res[k, l, "band"] <- fstat[[l]](sigmatrue, sigmaband)
					stat_res[k, l, "sample"] <- fstat[[l]](sigmatrue, sigmasample)
				}
			}
			for (m in method) {
				for (l in seq(length(fstat))) {
					data[i, j, m, l] <- mean(stat_res[, l, m])
					data_se[i, j, m, l] <- stats::sd(stat_res[, l, m])/sqrt(r)
				}
			}
		}
	}
	
	df <- data %>% as.tbl_cube(met_name = "data") %>% as_tibble()
	df$method <- as.factor(df$method)
	df_se <- data_se %>% as.tbl_cube(met_name = "data_se") %>% as_tibble()
	df$data_se <- df_se$data_se
	
	return(df)
}

plot_sigma_exp <- function(df, plot_title = "", plot_ylab = "") {
	
	lab_densities <- function(str) {
		return(paste0("Density = ", str, "/p"))
	}
	
	pl <- ggplot(df, aes(x = p, y = data, group = method)) +
		facet_grid(cols = vars(d), rows = vars(fstat),
		labeller = labeller(d = lab_densities, fstat = toupper),
							 scales = "free") +
		geom_line(aes(color = method)) +
		geom_point(aes(color = method)) +
		theme_bw() +
		theme(text = element_text(size = 20), legend.position = "bottom") +
		xlab("Number of nodes (p)") +
		ylab("") 
	
		pl <- pl +
			geom_ribbon(aes(ymin = data - data_se, ymax = data + data_se, fill = method),
									alpha = .2) +
			labs(fill = "Method", color = "Method")# +
			#scale_fill_discrete(labels = c("Banding", "Likelihood")) +
			#scale_color_discrete(labels = c("Banding", "Likelihood"))
	
	return(pl)
}

r <- 30 
p <- c(30, 100, 200, 500, 1000)
df <- get_statistics(p = p, r = r)
pl <- plot_sigma_exp(df)
ggplot2::ggsave(filename = "sigma_exp.pdf", plot = pl, device = "pdf", width = 11, height = 9,
								path = "../sparsecholeskycovariance/img/")
	
