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

get_statistics <- function(p, r) {
	fstat <- c("tpr" = stat_tpr,
						 "tnr" = stat_tnr,
						 "frob" = stat_frob,
						 "f1" = stat_f1)
	method <- c("lasso",
				"nestedlasso",
				"sparse_f")
	data <- array(
		dim = c(length(p), 3, length(method), length(fstat)),
		dimnames = list(p = p, d = 1:3, method = method, fstat = names(fstat))
	)
	data_sd <- array(
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
				error_lasso <- FALSE
				file_lasso <- paste0("l_exp/lasso_", p[i], "_", d[j], "_r", k, ".rds")
				if (file.exists(file_lasso)) {
					res_lasso <- readRDS(file = file_lasso)
				} else {
					error_lasso <- TRUE
				}
				error_nestedlasso <- FALSE
				file_nestedlasso <- paste0("l_exp/nestedlasso_", p[i], "_", d[j], "_r", k, ".rds")
				if (file.exists(file_nestedlasso)) {
					res_nestedlasso <- readRDS(file = file_nestedlasso)
				} else {
					error_nestedlasso <- TRUE
				}
				res_sparse_f <- readRDS(file = paste0("l_exp/sparse_f_", p[i], "_", d[j], "_r", k, ".rds"))
				res_true <- readRDS(file = paste0("l_exp/ltrue_", p[i], "_", d[j], "_r", k, ".rds"))
				for (l in seq(length(fstat))) {
					if (error_lasso) {
						stat_res[k, l, "lasso"] <- NA
					} else {
						stat_res[k, l, "lasso"] <- fstat[[l]](res_true, res_lasso) 
					}
					if (error_nestedlasso) {
						stat_res[k, l, "nestedlasso"] <- NA
					} else {
						stat_res[k, l, "nestedlasso"] <- fstat[[l]](res_true, res_nestedlasso) 
					}
					stat_res[k, l, "sparse_f"] <- fstat[[l]](res_true, res_sparse_f) 
				}
			}
			for (m in method) {
				for (l in seq(length(fstat))) {
					data[i, j, m, l] <- mean(stat_res[, l, m])
					data_sd[i, j, m, l] <- stats::sd(stat_res[, l, m])/sqrt(r)
				}
			}
		}
	}
	
	df <- data %>% as.tbl_cube(met_name = "data") %>% as_tibble()
	df$method<- as.factor(df$method)
	df_sd <- data_sd %>% as.tbl_cube(met_name = "data_sd") %>% as_tibble()
	df$data_sd <- df_sd$data_sd
	
	return(df)
}

plot_comparison <- function(df, plot_title = "", plot_ylab = "", method) {
	
	lab_densities <- function(str) {
		return(paste0("Density = ", str, "/p"))
	}
	
	pl <- ggplot(df, aes(x = p, y = data, group = method)) +
		facet_grid(cols = vars(d), rows = vars(fstat), 
							 labeller = labeller(fstat = toupper, d = lab_densities),
							 scales = "free") +
		geom_line(aes(color = method)) +
		geom_point(aes(color = method)) +
		theme_bw() +
		theme(text = element_text(size = 20), legend.position = "bottom") +
		xlab("Number of nodes (p)") +
		ylab("") 
	
		pl <- pl +
			geom_ribbon(aes(ymin = data - data_sd, ymax = data + data_sd, fill = method),
									alpha = .2) +
			labs(fill = "Method", color = "Method") #+
			#scale_fill_discrete(labels = c("Lasso", "Nested Lasso", "Likelihood")) +
			#scale_color_discrete(labels = c("Lasso", "Nested Lasso", "Likelihood"))
	
	return(pl)
}

r <- 50
p <- c(30, 100, 200, 500, 1000)
df <- get_statistics(p = p, r = r)
pl <- plot_comparison(df)
ggplot2::ggsave(filename = "l_exp.pdf", plot = pl, device = "pdf", width = 11, height = 9,
								path = "../sparsecholeskycovariance/img/")

	
