source("plot_lib.R")

get_statistics <- function(p, r) {
	method <- c("grad_lik", "grad_frob", "band", "sample")
	data <- array(
		dim = c(length(p), 3, length(method), length(fstat)),
		dimnames = list(p = p, sigma = 1:3, method = method, fstat = names(fstat))
	)
	data_se <- array(
		dim = c(length(p), 3, length(method), length(fstat)),
		dimnames = list(p = p, sigma = 1:3, method = method, fstat = names(fstat))
	)
	stat_res <- array(
		dim = c(r, length(fstat), length(method)),
		dimnames = list(r = 1:r, fstat = names(fstat), method = method)
	)
		
	for (i in 1:length(p)) {
		for (sigma in 1:3) {
			sigmatrue <- readRDS(file = paste0("rothman_exp/sigma", sigma, "true_", p[i], ".rds"))
			for (k in 1:r) {
				sigmasparse <- readRDS(file = paste0("rothman_exp/sigma", sigma, "sparse_", p[i], "_r", k, ".rds"))
				sigmasparse_f <- readRDS(file = paste0("rothman_exp/sigma", sigma, "sparse_f_", p[i], "_r", k, ".rds"))
				sigmaband <- readRDS(file = paste0("rothman_exp/sigma", sigma, "band_", p[i], "_r", k, ".rds"))
				sigmasample <- readRDS(file = paste0("rothman_exp/sigma", sigma, "sample_", p[i], "_r", k, ".rds"))
				for (l in seq(length(fstat))) {
					stat_res[k, l, "grad_lik"] <- fstat[[l]](sigmatrue, sigmasparse)
					stat_res[k, l, "grad_frob"] <- fstat[[l]](sigmatrue, sigmasparse_f)
					stat_res[k, l, "band"] <- fstat[[l]](sigmatrue, sigmaband)
					stat_res[k, l, "sample"] <- fstat[[l]](sigmatrue, sigmasample)
				}
			}
			for (m in method) {
				for (l in seq(length(fstat))) {
					data[i, sigma, m, l] <- mean(stat_res[, l, m])
					data_se[i, sigma, m, l] <- stats::sd(stat_res[, l, m])/sqrt(r)
				}
			}
		}
	}
	
	df <- data %>% as.tbl_cube(met_name = "data") %>% as_tibble()
	df$method <- as.factor(df$method)
	df$sigma <- as.factor(df$sigma)
	df$fstat <- as.factor(df$fstat)
	df_se <- data_se %>% as.tbl_cube(met_name = "data_se") %>% as_tibble()
	df$data_se <- df_se$data_se
	
	return(df)
}

plot_comparison <- function(df, plot_title = "", plot_ylab = "") {
	
	lab_sigmas <- function(str) {
		return(paste0("Sigma", str))
	}
	
	pl <- ggplot(df, aes(x = p, y = data, group = method)) +
		facet_grid(cols = vars(sigma), rows = vars(fstat), 
							 labeller = labeller(fstat = toupper, sigma = lab_sigmas),
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
		labs(fill = "Method", color = "Method") 
	
	return(pl)
}


r <- 50 
p <- c(30, 100, 200, 500, 1000)
df <- get_statistics(p = p, r = r)
pl <- plot_comparison(df)
ggplot2::ggsave(filename = "rothman_exp.pdf", plot = pl, device = "pdf", width = 11, height = 6,
								path = "../sparsecholeskycovariance/img/")

