library("ggplot2")
library("dplyr")

method <- c("lasso",
			"band",
			"grad_lik",
			"grad_frob")

stat_norm <- function(mtrue, mest) {
	return(log(norm(mest - mtrue)))
}
stat_f1 <- function(mtrue, mest) {
	p <- ncol(mtrue)
	tp <- sum(mtrue != 0 & mest != 0) - p
	fn <- sum(mtrue != 0 & mest == 0)
	fp <- sum(mtrue == 0 & mest != 0)
	return(2*tp/(2*tp + fp + fn))
}

fstat <- c("norm" = stat_norm,
		   "f1" = stat_f1)

fname_l_exp <- function(idx, p, m, rep, true = FALSE) {
	d <- idx/p
	if (true == TRUE) {
		return(paste0("l_exp/ltrue_", p, "_", d, "_r", rep, ".rds"))
	} else {
		return(paste0("l_exp/", m, "_", p, "_", d, "_r", rep, ".rds"))
	}
}

fname_rothman_exp <- function(idx, p, m, rep, true = FALSE) {
	if (true == TRUE) {
		return(paste0("rothman_exp/sigma", idx, "true_", p, ".rds"))
	} else {
		return(paste0("rothman_exp/sigma", idx, m, "_", p, "_r", rep, ".rds"))
	}
	return(file)
}

fname <- c("l_exp" = fname_l_exp,
			"rothman_exp" = fname_rothman_exp)

get_statistics <- function(p, r, ename) {
	data <- array(
		dim = c(length(p), 3, length(method), length(fstat)),
		dimnames = list(p = p, scenario = 1:3, method = method, fstat = names(fstat))
	)
	data_se <- array(
		dim = c(length(p), 3, length(method), length(fstat)),
		dimnames = list(p = p, scenario = 1:3, method = method, fstat = names(fstat))
	)
	stat_res <- array(
		dim = c(r, length(fstat), length(method)),
		dimnames = list(r = 1:r, fstat = names(fstat), method = method)
	)
 	error <- logical(length(method))
	names(error) <- method
	for (i in 1:length(p)) {
		for (j in 1:3) {
			for (k in 1:r) {
				for (m in method) {
					error[m] <- FALSE
					file_est <- fname[[ename]](idx = j, p = p[i], m = m, rep = k)
					if (file.exists(file_est)) {
						file_true <- fname[[ename]](idx = j, p = p[i], m = m,
													rep = k, true = TRUE)
						mest <- readRDS(file = file_est)
						mtrue <- readRDS(file = file_true)
						for (l in seq(length(fstat))) {
							stat_res[k, l, m] <- fstat[[l]](mtrue, mest) 
						}
					} else {
						error[m] <- TRUE
						stat_res[k, , m] <- NA
					}
				}
			}
			for (m in method) {
				for (l in seq(length(fstat))) {
					res <- stat_res[, l, m][!is.na(stat_res[, l, m])]

					if(length(res) == 0) {
						data[i, j, m, l] <- NA
						data_se[i, j, m, l] <- NA
					} else {
						data[i, j, m, l] <- mean(res)
						data_se[i, j, m, l] <- stats::sd(res)/sqrt(length(res))
					}
				}
			}
		}
	}
	
	df <- data %>% as.tbl_cube(met_name = "data") %>% as_tibble()	
	df$method<- as.factor(df$method)
	df_se <- data_se %>% as.tbl_cube(met_name = "data_se") %>% as_tibble()
	df$data_se <- df_se$data_se

	df <- stats::na.omit(df)
	
	return(df)
}

lab_densities <- function(str) {
	return(paste0("Density = ", str, "/p"))
}

lab_sigmas <- function(str) {
	return(paste0("Sigma", str))
}

lab <- c("l_exp" = lab_densities,
		"rothman_exp" = lab_sigmas)

plot_statistics <- function(df, ename) {
	
	pl <- ggplot(df, aes(x = p, y = data, group = method)) +
		facet_grid(cols = vars(scenario), rows = vars(fstat), 
							 labeller = labeller(fstat = toupper, 
							 scenario = lab[[ename]]),
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
for (ename in names(fname)) {
	df <- get_statistics(p = p, r = r, ename = ename)
	pl <- plot_statistics(df, ename = ename)
	ggplot2::ggsave(filename = paste0(ename, ".pdf"), plot = pl, device = "pdf", width = 11, height = 6,
								path = "../sparsecholeskycovariance/img/")
}
	
