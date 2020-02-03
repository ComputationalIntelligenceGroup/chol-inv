library("ggplot2")
library("dplyr")

r <- 10
p <- c(10, 30, 50, 100)

stat_tpr <- function(ltrue, lest) {
	p <- ncol(ltrue)
	return((sum(ltrue != 0 & lest != 0) - p)/(sum(ltrue != 0) - p))
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

get_statistics <- function(p, r, ename, fstat) {
	data <- array(
		dim = c(length(p), 3, length(ename)),
		dimnames = list(p = p, d = 1:3, ename = ename)
	)
	data_sd <- array(
		dim = c(length(p), 3, length(ename)),
		dimnames = list(p = p, d = 1:3, ename = ename)
	)
	
	for (i in 1:length(p)) {
		stat_res <- numeric(0)
		d <- c(1/p[i], 2/p[i], 3/p[i])
		for (j in seq_along(d)) {
			for (m in ename) {
				for (k in 1:r) {
					all_res <- readRDS(file = paste0(m, "/", p[i], "_", d[j], "_r", k, ".rds"))
					stat_res[k] <- fstat(all_res$ltrue, all_res$lest) 
				}
				data[i, j, m] <- mean(stat_res)
				data_sd[i, j, m] <- stats::sd(stat_res)
			}
		}
	}
	
	df <- data %>% as.tbl_cube(met_name = "data") %>% as_tibble()
	df$ename <- as.factor(df$ename)
	df_sd <- data_sd %>% as.tbl_cube(met_name = "data_sd") %>% as_tibble()
	df$data_sd <- df_sd$data_sd
	
	return(df)
}

plot_comparison <- function(df, show_sd = TRUE, plot_title = "", plot_ylab = "", ename) {
	pl <- ggplot(df, aes(x = p, y = data, group = interaction(as.factor(d), ename))) +
		geom_line(aes(color = as.factor(d))) +
		geom_point(aes(color = as.factor(d))) +
		theme(text = element_text(size = 20)) +
		xlab("Number of nodes") +
		ylab(plot_ylab) +
		ggtitle(plot_title) +
		labs(color = "Density") +
		scale_color_discrete(labels = c("1/p", "2/p", "3/p"))
	
	if (show_sd == TRUE) {
		pl <- pl +
			geom_ribbon(aes(ymin = data - data_sd, ymax = data + data_sd, fill = ename),
									alpha = .2) +
			labs(fill = "Method") +
			scale_fill_discrete(labels = c("Banding", "Likelihood"))
	}
	
	
	return(pl)
}

df <- get_statistics(p = p, r = r, ename = c("sparse_chol", "band_chol"), fstat = stat_tpr)
pl <- plot_comparison(df, ename = c("sparse_chol", "band_chol"), plot_ylab = "True positive rate")
ggplot2::ggsave(filename = "tpr.pdf", plot = pl, device = "pdf", width = 7, height = 5)

df <- get_statistics(p = p, r = r, ename = c("sparse_chol", "band_chol"), fstat = stat_acc)
pl <- plot_comparison(df, ename = c("sparse_chol", "band_chol"), plot_ylab = "Accuracy")
ggplot2::ggsave(filename = "acc.pdf", plot = pl, device = "pdf", width = 7, height = 5)

df <- get_statistics(p = p, r = r, ename = c("sparse_chol", "band_chol"), fstat = stat_f1)
pl <- plot_comparison(df, ename = c("sparse_chol", "band_chol"), plot_ylab = "F1 score")
ggplot2::ggsave(filename = "f1.pdf", plot = pl, device = "pdf", width = 7, height = 5)


################################################

	