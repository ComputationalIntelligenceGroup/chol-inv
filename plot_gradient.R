library("ggplot2")
library("dplyr")

r <- 10
p <- c(10, 30, 50, 100)

get_statistics <- function(p, r, ename) {
	data <- array(
		dim = c(length(p), 3, length(ename)),
		dimnames = list(p = p, d = 1:3, ename = ename)
	)
	data_sd <- array(
		dim = c(length(p), 3, length(ename)),
		dimnames = list(p = p, d = 1:3, ename = ename)
	)
	
	for (i in 1:length(p)) {
		tpr <- numeric(0)
		d <- c(1/p[i], 2/p[i], 3/p[i])
		for (j in seq_along(d)) {
			for (m in ename) {
				for (k in 1:r) {
					all_res <- readRDS(file = paste0(m, "/", p[i], "_", d[j], "_r", k, ".rds"))
					tpr[k] <- (sum(all_res$ltrue != 0 & all_res$lest != 0) - p[i])/(sum(all_res$ltrue != 0) - p[i])
				}
				data[i, j, m] <- mean(tpr)
				data_sd[i, j, m] <- stats::sd(tpr)
			}
		}
	}
	
	df <- data %>% as.tbl_cube(met_name = "data") %>% as_tibble()
	df$ename <- as.factor(df$ename)
	df_sd <- data_sd %>% as.tbl_cube(met_name = "data_sd") %>% as_tibble()
	df$data_sd <- df_sd$data_sd
	
	return(df)
}

plot_comparison <- function(df, show_sd = FALSE, plot_title = "", plot_ylab = "") {
	palette <- grDevices::colorRampPalette(colors = c("green4", "blue"))
	colors <- palette(length(ename))
	
	pl <- ggplot(df, aes(x = p, y = data, group = ename)) +
		geom_line(aes(color = ename)) +
		geom_point(aes(color = ename)) +
		scale_color_manual(labels = ename, values = colors) +
		theme(text = element_text(size = 20), legend.position = "bottom") +
		xlab("Number of nodes") +
		ylab(plot_ylab) +
		ggtitle(plot_title)
	
	if (show_sd == TRUE) {
		pl <- pl +
			geom_ribbon(aes(ymin = data - data_sd, ymax = data + data_sd, fill = ename),
									alpha = .3) +
			scale_fill_manual(labels = ename, values = colors)
	}
	
	
	return(pl)
}

df <- get_statistics(p = p, r = r, ename = c("sparse_chol", "band_chol"))
df %>% knitr::kable(format = 'latex', booktabs = TRUE)

################################################

	