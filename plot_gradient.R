library("ggplot2")
library("dplyr")

r <- 10
p <- c(10, 30, 50, 100)
d <- c(0.0025, 0.005, 0.025, 0.05, 0.25, 0.5)
N <- 10*p
ename <- "gradient"

data <- matrix(
		nrow = length(p), ncol = length(d),
		dimnames = list(p = p, d = d)
)
data_sd <- matrix(
		nrow = length(p), ncol = length(d),
		dimnames = list(p = p, d = d)
)
	
	
	for (i in 1:length(p)) {
		exp_res <- array(dim = c(p[i], p[i], N[i] * r))
		for (j in 1:length(d)) {
			for (rep in 1:r) {
					all_res <- readRDS(file = paste0(ename, "_r", rep, "/", p[i], "_", d[j], ".rds"))
					exp_res[, , ((rep - 1) * N[i] + 1):(rep * N[i])] <- 
						(sum(all_res$ltrue != 0 & all_res$lest != 0) - p[i])/(sum(all_res$ltrue != 0) - p[i])
			}
			data[i, j] <- mean(exp_res)
			data_sd[i, j] <- stats::sd(exp_res)
		}
	}
	
	wd <- getwd()
	dir.create(paste0(wd, "/plot_", r), showWarnings = FALSE)
	
	palette <- grDevices::colorRampPalette(colors = c("black", "red"))
	colors <- palette(length(d))
	
	df <- data %>% as.tbl_cube(met_name = "data") %>% as_tibble()
	df$d <- as.factor(df$d)
	df_sd <- data_sd %>% as.tbl_cube(met_name = "data_sd") %>% as_tibble()
	df$data_sd <- df_sd$data_sd
	
	pl <- ggplot(df, aes(x = p, y = data, group = d, color = d)) +
		geom_line() +
		geom_point() +
		theme(text = element_text(size = 20), legend.position = "bottom") +
		scale_color_manual(values = colors) +
		xlab("Number of nodes") + 
		ylab("True positive rate") +
		ggtitle("")
	
	plot(pl)
	
		pl <- pl +
			geom_ribbon(aes(ymin = data - data_sd, ymax = data + data_sd, fill = ename),
									alpha = .3) +
			scale_fill_manual(labels = ename, values = colors)  
plot(pl)

	