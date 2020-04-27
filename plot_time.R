library("ggplot2")
library("dplyr")

source("exp_lib.R")

p <- c(30, 100, 200, 500)
repetitions <- 10
dirname <- "time_exp/"

time <- array(dim = c(length(f_chol), length(p), 3, repetitions),
			dimnames = list(method = names(f_chol), p = p, d = 1:3, r =
			1:repetitions))
data <- array(dim = c(length(f_chol), length(p), 3),
			dimnames = list(method = names(f_chol), p = p, d = 1:3))
data_se <- array(dim = c(length(f_chol), length(p), 3),
			dimnames = list(method = names(f_chol), p = p, d = 1:3))

for (m in names(f_chol)) {
	for (i in 1:length(p)) {
		for (d in 1:3) {
			for (r in 1:repetitions) {
				time[m, i, d, r] <- 
					readRDS(
						file = paste0(dirname, m, "_", p[i], "_", d/p[i], 
						"_r", r, ".rds")
					)
			}
			value <- time[m, i, d, ]
			data[m, i, d] <- mean(value)
			data_se[m, i, d] <- stats::sd(value)/sqrt(length(value))
		}
	}
}

df <- data %>% as_tbl_cube(met_name = "data") %>% as_tibble()
df$method <- as.factor(df$method)
df_se <- data_se %>% as_tbl_cube(met_name = "data_se") %>% as_tibble()
df$data_se <- df_se$data_se

lab_densities <- function(str) {
	return(paste0("Density = ", str, "/p"))
}

pl <- ggplot(df, aes(x = p, y = data, group = method)) +
	facet_grid(cols = vars(d), labeller = labeller(d = lab_densities)) +
	geom_line(aes(color = method)) +
	geom_point(aes(color = method)) +
	geom_ribbon(aes(ymin = data - data_se, ymax = data +  data_se, fill =
	method) +
	theme_bw() +
	theme(text = element_text(size = 20), legend.position = "bottom") +
	xlab("Number of nodes(p)") +
	ylab("") +
	labs(fill = "Method", color = "Method")

ggsave(filename = "time.pdf", plot = pl, device = "pdf", width = 11, 
	   height = 3, path = "../sparsecholeskycovariance/img/")
