library("ggplot2")
library("dplyr")

stat_norm <- function(sigmatrue, sigmaest) {
	return(log(norm(sigmaest - sigmatrue)))
}
stat_f1 <- function(sigmatrue, sigmaest) {
	p <- ncol(sigmatrue)
	tp <- sum(sigmatrue != 0 & sigmaest != 0) - p
	fn <- sum(sigmatrue != 0 & sigmaest == 0)
	fp <- sum(sigmatrue == 0 & sigmaest != 0)
	return(2*tp/(2*tp + fp + fn))
}

fstat <- c("norm" = stat_norm,
		   "f1" = stat_f1)

