library("ggplot2")
library("dplyr")

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

