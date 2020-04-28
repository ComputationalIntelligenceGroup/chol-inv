exp <- "robot" # one of sonar or robot

filename <- c(sonar = "sonar.all-data",
			  robot = "sensor_readings_24.data")

data <- read.table(paste0("data/", filename[exp]), sep = ",", header = FALSE)
	
class <- ncol(data)
type <- levels(data[, class])
