source("exp_lib.R")

exp <- "robot" # one of sonar or robot

filename <- c(sonar = "sonar.all-data",
			  robot = "sensor_readings_24.data")

data <- read.table(paste0("data/", filename[exp]), sep = ",", header = FALSE)
gen_chols <- function(data, dirname) {
	class <- ncol(data)
	type <- levels(data[, class])

	for (t in type) {
		X <- data[data[, class] == t,][, -class]
		ntrain <- floor(nrow(X)/2)
		
		for (m in names(f_chol)) {
			saveRDS(f_chol[[m]](ntrain, X),
					file = paste0(dirname, m, "_", t, ".rds"))
		}
	}
}
prediction <- function(data) {
	preds <- matrix(nrow = nrow(data), ncol = length(f_chol) + 1,
				dimnames = list(sample = 1:nrow(data), est = c(names(f_chol), "true")))
	
	class <- ncol(data)
	
	type <- levels(data[, class])
	
	ll <- numeric(length(type))
	names(ll) <- type
	
	
	for (n in 1:nrow(data)) {
		data_train <- data[-n, ]
		x <- as.matrix(data[n, ][, -class])
		preds[n, "true"] <- data[n, class]
		
		subdirname <- paste0(dirname, n, "/")
		dir.create(subdirname, showWarnings = FALSE)
		gen_chols(data = data_train, dirname = subdirname)
		
		for (t in type) {
			train_t <- data_train[data_train[, class] == t, ][, -class]
			
			p_hat <- nrow(train_t) / nrow(data_train)

			mean_hat <- colMeans(train_t)
		
			for (m in names(f_chol)) {
			
				chol_hat <- readRDS(file = paste0(dirname, n, "/",  m, "_", t, ".rds"))
				
				h_hat <- forwardsolve(chol_hat, t(x - mean_hat))
			
				ll[t] <-  -sum(log(diag(chol_hat)))  - 0.5 * crossprod(h_hat) +
				log(p_hat)
			}
		}
	
		
		for (m in names(f_chol)) {
			ll_max <- type[1]
			preds[n, m] <- 1
			for (i in 2:length(type)) {
				t <- type[i]
				if (ll[ll_max] < ll[t]) {
					preds[n, m] <- match(t, type)
					ll_max <- t
				}
			}
		}
	}

	return(preds)
}


### Cholesky factor for heatmaps
#dirname <- paste0(exp, "_heat_exp/")
#dir.create(dirname, showWarnings = FALSE)
#gen_chols(data = data, dirname = dirname)

### Prediction: small sample size --> leave one out CV instead of train/test
dirname <- paste0(exp, "_pred_exp/")
dir.create(dirname, showWarnings = FALSE)
preds <- prediction(data = data)
saveRDS(preds, file = paste0("data/", exp, "_preds.rds"))

