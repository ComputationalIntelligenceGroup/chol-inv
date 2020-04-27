source("exp_lib.R")

sonar_data <- read.table("data/sonar.all-data", sep = ",", header = FALSE)
gen_chols <- function(data, dirname) {
	type <- levels(data$V61)

	for (t in type) {
		X <- data[data$V61 == t,][, -61]
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
	for (n in 1:nrow(data)) {
		data_train <- data[-n, ]
		preds[n, "true"] <- data[n, 61]

		train_r <- data_train[data_train$V61 == "R", ][, -61]
		train_m <- data_train[data_train$V61 == "M", ][, -61]

		p_r <- nrow(train_r) / nrow(data_train)
		p_m <- nrow(train_m) / nrow(data_train)

		mean_r <- colMeans(train_r)
		mean_m <- colMeans(train_m)
		
		subdirname <- paste0(dirname, n, "/")
		dir.create(subdirname, showWarnings = FALSE)
		gen_chols(data = data_train, dirname = subdirname)
		for (m in names(f_chol)) {
			l_r <- readRDS(file = paste0(dirname, n, "/",  m, "_R.rds"))
			l_m <- readRDS(file = paste0(dirname, n, "/",  m, "_M.rds"))
			
			x <- as.matrix(data[n, ][, -61])
			h_r <- forwardsolve(l_r, t(x - mean_r))
			h_m <- forwardsolve(l_m, t(x - mean_m))
			ll_r <-  -sum(log(diag(l_r)))  - 0.5 * crossprod(h_r) + log(p_r)
			ll_m <-  -sum(log(diag(l_m)))  - 0.5 * crossprod(h_m) + log(p_m)
 
			if (ll_r < ll_m) {
				preds[n, m] <- 1 # Mine
	 		} else {
				preds[n, m] <- 2 # Rock
			}
		}
	}

	return(preds)
}


### Cholesky factor for heatmaps
#dirname <- "sonar_heat_exp/"
#dir.create(dirname, showWarnings = FALSE)
#gen_chols(data = sonar_data, dirname = dirname)

### Prediction: small sample size --> leave one out CV instead of train/test
dirname <- "sonar_pred_exp/"
dir.create(dirname, showWarnings = FALSE)
preds <- prediction(data = sonar_data)
saveRDS(preds, file = "data/preds.rds")

