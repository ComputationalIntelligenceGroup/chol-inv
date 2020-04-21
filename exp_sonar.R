source("exp_lib.R")

devtools::install_github("irenecrsn/covchol")

type <- c("R", "M")

get_covs <- function(data, est) {

	covs <- array(dim = c(length(type), length(est), 60, 60),
			dimnames = list(type = type, est = est, rows = 1:60, cols =
			1:60))

	for (t in type) {
		X <- data[data$V61 == t,][, -61]
		ntrain <- floor(nrow(X)/2)
		
		for (e in est) {
			covs[t, e, ,] <- f_chol[[e]](ntrain, X)
		}
	}

	return(covs)
}

### Whole covariances
#data <- read.table("data/sonar.all-data", sep = ",", header = FALSE)
#covs <- get_covs(data = data, est = names(f_cov))
#saveRDS(covs, file = "data/covs.rds")

### Prediction: small sample size --> leave one out CV instead of train/test
dirname <- "sonar_exp/"

# Covariance generation
get_covs_sonar <- function(data) {
	dir.create(dirname, showWarnings = FALSE)

	for (n in 1:nrow(data)) {
		data_train <- data[-n, ]

		for (t in type) {
			X <- data_train[data_train$V61 == t,][, -61]
			ntrain <- floor(nrow(X)/2)
		
			for (m in names(f_cov)) {
				L <- f_chol[[m]](ntrain, X)
				saveRDS(L, file = paste0(dirname, m, "_", t, "_", n, ".rds"))
			}
		}
	}
}
data <- read.table("data/sonar.all-data", sep = ",", header = FALSE)
get_covs_sonar(data)

preds <- matrix(nrow = nrow(data), ncol = length(f_cov) + 1,
			dimnames = list(sample = 1:nrow(data), est = c(names(f_cov), "true")))
for (n in 1:nrow(data)) {
	data_train <- data[-n, ]
	preds[n, "true"] <- data[n, 61]
	train_r <- data_train[data_train$V61 == "R", ][, -61]
	train_m <- data_train[data_train$V61 == "M", ][, -61]

	p_r <- nrow(train_r) / nrow(data_train)
	p_m <- nrow(train_m) / nrow(data_train)

	mean_r <- colMeans(train_r)
	mean_m <- colMeans(train_m)

	for (m in names(f_cov)) {
		tryCatch(
		{
			l_r <- readRDS(file = paste0(dirname, m, "_R_", n, ".rds"))
			l_m <- readRDS(file = paste0(dirname, m, "_M_", n, ".rds"))
			
			x <- as.matrix(data[n, ][, -61])
			h_r <- forwardsolve(l_r, t(x - mean_r))
			h_m <- forwardsolve(l_m, t(x - mean_m))
			ll_r <-  -sum(log(diag(l_r)))  - 0.5 * crossprod(h_r)  + log(p_r)
			ll_m <-  -sum(log(diag(l_m)))  - 0.5 * crossprod(h_m)  + log(p_m)
 
			if (ll_r < ll_m) {
				preds[n, m] <- "M"
	 		} else {
				preds[n, m] <- "R"
			}
		},
		error = function(error_type) {
			preds[n, m] <- NA
		}
		)
	}
}

saveRDS(preds, file = "data/preds.rds")

