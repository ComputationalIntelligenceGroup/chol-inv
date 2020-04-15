source("lib.R")

devtools::install_github("irenecrsn/covchol")

type <- c("R", "M")
est <- c("sample", "sparse", "sparse_f", "band")

get_covs <- function(data) {
	covs <- array(dim = c(length(type), length(est), 60, 60),
			dimnames = list(type = type, est = est, rows = 1:60, cols = 1:60))

	for (t in type) {
		X <- data[data$V61 == t,][, -61]
		ntrain <- floor(nrow(X)/2)
		Xtrain <- X[1:ntrain, ]
		Xtest <- X[(ntrain + 1):nrow(X), ]

		covs[t, "sample", ,] <- cov(Xtest)
	
		Lest <- gradient_est(ntrain = ntrain, X)
		covs[t, "sparse", ,] <- Lest %*% t(Lest)
	
		Lest <- gradient_est_f(ntrain = ntrain, X)
		covs[t, "sparse_f", ,] <- Lest %*% t(Lest)

		#Lest <- lasso_est(ntrain = ntrain, X)
		#covs[t, "lasso", ,] <- Lest %*% t(Lest)

		covs[t, "band", , ] <- band_est(ntrain = ntrain, X)
	}

	return(covs)
}

### Whole covariances
#data <- read.table("data/sonar.all-data", sep = ",", header = FALSE)
#covs <- get_covs(data)
#saveRDS(covs, file = "data/covs.rds")

### Prediction
data <- read.table("data/sonar.all-data", sep = ",", header = FALSE)
ntrain <- floor(nrow(data)/2)
data_train <- data[1:ntrain, ]
data_test <- data[(ntrain + 1):nrow(data), ]
#covs <- get_covs(data_train)
#saveRDS(covs, file = "data/covs_preds.rds")

covs <- readRDS(file = "data/covs_preds.rds")
train_r <- data_train[data_train$V61 == "R", ]
train_m <- data_train[data_train$V61 == "M", ]

p_r <- nrow(train_r) / nrow(data_train)
p_m <- nrow(train_m) / nrow(data_train)

mean_r <- colMeans(train_r[, -61])
mean_m <- colMeans(train_m[, -61])

preds <- matrix(nrow = nrow(data_test), ncol = length(est) + 1,
			dimnames = list(sample = 1:nrow(data_test), est = c(est, "true")))
preds[, "true"] <- data_test[, 61]

est <- c("sparse_f", "band")
for (e in est) {
	O_r <- solve(covs["R", e, ,])
	O_m <- solve(covs["M", e, ,])
	
	preds[, e] <- apply(data_test[,-61], MARGIN = 1, function(x){
		ll_rocks <-  -0.5*log(det(O_r))  - 0.5 * 
        	       t(x - mean_r) %*% O_r %*% (x - mean_r)  + log(p_r)
 
		ll_mines <-  -0.5*log(det(O_m))  - 0.5 * 
  					t(x - mean_m) %*% O_m %*% (x - mean_m)  + log(p_m)
		if (ll_rocks < ll_mines) return("M")
	 	else return("R")
	})
}

saveRDS(preds, file = "data/preds.rds")

