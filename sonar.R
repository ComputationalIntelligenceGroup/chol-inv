source("lib.R")

devtools::install_github("irenecrsn/covchol")

type <- c("R", "M")

cov_sparse <- function(ntrain, X) {
	Lest <- gradient_est(ntrain = ntrain, X)
	return(Lest %*% t(Lest))
}

cov_sparse_f <- function(ntrain, X) {
	Lest <- gradient_est_f(ntrain = ntrain, X)
	return(Lest %*% t(Lest))
}

cov_sample <- function(ntrain, X) {
	return(cov(X))
}

cov_band <- function(ntrain, X) {
	return(band_est(ntrain = ntrain, X = X))
}

cov_lasso <- function(ntrain, X) {
	Lest <- lasso_est(ntrain = ntrain, X)
	return(Lest %*% t(Lest))
}

f_cov <- c("grad_lik" = cov_sparse,
			"grad_frob" = cov_sparse_f,
			"sample" = cov_sample,
			"band" = cov_band,
			"lasso" = cov_lasso)

get_covs <- function(data, est) {

	covs <- array(dim = c(length(type), length(est), 60, 60),
			dimnames = list(type = type, est = est, rows = 1:60, cols =
			1:60))

	for (t in type) {
		X <- data[data$V61 == t,][, -61]
		ntrain <- floor(nrow(X)/2)
		
		for (e in est) {
			covs[t, e, ,] <- f_cov[[e]](ntrain, X)
		}
	}

	return(covs)
}

### Whole covariances
data <- read.table("data/sonar.all-data", sep = ",", header = FALSE)
est <- c("sample", "grad_lik", "grad_frob", "band")
covs <- get_covs(data = data, est = est)
saveRDS(covs, file = "data/covs.rds")

### Prediction
data <- read.table("data/sonar.all-data", sep = ",", header = FALSE)
ntrain <- floor(nrow(data)/2)
data <- data[sample.int(nrow(data)), ]  # Shuffle
data_train <- data[1:ntrain, ]
data_test <- data[(ntrain + 1):nrow(data), ]
train_r <- data_train[data_train$V61 == "R", ]
train_m <- data_train[data_train$V61 == "M", ]

est <- c("grad_frob", "band") # sparse and sample yield non positive definite
covs <- get_covs(data_train, est)

p_r <- nrow(train_r) / nrow(data_train)
p_m <- nrow(train_m) / nrow(data_train)

mean_r <- colMeans(train_r[, -61])
mean_m <- colMeans(train_m[, -61])

preds <- matrix(nrow = nrow(data_test), ncol = length(est) + 1,
			dimnames = list(sample = 1:nrow(data_test), est = c(est, "true")))
preds[, "true"] <- data_test[, 61]

for (e in est) {
	O_r <- solve(covs["R", e, ,])
	O_m <- solve(covs["M", e, ,])
	
	preds[, e] <- apply(data_test[,-61], MARGIN = 1, function(x){
		ll_rocks <-  -0.5*log(det(O_r))  - 0.5 * 
        	       t(x - mean_r) %*% O_r %*% (x - mean_r)  + log(p_r)
 
		ll_mines <-  -0.5*log(det(O_m))  - 0.5 * 
  					t(x - mean_m) %*% O_m %*% (x - mean_m)  + log(p_m)
		if (ll_rocks < ll_mines) {
			return("M") 
	 	} else {
			return("R")
		}
	})
}

saveRDS(preds, file = "data/preds.rds")

