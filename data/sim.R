library("MASS")
library("Matrix")

# Basic initialization
cmd_args <- commandArgs(trailingOnly = TRUE)

if (length(cmd_args) == 5) {
	cov_dim <- as.integer(cmd_args[[1]])
	cov_sparsity <- as.double(cmd_args[[2]])
	sample_size <- as.integer(cmd_args[[3]])
	data_file <- cmd_args[[4]]
	cov_file <- cmd_args[[5]]
} 

### Generate a random sparse covariance matrix
## Generate random sparse cholesky factor 
# Strictly lower triangular sparse matrix
rsl <- Matrix::tril(Matrix::rsparsematrix(cov_dim, cov_dim, cov_sparsity), -1)
# Ones along the diagonal
diag(rsl) <- 1

## Generate random conditional covariance matrix
# Random non-sparse diagonal matrix
rsd <- Matrix::Diagonal(cov_dim, runif(cov_dim, 0, 1))

## Compute the product --> sparse covariance matrix
rsc <- rsl %*% rsd %*% t(rsl)
# Save true covariance matrix for future comparison
saveRDS(rsc, cov_file)

## Simulate from a N(0, rsc)
rdata <- MASS::mvrnorm(n = sample_size, rep(0, cov_dim), rsc)
# Save the simulated data for future estimation
saveRDS(rdata, data_file)

