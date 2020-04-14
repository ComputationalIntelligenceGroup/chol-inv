devtools::install_github("irenecrsn/covchol")

data <- read.table("data/sonar.all-data", sep = ",", header = FALSE)

type <- c("R", "M")

pdf(file = "../sparsecholeskycovariance/img/sonar.pdf")

for (t in type) {
	X <- data[data$V61 == t,][, -61]
	ntrain <- floor(nrow(X)/2)
	Xtrain <- X[1:ntrain, ]
	Xtest <- X[(ntrain + 1):nrow(X), ]

	Covtest <- cov(Xtest)
	filled.contour((Covtest[60:1,]), color.palette = rainbow, asp = 1,
		plot.title = paste0("Sample covariance test data. Type = ", t))
	
	Lest <- gradient_est_f(ntrain = ntrain, X)
	Sigmaest <- Lest %*% t(Lest)
	filled.contour((Sigmaest[60:1,]), color.palette = rainbow, asp = 1,
		plot.title = paste0("Estimated covariance matrix. Type = ", t))

	means <- colMeans(Xtest)
}


#rng <- range(c(S_rocks, S_mines))
#filled.contour((S_mines[60:1,]) - S_rocks[60:1,], 
#               color.palette = rainbow, asp = 1 )

############ PREDICTION


######### WE HAVE TO DO CROSS_VAL

#S_r <- S_rocks
#S_m <- S_mines
#O_r <- solve(S_r)
#O_m <- solve(S_m)
#p_m <- nrow(mines) / (nrow(mines) + nrow(rocks)) 
#p_r <- 1 - p_m
#pred <- apply(mines[,-61], MARGIN = 1, function(x){
#  ll_rocks <-   -sum(log(diag(res_rocks$L)))  - 0.5 * 
#                t(x - m_rocks) %*% O_r %*% (x - m_rocks)  + log(p_r)
#  
#  ll_mines <-  -sum(log(diag(res_mines$L)))  - 0.5 * 
#    t(x - m_mines) %*% O_m %*% (x - m_mines)  + log(p_m)
#  print(ll_rocks - ll_mines)
#  if (ll_rocks < ll_mines) return("M")
#  else return("R")
#})

#sum(pred == "R")
#sum(pred == "M")
