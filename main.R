args = commandArgs(trailingOnly=TRUE)

source("gradient.R")


execute(r = args[1], ename = "sigma_exp", emethod = sigma_exp_sparse, nodes = nodes, n = 200, ntrain = 100)
