args = commandArgs(trailingOnly=TRUE)

source("exp_sim.R")

nodes <- c(30, 100, 200, 500, 1000)
execute(r = args[1], ename = "l_exp", emethod = l_exp, m = "lasso", nodes = nodes, n = 200, ntrain = 100)
