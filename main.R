args = commandArgs(trailingOnly=TRUE)

source("exp_sim.R")

nodes <- c(1000)
execute(r = args[1], ename = "l_exp", emethod = l_exp_lasso, nodes = nodes, n = 200, ntrain = 100)
