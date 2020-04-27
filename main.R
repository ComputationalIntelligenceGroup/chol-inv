args = commandArgs(trailingOnly=TRUE)

source("exp_sim.R")

nodes <- c(30, 100, 200, 500)
rothman_exp(repetition = args[1], m = "grad_lik", nodes = nodes, n = 200, ntrain = 100)
