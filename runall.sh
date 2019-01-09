#!/bin/bash

echo p = 5, N = 10, s = 0.25
Rscript main.R data/sim/data_5_10.rds data/sim/cov_5_10.rds

echo p = 10, N = 100, s = 0.25
Rscript main.R data/sim/data_10_100.rds data/sim/cov_10_100.rds

echo p = 100, N = 100, s = 0.25
Rscript main.R data/sim/data_100_100.rds data/sim/cov_100_100.rds

echo p = 100, N = 500, s = 0.25
Rscript main.R data/sim/data_100_500.rds data/sim/cov_100_500.rds

echo p = 100, N = 1000, s = 0.25
Rscript main.R data/sim/data_100_1000.rds data/sim/cov_100_1000.rds

echo p = 1000, N = 100, s = 0.25
Rscript main.R data/sim/data_1000_100.rds data/sim/cov_1000_100.rds

