#!/bin/bash

echo simulating p = 5, N = 10, s = 25 
Rscript sim.R 5 0.25 10 data_5_10.rds cov_5_10.rds

echo simulating p = 10, N = 100, s = 25 
Rscript sim.R 10 0.25 100 data_10_100.rds cov_10_100.rds 

echo simulating p = 100, N = 100, s = 25 
Rscript sim.R 100 0.25 100 data_100_100.rds cov_100_100.rds

echo simulating p = 100, N = 500, s = 25 
Rscript sim.R 100 0.25 500 data_100_500.rds cov_100_500.rds

echo simulating p = 100, N = 1000, s = 25 
Rscript sim.R 100 0.25 1000 data_100_1000.rds cov_100_1000.rds

echo simulating p = 1000, N = 100, s = 25 
Rscript sim.R 1000 0.25 100 data_1000_100.rds cov_1000_100.rds


mkdir -p sim/

mv *.rds sim/

