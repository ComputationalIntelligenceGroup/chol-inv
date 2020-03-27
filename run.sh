#! /bin/bash

for i in {21..50}; do 
tmux new-session -d "Rscript main.R ${i}"
done 
