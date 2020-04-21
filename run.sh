#! /bin/bash

for i in {1..10}; do 
tmux new-session -d "Rscript main.R ${i}"
done 
