#!/bin/bash
K=8

for fname in ../generate_input/grid/* 
do
	echo $fname >> ../results/local_search_grid.out
	echo "K=8" >> ../results/local_search_grid.out
	./local_search.exe $fname $K >> ../results/local_search_grid.out
	
done

for fname in ../generate_input/bcc/* 
do
	echo $fname >> ../results/local_search_bcc.out
	echo "K=8" >> ../results/local_search_bcc.out
	./local_search.exe $fname $K >> ../results/local_search_bcc.out
done
