#!/bin/bash

# Timestamp function
timestamp() {
	date +"%D %T"
}

# Run all CH and BCH instances with K = 2 zones
K=2

for fname in ../generate_input/grid/* 
do
	echo $fname >> ../results/local_search_grid.out
	echo "K=$K" >> ../results/local_search_grid.out
	echo "$(timestamp)" >> ../results/local_search_grid.out
	#./local_search.exe $fname $K >> ../results/local_search_grid.out
done

echo "" >> ../results/local_search_grid.out

echo "$(timestamp)"

for fname in ../generate_input/bcc/* 
do
	echo $fname >> ../results/local_search_bcc.out
	echo "K=$K" >> ../results/local_search_bcc.out
	echo "$(timestamp)" >> ../results/local_search_bcc.out
	#./local_search.exe $fname $K >> ../results/local_search_bcc.out
done

echo "" >> ../results/local_search_bcc.out

echo "$(timestamp)"

# Repeat all instances for K = 25 zones
K=25

for fname in ../generate_input/grid/* 
do
	echo $fname >> ../results/local_search_grid.out
	echo "K=$K" >> ../results/local_search_grid.out
	echo "$(timestamp)" >> ../results/local_search_grid.out
	#./local_search.exe $fname $K >> ../results/local_search_grid.out
done

echo "" >> ../results/local_search_grid.out

echo "$(timestamp)"

for fname in ../generate_input/bcc/* 
do
	echo $fname >> ../results/local_search_bcc.out
	echo "K=$K" >> ../results/local_search_bcc.out
	echo "$(timestamp)" >> ../results/local_search_bcc.out
	./local_search.exe $fname $K >> ../results/local_search_bcc.out
done

echo "" >> ../results/local_search_bcc.out

echo "$(timestamp)"
