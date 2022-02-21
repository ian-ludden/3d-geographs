#!/bin/bash
echo '$0 = ' $0
echo '$1 = ' $1
filepath=$1

K=8

for fname in ../generate_input/grid/* 
do
	echo $fname >> ../results/local_search_grid.out
	echo "K=8" >> ../results/local_search_grid.out
	./local_search.exe $fname $K >> ../results/local_search_grid.out
	# echo "10k flips,K=8,attempt_flip_BFS" >> ../results/local_search_grid.out
	# ./local_search_BFS.exe $fname $K >> ../results/local_search_grid.out
	
done

for fname in ../generate_input/bcc/* 
do
	echo $fname >> ../results/local_search_bcc.out
	echo "K=8" >> ../results/local_search_bcc.out
	./local_search.exe $fname $K >> ../results/local_search_bcc.out
	# echo "10k flips,K=8,attempt_flip_BFS" >> ../results/local_search_bcc.out
	# ./local_search_BFS.exe $fname $K >> ../results/local_search_bcc.out
done
