#!/bin/bash
echo '$0 = ' $0
echo '$1 = ' $1
filepath=$1

K=8

for fname in ../generate_input/grid/* 
do
	echo $fname >> local_search_grid.out
	echo "10k flips,K=8,attempt_flip" >> local_search_grid.out
	./local_search.exe $fname $K >> local_search_grid.out
	echo "10k flips,K=8,attempt_flip_BFS" >> local_search_grid.out
	./local_search_BFS.exe $fname $K >> local_search_grid.out
	
done

for fname in ../generate_input/bcc/* 
do
	echo $fname >> local_search_bcc.out
	echo "10k flips,K=8,attempt_flip" >> local_search_bcc.out
	./local_search.exe $fname $K >> local_search_bcc.out
	echo "10k flips,K=8,attempt_flip_BFS" >> local_search_bcc.out
	./local_search_BFS.exe $fname $K >> local_search_bcc.out
done
