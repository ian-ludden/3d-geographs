#!/bin/bash

# Timestamp function
timestamp() {
	date +"%D %T"
}

# Write CSV headers
for type in grid bcc
do
	echo "number of zones,number of cells,flip type,status,average microseconds,number of samples" >> ../results/data_${type}.csv
done

# Run random local search experiments and write results
for K in {2,8,25}
do
	for type in grid bcc
	do
		echo "Start K=${K} for ${type}: $(timestamp)"
		for fname in ../generate_input/${type}/*
		do
			fname_short="$(basename -s .csv ${fname})"
			./local_search.exe $fname $K --partition ../results/partition_${fname_short}_K${K}.csv >> ../results/data_${type}.csv
		done
		echo "End K=${K} for ${type}: $(timestamp)"
		echo ""
	done
done

