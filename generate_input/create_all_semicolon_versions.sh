#!/bin/bash

for type in grid bcc
do
	for fname in ${type}/*
	do
		newfname="../tmp/$(basename -s .csv ${fname})_sc.csv"
		echo "${fname} -> ${newfname}"
		python convert_commas_to_semicolons.py $fname $newfname
	done
done