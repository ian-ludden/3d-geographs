#!/bin/bash

# Timestamp function
timestamp() {
	date +"%D %T"
}

echo "Start of generate_all_inputs.sh"
echo "$(timestamp)"

# Make directories, if missing
for type in grid bcc
do
	mkdir -p $type
done

# Generate all CH (cubic honeycomb, 
# that is, uniform grid of unit cubes) 
# and BCH (bitruncated cubic honeycomb, 
# from body-centered cubic lattice) inputs
for sidelength in {5,10,15,20,25,30,35,40}
do
	type="grid"
	fname="${type}/${type}_${sidelength}.csv"
	echo $fname
	./uniform_grid.exe -s $sidelength -o $fname
	
	type="bcc"
	fname="${type}/${type}_${sidelength}.csv"
	echo $fname
	./bcc_lattice.exe -s $sidelength -o $fname
done

echo "End of generate_all_inputs.sh"
echo "$(timestamp)"
