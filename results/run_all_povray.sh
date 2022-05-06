#!/bin/bash

mkdir -p img

for pov_fname in *.pov
do
	png_fname="$(basename -s .pov ${pov_fname}).png"
	pvengine64.exe /EXIT /RENDER ${pov_fname}
	mv ${png_fname} img/
done
