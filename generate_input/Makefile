# Makefile for input generators
include ../config.mk

EXECS=uniform_grid random_points_box #conversion_utils

all: $(EXECS)

uniform_grid: uniform_grid.cc
	$(COMPILER) $(CFLAGS_VORO) $(INC_VORO) $(LIB_VORO) -o uniform_grid uniform_grid.cc -lvoro++

random_points_box: random_points_box.cc
	$(COMPILER) $(CFLAGS_VORO) $(INC_VORO) $(LIB_VORO) -o random_points_box random_points_box.cc -lvoro++

# # Requires vector_utils to compile
# conversion_utils: conversion_utils.cc
# 	$(COMPILER) $(CFLAGS) $(INC_VORO) $(LIB_VORO) -o conversion_utils conversion_utils.cc vector_utils.cc

clean:
	rm -f $(EXECS)

.PHONY: all clean
