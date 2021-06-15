# Makefile for input generators
include ../config.mk

EXECS=uniform_grid random_points_box find_augmented_neighbors

all: $(EXECS)

uniform_grid: uniform_grid.cc
	$(COMPILER) $(CFLAGS_VORO) $(INC_VORO) $(LIB_VORO) -o uniform_grid uniform_grid.cc -lvoro++

random_points_box: random_points_box.cc
	$(COMPILER) $(CFLAGS_VORO) $(INC_VORO) $(LIB_VORO) -o random_points_box random_points_box.cc -lvoro++

# Requires vector_utils to compile
find_augmented_neighbors: find_augmented_neighbors.cc
	$(COMPILER) $(CFLAGS) $(INC_VORO) $(LIB_VORO) -o find_augmented_neighbors find_augmented_neighbors.cc vector_utils.cc

clean:
	rm -f $(EXECS)

.PHONY: all clean
