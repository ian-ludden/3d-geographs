# Makefile for input generators
include ../config.mk

EXECS=deterministic_points uniform_grid random_points_box

all: $(EXECS)

deterministic_points: deterministic_points.cc
	$(COMPILER) $(CFLAGS) $(INC) $(LIB) -o deterministic_points deterministic_points.cc -lvoro++

uniform_grid: uniform_grid.cc
	$(COMPILER) $(CFLAGS) $(INC) $(LIB) -o uniform_grid uniform_grid.cc -lvoro++

random_points_box: random_points_box.cc
	$(COMPILER) $(CFLAGS) $(INC) $(LIB) -o random_points_box random_points_box.cc -lvoro++

clean:
	rm -f $(EXECS)

.PHONY: all clean
