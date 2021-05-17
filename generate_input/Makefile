# Makefile for input generators
include ../config.mk

EXECS=deterministic_points uniform_grid

all: $(EXECS)

deterministic_points: deterministic_points.cc
	$(COMPILER) $(CFLAGS) $(INC) $(LIB) -o deterministic_points deterministic_points.cc -lvoro++

uniform_grid: uniform_grid.cc
	$(COMPILER) $(CFLAGS) $(INC) $(LIB) -o uniform_grid uniform_grid.cc -lvoro++

clean:
	rm -f $(EXECS)

.PHONY: all clean
