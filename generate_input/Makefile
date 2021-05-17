# Makefile for input generators
include ../config.mk

EXECS=deterministic_points

all: $(EXECS)

deterministic_points: deterministic_points.cc
	$(COMPILER) $(CFLAGS) $(INC) $(LIB) -o deterministic_points deterministic_points.cc -lvoro++

clean:
	rm -f $(EXECS)

.PHONY: all clean
