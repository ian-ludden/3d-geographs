# Makefile for input generators
include ../config.mk

EXECS=geograph3d

all: $(EXECS)

geograph3d: geograph3d.cc
	$(COMPILER) $(CFLAGS) $(INC) $(LIB) -o geograph3d geograph3d.cc static_graph.cc ../generate_input/vector_utils.cc

clean:
	rm -f $(EXECS)

.PHONY: all clean
