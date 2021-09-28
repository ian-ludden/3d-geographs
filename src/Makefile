# Makefile for input generators
include ../config.mk

EXECS=local_search local_search_BFS

all: $(EXECS)

# geograph3d: geograph3d.cc
# 	$(COMPILER) $(CFLAGS) $(INC) $(LIB) -o geograph3d geograph3d.cc static_graph.cc ../generate_input/vector_utils.cc

local_search: local_search.cc
	$(COMPILER) $(CFLAGS) $(INC) $(LIB) -o local_search local_search.cc geograph3d.cc static_graph.cc ../generate_input/vector_utils.cc	../generate_input/conversion_utils.cc

local_search_BFS: local_search_BFS.cc
	$(COMPILER) $(CFLAGS) $(INC) $(LIB) -o local_search_BFS local_search_BFS.cc geograph3d.cc static_graph.cc ../generate_input/vector_utils.cc	../generate_input/conversion_utils.cc

clean:
	rm -f $(EXECS)

.PHONY: all clean
