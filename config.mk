### Global configurations for Makefiles ###

# C++ compiler
COMPILER=g++

# Compiler flags
CFLAGS_VORO=-Wall -ansi -pedantic -O3
CFLAGS=-Wall -std=c++17 -pedantic -O3

# Include and library paths for Voro++
INC_VORO=-I/usr/local/include/voro++
LIB_VORO=-L/usr/local/lib
INC=-I../generate_input
LIB=-L../generate_input
