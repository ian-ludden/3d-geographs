/** 
 * \file geograph3d.cc
 * \brief Implementation of the 3-D geo-graph. */
#include "geograph3d.hh"
#include "../generate_input/find_augmented_neighbors.hh"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using std::string;
using std::vector;

namespace gg3d {
/** Constructor for 3-D geo-graph. 
 * Builds representation of 3-D geo-graph from 
 * output of Voro++ representing a Voronoi tesselation 
 * inside of a 3-D box (rectangular prism). 
 * 
 * In particular, the constructor: 
 *  - Builds the cell adjacency graph. 
 *  - Builds dual graphs of each cell's surface. 
 *  - Builds augmented neighborhood graphs for each cell. 
 * 
 * \param[in] (in_filename) String name of the input CSV file, 
 *                          such as that produced by find_augmented_neighbors.
 * */
geograph3d::geograph3d(const string in_filename) {
    csv_row row;
    std::ifstream in_file(in_filename);
    in_file >> row; // Read first row to get number of cells

    N = stoi(row[0]);
}
}