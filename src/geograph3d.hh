/** 
 * \file geograph3d.hh
 * \brief Header file for geograph3d. */

#ifndef GEOGRAPH3D_HH
#define GEOGRAPH3D_HH
#include "staticgraph.hh"
#include <string>
#include <vector>
using std::string;
using std::vector;

namespace gg3d {
/** \brief Class representing a complete 3-D geo-graph. 
 * 
 * This is a class representing a 3-D geo-graph, including:
 * - An undirected graph of cells with 2-D face adjacencies
 * - For each cell, an undirected graph of its surface dual
 * - For each cell, the induced subgraph of its augmented neighborhood
 */
class geograph3d {
public:
    /** Number of cells */
    int N;
    /** A map of cell IDs to part assignments */
    vector<int> assignment;
    /** 
     * The undirected graph of cells with 2-D face adjacencies, 
     * i.e., the dual graph of the cell complex 
     * with the external dummy vertex omitted. 
     */
    static_graph g;
    /** The cells' surface dual graphs, indexed by cell ID */
    vector<static_graph> surface_dual;
    /** The induced subgraphs of the cells' augmented neighborhoods, indexed by cell ID */
    vector<static_graph> aug_neighbor_graph;

    geograph3d(const string in_filename);
};
}

#endif /* GEOGRAPH3D_HH */
