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
private:
    /** Number of cells */
    int N;
    /** Number of parts */
    const int K;
    /** A map of cell IDs to part assignments, indexed 1 to K */
    vector<int> assignment;

public:
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

    /** Constructor for 3-D geo-graph. 
     * Builds representation of 3-D geo-graph from 
     * output of Voro++ representing a Voronoi tesselation 
     * inside of a 3-D box (rectangular prism). 
     * 
     * In particular, the constructor: 
     *  - Builds the cell adjacency graph. 
     *  - Builds dual graphs of each cell's surface. 
     *  - Builds augmented neighborhood induced subgraphs for each cell. 
     * 
     * \param[in] (in_filename) String name of the input CSV file, 
     *                          such as that produced by find_augmented_neighbors.
     *                          NB: This constructor assumes rows are sorted by cell ID, ascending. 
     * \param[in] (num_parts) The number of parts in the partition (stored as K). 
     */
    geograph3d(const string in_filename, const int num_parts);

    /**
     * Setter for assignment member variable. 
     * 
     * \param[in] (assignment) New assignments, as a vector of part IDs (1 to K) 
     *                         indexed by cell ID (0 to N-1)
     */
    void set_assignment(vector<int> assignment) {
        this->assignment = assignment;
    }

    /** Getter for assignment member variable. */
    vector<int> get_assignment() { return assignment; }

    /** Returns the number of parts in the partition (i.e., member variable K). */
    int num_parts() { return K; }

    /** Returns the number of cells in the graph (i.e., member variable N). */
    int num_cells() { return N; }
};
}

#endif /* GEOGRAPH3D_HH */
