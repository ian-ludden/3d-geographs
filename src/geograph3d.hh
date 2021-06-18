/** 
 * \file geograph3d.hh
 * \brief Header file for geograph3d. */

#ifndef GEOGRAPH3D_HH
#define GEOGRAPH3D_HH
#include "static_graph.hh"
#include <iostream>
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
    /** List of lists of neighbor IDs; primary index is cell ID. Includes wall neighbors (-1 through -6)*/
    vector<vector<int>> neighbors;
    /** List of lists of augmented neighbor (but not neighbor) IDs; primary index is cell ID */
    vector<vector<int>> aug_neighbors;
    /** 
     * List of lists of faces, stored as uint64_t to support bit operations for detecting when two faces share an edge. 
     * NB: Breaks if any cell has more than 64 vertices. 
     */
    vector<vector<uint64_t>> faces;

    /**
     * Generates an initial assignment that achieves spherical zones, 
     * which are a pre-condition for 
     * the correctness of the attempt_flip function. 
     * The zones are, in general, not population-balanced, 
     * since all but one zone is a singleton (only one cell). 
     * 
     * The procedure is deterministic, 
     * selecting for each new zone the cell that:
     *   1. is not a neighbor of other selected singleton zones; 
     *   2. has the most wall neighbors among cells satisfying condition 1; 
     *   3. satisfies condition (2) of attempt_flip w.r.t. the default (remainder) zone; and
     *   4. has the least index among cells satisfying conditions 1-3.  
     * 
     * TODO: Add theorem/proof to paper that this procedure produces spherical zones.
     */
    void generate_initial_assignment() {
        vector<bool> is_adjacent_to_singleton; // Indicators of whether each cell is adjacent to a singleton zone
        vector<int> num_wall_neighbors; // Number of neighbors of each cell that are walls
        vector<int> singleton_zone_cells; // Indices of cells used to create singleton zones

        is_adjacent_to_singleton.resize(N, false);
        num_wall_neighbors.resize(N, 0);
        for (int i = 0; i < N; ++i) {
            for (auto & i_neighbor : this->neighbors[i]) {
                if (i_neighbor < 0) ++num_wall_neighbors[i];
            }
        }
        
        while (singleton_zone_cells.size() < (size_t) K - 1) {
            int new_singleton_index = -1;
            for (int i = 0; i < N; ++i) {
                if (!is_adjacent_to_singleton[i]) {
                    if (new_singleton_index < 0 
                        || num_wall_neighbors[i] > num_wall_neighbors[new_singleton_index]) {
                        // Check condition (2) of attempt_flip w.r.t. default (remainder) zone. 
                        // All nonnegative neighbors are currently in the default zone. 
                        vector<int> old_part_face_ids;
                        vector<int> old_complement_face_ids;
                        for (int face_id = 0; face_id < surface_dual[i].size; ++face_id) {
                            string neighbor_name = surface_dual[i].vertex_name[face_id];
                            int neighbor_id = stoi(neighbor_name);
                            if (neighbor_id >= 0) {
                                old_part_face_ids.push_back(face_id);
                            } else {
                                old_complement_face_ids.push_back(face_id);
                            }
                        }
                        if (surface_dual[i].is_connected_subgraph(old_part_face_ids) 
                            && surface_dual[i].is_connected_subgraph(old_complement_face_ids)) {
                            new_singleton_index = i; // Pick this singleton for now. 
                                                     // Could replace if a later cell 
                                                     // has more wall neighbors. 
                        }
                    }
                }
            }

            if (new_singleton_index < 0) {
                return;
                // Alternative: throw "Failed to find a valid new singleton zone.";
            }
            singleton_zone_cells.push_back(new_singleton_index);
            is_adjacent_to_singleton[new_singleton_index] = true; // Consider singleton adjacent to itself to remove from candidates
            for (auto & neighbor_id : g.adjacency_list[new_singleton_index]) {
                is_adjacent_to_singleton[neighbor_id] = true; // Remove non-wall neighbors from candidates
            }
        }

        // Populate & set initial assignment
        vector<int> init_assignment;
        init_assignment.resize(N, K);
        int zone_index = 1;
        for (auto & singleton : singleton_zone_cells) {
            init_assignment[singleton] = zone_index++;
        }
        set_assignment(init_assignment);
    };

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
     * Checks whether flipping a cell to a new part maintains spherical zones (parts). 
     * If the flip is valid, it is made. 
     * 
     * \param[in] (cell_id) The ID of the cell to be flipped
     * \param[in] (new_part) The part to which cell_id will be moved. 
     *                       Enforced to be between 1 and K, inclusive; 
     *                       otherwise, throws an invalid_argument exception.
     */
    bool attempt_flip(int &cell_id, int &new_part);

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
