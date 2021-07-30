/** 
 * \file geograph3d.hh
 * \brief Header file for geograph3d. */

#ifndef GEOGRAPH3D_HH
#define GEOGRAPH3D_HH
#include "static_graph.hh"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
using std::string;
using std::vector;

namespace gg3d {
/** \brief Class representing an augmented neighbor. 
 * 
 * This is a class representing an augmented neighbor of 
 * a (3-D) cell in a 3-D geo-graph. 
 * It contains the augmented neighbor's ID and 
 * the shared faces of each dimension. 
 */
class augmented_neighbor {
private:
    /** Cell ID of the augmented neighbor */
    const size_t id;
    
public:
    /** IDs of the vertices shared with the augmented_neighbor */
    vector<size_t> shared_vertices;
    /** IDs of the edges shared with the augmented_neighbor */
    vector<size_t> shared_edges;
    /** IDs of the faces shared with the augmented_neighbor; should never have more than one */
    vector<size_t> shared_faces;

    /** Constructor for augmented_neighbor. 
     * 
     * \param[in] (id) The cell id of the augmented neighbor.
     */
    augmented_neighbor(const size_t id);
};

/** \brief Class representing a vertex (0-face) in 
 * the polyhedral cell complex. 
 * 
 * This is a class representing a vertex, i.e., 0-face, 
 * in the polyhedral cell complex from which 
 * a 3-D geo-graph is constructed. 
 * It contains an ID for the vertex and 
 * the IDs of the edges (1-faces) and faces (2-faces)
 * in which the vertex participates. 
 */
class cell_vertex {
private:
    /** Tolerance for comparing vertex coordinates. See operator == */
    const double tol = 1.0e-10;

public:
    /** Unique ID of the vertex (0-face) */
    const size_t id;
    /** x-y-z coordinates of the vertex */
    const double pos[3];
    /** IDs of incident edges (1-faces), in increasing ID order */
    vector<size_t> edges;
    /** IDs of incident faces (2-faces), in increasing ID order */
    vector<size_t> faces;

    /** Constructor for cell_vertex. 
     * 
     * \param[in] (id) The id of the cell vertex. 
     * \param[in] (x) The vertex's x-coordinate. 
     * \param[in] (y) The vertex's y-coordinate. 
     * \param[in] (z) The vertex's z-coordinate. 
     */
    cell_vertex(size_t id, double x, double y, double z);

    bool operator == (const cell_vertex& v) const {
        for (size_t i = 0; i < sizeof(pos) / sizeof(pos[0]); ++i) {
            double diff = pos[i] - v.pos[i];
            if (diff > tol || diff < - tol) {
                return false;
            }
        }
        return true;
    }
};

/** \brief Class representing an edge (1-face) in 
 * the polyhedral cell complex. 
 * 
 * This is a class representing an edge, i.e., 1-face, 
 * in the polyhedral cell complex from which 
 * a 3-D geo-graph is constructed. 
 * It contains an ID for the edge and 
 * the IDs of the two endpoint vertices (0-faces) and 
 * any faces (2-faces) in which the edge participates. 
 */
class cell_edge {
public:
    /** Unique ID of the edge (1-face) */
    const size_t id;
    /** IDs of endpoint vertices (0-faces), in increasing ID order; should always have two */
    vector<size_t> vertices;
    /** IDs of incident faces (2-faces), in increasing ID order; should always have two */
    vector<size_t> faces;

    /** Constructor for cell_edge. 
     * 
     * \param[in] (id) The id of the cell edge. 
     */
    cell_edge(size_t id);

    bool operator == (const cell_edge& e) const {
        if (vertices.size() != 2) return false;
        return (vertices[0] == e.vertices[0]) && (vertices[1] == e.vertices[1]);
    }
};

/** \brief Class representing a face (2-face) in 
 * the polyhedral cell complex. 
 * 
 * This is a class representing a face, i.e., 2-face, 
 * in the polyhedral cell complex from which 
 * a 3-D geo-graph is constructed. 
 * It contains an ID for the face and 
 * the IDs of the edges (1-faces) and 
 * vertices (0-faces) comprising the face. 
 */
class cell_face {
public:
    /** Unique ID of the face (2-face) */
    const size_t id;
    /** IDs of edges (1-faces) on the face, in increasing ID order */
    vector<size_t> edges;
    /** IDs of vertices (0-faces) on the face, in increasing ID order */
    vector<size_t> vertices;

    /** Constructor for cell_face. 
     * 
     * \param[in] (id) The id of the cell face. 
     */
    cell_face(size_t id);

    bool operator == (const cell_face& f) const {
        return vertices.size() == count_shared_vertices(f);
    }

    size_t count_shared_vertices(const cell_face& f) const {
        std::vector<size_t> intersection;
        std::set_intersection(vertices.begin(), vertices.end(), 
                              f.vertices.begin(), f.vertices.end(), 
                              std::back_inserter(intersection));
        return intersection.size();
    }
};

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
    size_t N;
    /** Number of parts */
    const size_t K;
    /** A map of cell IDs to part assignments, indexed 1 to K */
    vector<size_t> assignment;
    /** List of lists of neighbor IDs; primary index is cell ID. Includes wall neighbors (-1 through -6)*/
    vector<vector<size_t>> neighbors;
    /** List of lists of augmented neighbors; primary index is cell ID */
    vector<vector<augmented_neighbor>> aug_neighbors;
    /** List of cell vertices */
    vector<cell_vertex> cell_vertices;
    /** List of cell edges */
    vector<cell_edge> cell_edges;
    /** List of cell faces */
    vector<cell_face> cell_faces;

    /** 
     * List of lists of faces, stored as uint64_t to support bit operations for detecting when two faces share an edge. 
     * NB: Breaks if any cell has more than 64 vertices. 
     * TODO: Should replace this during rework
     */
    // vector<vector<uint64_t>> faces;

    /**
     * Generates an initial assignment that achieves spherical zones, 
     * which are a pre-condition for 
     * the correctness of the attempt_flip function. 
     * The zones are, in general, not population-balanced, 
     * since all but one zone is a singleton (only one cell). 
     * 
     * The procedure is deterministic, 
     * selecting for each new zone a cell that
     *   1. is not a neighbor of other selected singleton zones; 
     *   2. satisfies condition (1) of attempt_flip; and 
     *   3. satisfies condition (2) of attempt_flip w.r.t. the default (remainder) zone. 
     * 
     * Among cells with maximum number of wall neighbors satisfying the conditions, 
     * the cell with least index is chosen. 
     */
    // void generate_initial_assignment() {
    //     vector<bool> is_adjacent_to_singleton; // Indicators of whether each cell is adjacent to a singleton zone
    //     vector<int> num_wall_neighbors; // Number of neighbors of each cell that are walls
    //     vector<int> singleton_zone_cells; // Indices of cells used to create singleton zones

    //     is_adjacent_to_singleton.resize(N, false);
    //     num_wall_neighbors.resize(N, 0);
    //     for (int i = 0; i < N; ++i) {
    //         for (auto & i_neighbor : this->neighbors[i]) {
    //             if (i_neighbor < 0) ++num_wall_neighbors[i];
    //         }
    //     }
        
    //     while (singleton_zone_cells.size() < (size_t) K - 1) {
    //         int new_singleton_index = -1;
    //         for (int i = 0; i < N; ++i) {
    //             if (!is_adjacent_to_singleton[i]) {
    //                 if (new_singleton_index < 0 
    //                     || num_wall_neighbors[i] > num_wall_neighbors[new_singleton_index]) {
    //                     // Check condition (1) of attempt_flip
    //                     static_graph subgraph = aug_neighbor_graph[i];
    //                     vector<int> old_part_neighbor_new_ids; // New IDs (in aug neighborhood subgraph) of other cells from old part
    //                     for (int i = 0; i < subgraph.size; ++i) {
    //                         int vertex_id = stoi(subgraph.vertex_name[i]);
    //                         // Make sure the old part aug neighbor is not a singleton zone
    //                         vector<int>::iterator it = std::find(singleton_zone_cells.begin(), singleton_zone_cells.end(), vertex_id);
    //                         if (it == singleton_zone_cells.end()) {
    //                             old_part_neighbor_new_ids.push_back(i);
    //                         }
    //                     }
    //                     if (!subgraph.is_connected_subgraph(old_part_neighbor_new_ids)) {
    //                         continue;
    //                     }

    //                     // Check condition (2) of attempt_flip w.r.t. default (remainder) zone. 
    //                     // All nonnegative neighbors are currently in the default zone. 
    //                     vector<int> old_part_face_ids;
    //                     vector<int> old_complement_face_ids;
    //                     for (int face_id = 0; face_id < surface_dual[i].size; ++face_id) {
    //                         string neighbor_name = surface_dual[i].vertex_name[face_id];
    //                         int neighbor_id = stoi(neighbor_name);
    //                         if (neighbor_id >= 0) {
    //                             old_part_face_ids.push_back(face_id);
    //                         } else {
    //                             old_complement_face_ids.push_back(face_id);
    //                         }
    //                     }
    //                     if (surface_dual[i].is_connected_subgraph(old_part_face_ids) 
    //                         && surface_dual[i].is_connected_subgraph(old_complement_face_ids)) {
    //                         new_singleton_index = i; // Pick this singleton for now. 
    //                                                  // Could replace if a later cell 
    //                                                  // has more wall neighbors. 
    //                     }
    //                 }
    //             }
    //         }

    //         if (new_singleton_index < 0) {
    //             std::cout << "Failed to find a valid new singleton zone. Found " << singleton_zone_cells.size() << " so far:";
    //             for (auto & singleton : singleton_zone_cells) {
    //                 std::cout << "\t" << singleton << "\n";
    //             }
    //             std::cout << "\n";
    //             return;
    //             // Alternative: throw "Failed to find a valid new singleton zone.";
    //         }
    //         singleton_zone_cells.push_back(new_singleton_index);
    //         is_adjacent_to_singleton[new_singleton_index] = true; // Consider singleton adjacent to itself to remove from candidates
    //         for (auto & neighbor_id : g.adjacency_list[new_singleton_index]) {
    //             is_adjacent_to_singleton[neighbor_id] = true; // Remove non-wall neighbors from candidates
    //         }
    //     }

    //     // Populate & set initial assignment
    //     vector<int> init_assignment;
    //     init_assignment.resize(N, 1);
    //     int zone_index = 2;
    //     for (auto & singleton : singleton_zone_cells) {
    //         init_assignment[singleton] = zone_index++;
    //     }
    //     set_assignment(init_assignment);
    // };

public:
    /** 
     * The undirected graph of cells with 2-D face adjacencies, 
     * i.e., the 1-skeleton of the dual of the cell complex, 
     * omitting the vertex representing the infinite exterior. 
     */
    static_graph g;
    /** The cells' surface poset graphs, indexed by cell ID */
    vector<static_graph> surface_poset_graphs;

    /** Constructor for 3-D geo-graph. 
     * Builds representation of 3-D geo-graph from 
     * output of Voro++ representing a Voronoi tesselation 
     * inside of a 3-D box (rectangular prism). 
     * 
     * In particular, the constructor: 
     *  - Builds the cell adjacency graph. 
     *  - Builds surface poset graphs for each cell. 
     * 
     * \param[in] (in_filename) String name of the input CSV file, 
     *                          such as that produced by conversion_utils.
     *                          NB: This constructor assumes rows are sorted by cell ID, ascending. 
     * \param[in] (num_parts) The number of parts in the partition (stored as K). 
     */
    geograph3d(const string in_filename, const size_t num_parts);

    /**
     * Checks whether flipping a cell to a new part maintains spherical zones (parts). 
     * If the flip is valid, it is made. 
     * 
     * \param[in] (cell_id) The ID of the cell to be flipped
     * \param[in] (new_part) The part to which cell_id will be moved. 
     *                       Enforced to be between 1 and K, inclusive; 
     *                       otherwise, throws an invalid_argument exception.
     */
    bool attempt_flip(size_t &cell_id, size_t &new_part);

    /**
     * Setter for assignment member variable. 
     * 
     * \param[in] (assignment) New assignments, as a vector of part IDs (1 to K) 
     *                         indexed by cell ID (0 to N-1)
     */
    void set_assignment(vector<size_t> assignment) {
        this->assignment = assignment;
    }

    /** Getter for assignment member variable. */
    vector<size_t> get_assignment() { return assignment; }

    /** Returns the number of parts in the partition (i.e., member variable K). */
    size_t num_parts() { return K; }

    /** Returns the number of cells in the graph (i.e., member variable N). */
    size_t num_cells() { return N; }
};
}

#endif /* GEOGRAPH3D_HH */
