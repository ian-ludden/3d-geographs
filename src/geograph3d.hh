/** 
 * \file geograph3d.hh
 * \brief Header file for geograph3d. */
#ifndef GEOGRAPH3D_HH
#define GEOGRAPH3D_HH
#include "static_graph.hh"
#include <algorithm>
#include <iostream>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using std::string;
using std::vector;

namespace gg3d {
/** \brief Enum for result of geograph3d::attempt_flip. 
 * 
 * This enum comprises the possible statuses of a flip attempt.  
 * The status "fail_[i]" indicates that condition (i) has failed, 
 * for (i) in 1 through 5. 
 * The status "success" means the flip has succeeded and 
 * the geograph3d object has been updated accordingly. 
 */
enum flip_status { fail_1, fail_2, fail_3, fail_4, fail_5, success };

/** Convert flip_status enum value to string representation */
string flip_status_string(flip_status status);

/** \brief Class representing an augmented neighbor. 
 * 
 * This is a class representing an augmented neighbor of 
 * a (3-D) cell in a 3-D geo-graph. 
 * It contains the augmented neighbor's ID and 
 * the shared faces of each dimension. 
 */
class augmented_neighbor {
public:
    /** Cell ID of the augmented neighbor */
    const size_t id;
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
    /** Flag for whether the vertex lies on a zone boundary */
    bool is_boundary;
    /** Flag for whether the vertex lies on the outer boundary (i.e., a bounding wall) */
    bool is_outer_boundary;

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
    /** Setter for is_boundary private member variable */
    void set_is_boundary(bool new_is_boundary);
    /** Getter for is_boundary private member variable */
    bool get_is_boundary();
    /** Setter for is_outer_boundary private member variable */
    void set_is_outer_boundary(bool new_is_outer_boundary);
    /** Getter for is_outer_boundary private member variable */
    bool get_is_outer_boundary();
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
private:
    /** Flag for whether the edge lies on a zone boundary */
    bool is_boundary;
    /** Flag for whether the edge lies on the outer boundary (i.e., a bounding wall) */
    bool is_outer_boundary;
public:
    /** Unique ID of the edge (1-face) */
    const size_t id;
    /** IDs of endpoint vertices (0-faces), in increasing ID order; should always have two */
    vector<size_t> vertices;
    /** IDs of incident faces (2-faces), in increasing ID order; should always have at least two */
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

    /** Setter for is_boundary private member variable */
    void set_is_boundary(bool new_is_boundary);
    /** Getter for is_boundary private member variable */
    bool get_is_boundary();
    /** Setter for is_outer_boundary private member variable */
    void set_is_outer_boundary(bool new_is_outer_boundary);
    /** Getter for is_outer_boundary private member variable */
    bool get_is_outer_boundary();
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
private:
    /** Flag for whether the face lies on a zone boundary */
    bool is_boundary;
    /** Flag for whether the face lies on the outer boundary (i.e., a bounding wall) */
    bool is_outer_boundary;
public:
    /** Unique ID of the face (2-face) */
    const size_t id;
    /** IDs of edges (1-faces) on the face, in increasing ID order */
    vector<size_t> edges;
    /** IDs of vertices (0-faces) on the face, in increasing ID order */
    vector<size_t> vertices;
    /** ID of lesser-indexed cell */
    size_t first_cell;
    /** ID of greater-indexed cell, if any (0 for outer boundary faces) */
    size_t second_cell;

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

    /** Setter for is_boundary private member variable */
    void set_is_boundary(bool new_is_boundary);
    /** Getter for is_boundary private member variable */
    bool get_is_boundary();
    /** Setter for is_outer_boundary private member variable */
    void set_is_outer_boundary(bool new_is_outer_boundary);
    /** Getter for is_outer_boundary private member variable */
    bool get_is_outer_boundary();
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
    /** Part sizes, indexed 1 to K (part_sizes[0] is set to 0) */
    vector<size_t> part_sizes;
    /** A map of cell IDs to part assignments, indexed 1 to K */
    vector<size_t> assignment;
    /** List of lists of neighbor IDs; primary index is cell ID. Includes wall neighbors (-1 through -6)*/
    vector<vector<int>> neighbors;
    /** List of lists of augmented neighbors; primary index is cell ID */
    vector<vector<augmented_neighbor>> aug_neighbors;
    /** List of cell vertices */
    vector<cell_vertex> cell_vertices;
    /** List of cell edges */
    vector<cell_edge> cell_edges;
    /** List of cell faces */
    vector<cell_face> cell_faces;
    /** List of cell IDs */
    vector<size_t> cell_ids;

    /**
     * Generates an initial assignment that achieves spherical zones, 
     * which are a pre-condition for 
     * the correctness of the attempt_flip function. 
     * The zones are, in general, not population-balanced, 
     * since all but one zone is a singleton (only one cell). 
     * 
     * The procedure is deterministic, 
     * selecting for each new zone a cell that
     * is not a neighbor of other selected singleton zones. 
     * 
     * Among cells with maximum number of wall neighbors satisfying the conditions, 
     * the cell with least index is chosen. 
     * For a uniform grid of unit cubes, this prioritizes corners. 
     */
    void generate_initial_assignment() {
        vector<bool> is_adjacent_to_singleton; // Indicators of whether each cell is adjacent to a singleton zone, indexed by cell ID
        vector<int> num_wall_neighbors; // Number of neighbors of each cell that are walls, indexed by cell ID
        vector<int> singleton_zone_cells; // IDs of cells used to create singleton zones

        is_adjacent_to_singleton.resize(N, false);
        num_wall_neighbors.resize(N, 0);
        for (size_t i = 0; i < N; ++i) {
            size_t cell_id = cell_ids[i];
            for (auto &neighbor_id : this->neighbors[cell_id]) {
                if (neighbor_id < 0) ++num_wall_neighbors[cell_id];
            }
        }
        
        while (singleton_zone_cells.size() < (size_t) K - 1) {
            int new_singleton_index = -1;
            for (size_t i = 0; i < N; ++i) {
                size_t cell_id = cell_ids[i];
                if (!is_adjacent_to_singleton[cell_id]) {
                    if (new_singleton_index < 0 
                        || num_wall_neighbors[cell_id] > num_wall_neighbors[cell_ids[new_singleton_index]]) {
                        new_singleton_index = i;
                    }
                }
            }

            if (new_singleton_index < 0) {
                throw std::runtime_error("Failed to find a valid new singleton zone.");
            }
            singleton_zone_cells.push_back(new_singleton_index);
            is_adjacent_to_singleton[cell_ids[new_singleton_index]] = true; // Consider singleton adjacent to itself to remove from candidates
            for (auto & neighbor_id : g.adjacency_list[cell_ids[new_singleton_index]]) {
                is_adjacent_to_singleton[neighbor_id] = true; // Remove non-wall neighbors from candidates
            }
        }

        // Populate & set initial assignment
        vector<size_t> init_assignment;
        init_assignment.resize(N, 1);
        size_t zone_index = 2;
        for (auto & singleton_cell_index : singleton_zone_cells) {
            size_t singleton_id = cell_ids[singleton_cell_index];
            part_sizes[zone_index] = 1;
            init_assignment[singleton_id] = zone_index++;
        }
        part_sizes[1] = N - singleton_zone_cells.size();
        set_assignment(init_assignment);

        /** (inefficiently) update is_boundary flags of cell vertices, edges, and faces */
        update_boundary_flags();
    };

    /**
     * Computes the shared elements (vertices, edges, faces) of 
     * a given cell with augmented neighbors in a given part. 
     * 
     * \param[in] (cell_id) The ID of the cell
     * \param[in] (part) The ID of the part (1 to K) with which we want the shared elements
     * \param[in] (dimension) 0 for vertices, 1 for edges, 2 for faces
     * 
     * \return A std::set<string> of the descriptive names of the elements, i.e., 
     *         vXX for vertices, eXX for edges, fXX for faces, where XX is the element's ID
     */
    std::set<string> shared_elements_with_part(size_t cell_id, size_t part, uint8_t dimension) {
        std::set<string> shared_elements;
        for (auto &aug_neighbor : aug_neighbors[cell_id]) {
            if (assignment[aug_neighbor.id] != part) continue;
            vector<size_t> shared_elements_indices;
            if (dimension == 0) {
                shared_elements_indices = aug_neighbor.shared_vertices;
            } else if (dimension == 1) {
                shared_elements_indices = aug_neighbor.shared_edges;
            } else if (dimension == 2) {
                shared_elements_indices = aug_neighbor.shared_faces;
            } else {
                throw std::invalid_argument("Invalid argument: dimension must be 0, 1, or 2.");
            }
            for (auto &shared_element_index : shared_elements_indices) {
                string prefixes[3] = { "v", "e", "f" };
                string name_str = prefixes[dimension] + std::to_string(shared_element_index);
                shared_elements.insert(name_str);
            }
        }
        return shared_elements;
    }

    /**
     * Computes Y_v, the vertices and edges on the boundary of the shared surface between 
     * a cell and a part. 
     * 
     * The function starts by marking all nodes in the cell's surface poset graph, G_v, that are in S_1, 
     * the shared surface. It then performs a breadth-first search over these shared nodes, 
     * computing for each edge-node its degree in the induced subgraph G_v[S_1]. 
     * Edges with degree 3 (two vertices from the endpoints, and one of the two faces) 
     * must be on the boundary. 
     * 
     * Vertices on the boundary are then found by taking the union of the endpoints of the boundary edges. 
     */
    vector<string> boundary_vertices_and_edges_of_shared_surface(size_t &cell_id, vector<string> &shared_element_names) {
        if (shared_element_names.size() <= 0) throw std::invalid_argument("shared_element_names cannot be empty.\n");

        static_graph g_v = surface_poset_graphs[cell_id];
        vector<string> boundary_edges;
        vector<bool> is_found(shared_element_names.size(), false); // Indicator vector for visited nodes
        vector<bool> is_boundary(shared_element_names.size(), false); // Indicator vector for visited nodes
        vector<bool> is_available(shared_element_names.size(), false); // Indicator vector for nodes in the induced subgraph
        vector<size_t> node_id(shared_element_names.size());
        for (size_t i = 0; i < shared_element_names.size(); ++i) {
            node_id[i] = g_v.get_index_of_name(shared_element_names[i]);
            is_available[node_id[i]] = true;
        }

        // Run BFS starting from the first node
        std::queue<size_t> frontier;
        is_found[node_id[0]] = true;
        frontier.push(node_id[0]);

        while (!frontier.empty()) {
            size_t current_node = frontier.front();
            frontier.pop();

            string node_name = g_v.vertex_name[current_node];
            
            size_t count_neighbors_in_shared_surface = 0;

            // Add unvisited neighbors to the queue if in the induced subgraph
            for (auto &neighbor_node: g_v.adjacency_list[current_node]) {
                if (is_available[neighbor_node] && !is_found[neighbor_node]) {
                    is_found[neighbor_node] = true;
                    frontier.push(neighbor_node);
                }
                if (is_available[neighbor_node]) {
                    count_neighbors_in_shared_surface++;
                }
            }

            // Add degree-3 edge-nodes in G_v[S_1] to boundary_edges
            if (node_name[0] == 'e' && count_neighbors_in_shared_surface == 3) {
                boundary_edges.push_back(node_name);
            }
        }

        std::sort(boundary_edges.begin(), boundary_edges.end());

        // Compute boundary_vertices
        std::set<string> boundary_vertices_set;
        for (auto &e : boundary_edges) {
            size_t e_id = std::stoi(e.substr(1, e.size()));
            for (auto &vertex_neighbor : cell_edges[e_id].vertices) {
                string vertex_neighbor_name = "v" + std::to_string(vertex_neighbor);
                boundary_vertices_set.insert(vertex_neighbor_name);
            }
        }

        // Return union of boundary vertices and edges
        vector<string> boundary_vertices{boundary_vertices_set.begin(), boundary_vertices_set.end()};
        std::sort(boundary_vertices.begin(), boundary_vertices.end());
        vector<string> boundary_vertices_and_edges;
        for (auto &boundary_edge : boundary_edges) boundary_vertices_and_edges.push_back(boundary_edge);
        for (auto &boundary_vertex : boundary_vertices) boundary_vertices_and_edges.push_back(boundary_vertex);
        return boundary_vertices_and_edges;
    }

    /** 
     * Updates the is_boundary flags for cell faces after a flip succeeds. 
     * 
     * \param[in] (cell_id) the ID of the flipped cell
     * \param[in] (old_part) the ID of the part (1 to K) **from** which the cell was flipped
     * \param[in] (new_part) the ID of the part (1 to K) **to** which the cell was flipped
     */
    void update_boundary_faces(const size_t cell_id, const size_t old_part, const size_t new_part) {
        for (auto &aug_neighbor : aug_neighbors[cell_id]) {
            if (aug_neighbor.shared_faces.empty()) continue;

            size_t neighbor_id = aug_neighbor.id;
            size_t shared_face_index = aug_neighbor.shared_faces[0];
            
            if (assignment[neighbor_id] == old_part) {
                cell_faces[shared_face_index].set_is_boundary(true);
            } else if (assignment[neighbor_id] == new_part) {
                cell_faces[shared_face_index].set_is_boundary(false);
            } else {
                (void) 0; // Pass
            }
        }
    }

    /** 
     * Update the part sizes after a flip. 
     * 
     * \param[in] (old_part) the old part (1 to K) of the flipped unit
     * \param[in] (new_part) the new part (1 to K) of the flipped unit
     */
    void update_part_sizes(size_t old_part, size_t new_part) {
        part_sizes[old_part]--;
        part_sizes[new_part]++;
    }

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
    flip_status attempt_flip(size_t &cell_id, size_t &new_part);

    /**
     * Checks whether flipping a cell to a new part 
     * maintains contiguous, but not necessarily spherical, zones (parts). 
     * If the flip does maintain contiguous zones, 
     * this->assignment is updated to indicate the flip. 
     * 
     * \param[in] (cell_id) The ID of the cell to be flipped
     * \param[in] (new_part) The part to which cell_id will be moved. 
     *                       Enforced to be between 1 and K, inclusive; 
     *                       otherwise, throws an invalid_argument exception. 
     * 
     * \return true if the flip succeeds, false if it is rejected
     */
    bool attempt_flip_BFS(size_t &cell_id, size_t &new_part);

    /**
     * Checks whether flipping a cell to a new part 
     * maintains contiguous, but not necessarily spherical, zones (parts). 
     * Unlike `attempt_flip_BFS`, does not actually perform the flip 
     * if it would be successful. 
     * 
     * \param[in] (cell_id) The ID of the cell to be flipped
     * \param[in] (new_part) The part to which cell_id would be moved. 
     *                       Enforced to be between 1 and K, inclusive; 
     *                       otherwise, throws an invalid_argument exception. 
     * 
     * \return true if the flip would succeed, false if it is rejected
     */
    bool check_flip_BFS(size_t &cell_id, size_t &new_part);

    /**
     * Setter for assignment member variable. 
     * 
     * \param[in] (assignment) New assignments, as a vector of part IDs (1 to K) 
     *                         indexed by cell ID (0 to N-1)
     */
    void set_assignment(vector<size_t> assignment) {
        this->assignment = assignment;
        this->update_boundary_flags();
    }

    /** 
     * Updates is_boundary flags for all vertices and edges. 
     * Inefficient, since it examines every vertex and edge; 
     * typically used once after a call to set_assignment 
     * for a specific initial partition.
     */
    void update_boundary_flags();

    /**
     * Gets the size of the part with the given index (1 to K).
     * 
     * \param[in] (part_index) the index of the part (1 to K)
     * 
     * \return the number of cells in the part with the given index
     */
    size_t get_part_size(size_t part_index) {
        return this->part_sizes[part_index];
    }

    /** Getter for assignment member variable. */
    vector<size_t> get_assignment() { return assignment; }

    /** Getter for specific element of assignment member variable. */
    size_t get_assignment(size_t cell_id) {
        if (cell_id > N) throw std::invalid_argument("Out-of-bounds cell_id passed to get_assignment().");
        return this->assignment[cell_id];
    }

    /** Getter for cell_faces. */
    vector<cell_face> get_cell_faces() { return cell_faces; }

    /** Returns the number of parts in the partition (i.e., member variable K). */
    size_t num_parts() { return K; }

    /** Returns the number of cells in the graph (i.e., member variable N). */
    size_t num_cells() { return N; }

    /** Returns the size of the augmented neighborhood of a given cell. */
    size_t count_aug_neighbors(size_t cell_id) {
        return aug_neighbors[cell_id].size();
    }
};
}

#endif /* GEOGRAPH3D_HH */
