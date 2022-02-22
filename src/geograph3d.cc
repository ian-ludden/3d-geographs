/** 
 * \file geograph3d.cc
 * \brief Implementation of the 3-D geo-graph. */
#include "geograph3d.hh"
#include "static_graph.hh"
#include "../generate_input/conversion_utils.hh"
#include "../generate_input/vector_utils.hh"
#include <fstream>
#include <iostream>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
using std::string;
using std::vector;
\
namespace gg3d {
/** Convert flip_status enum value to string representation */
string flip_status_string(flip_status status) {
    switch (status)
    {
    case flip_status::fail_1:
        return "fail_1";
    case flip_status::fail_2:
        return "fail_2";
    case flip_status::fail_3:
        return "fail_3";
    case flip_status::fail_4:
        return "fail_4";
    case flip_status::fail_5:
        return "fail_5";
    case flip_status::success:
        return "success";
    default:
        return "unknown status";
    }
}

/** Constructor for augmented_neighbor. 
 * 
 * \param[in] (id) The cell id of the augmented neighbor.
 */
augmented_neighbor::augmented_neighbor(const size_t id) : id{id} {}

/** Constructor for cell_vertex. 
 * 
 * \param[in] (id) The id of the cell vertex. 
 */
cell_vertex::cell_vertex(size_t id, double x, double y, double z) : id{id}, pos{x, y, z} { 
    is_boundary = false; 
    is_outer_boundary = false; 
}

/** Setter for is_boundary private member variable.
 * Replaces this->is_boundary with new_is_boundary, 
 * unless this->is_outer_boundary is true, 
 * in which case this->is_boundary is left as true. 
 */
void cell_vertex::set_is_boundary(bool new_is_boundary) {
    this->is_boundary = new_is_boundary || this->is_outer_boundary;
}
/** Getter for is_boundary private member variable */
bool cell_vertex::get_is_boundary() {
    return this->is_boundary;
}
/** Setter for is_outer_boundary private member variable */
void cell_vertex::set_is_outer_boundary(bool new_is_outer_boundary) {
    this->is_outer_boundary = new_is_outer_boundary;
}
/** Getter for is_outer_boundary private member variable */
bool cell_vertex::get_is_outer_boundary() {
    return this->is_outer_boundary;
}

/** Constructor for cell_edge. 
 * 
 * \param[in] (id) The id of the cell edge. 
 */
cell_edge::cell_edge(size_t id) : id{id} { 
    is_boundary = false; 
    is_outer_boundary = false;
}

/** Setter for is_boundary private member variable.
 * Replaces this->is_boundary with new_is_boundary, 
 * unless this->is_outer_boundary is true, 
 * in which case this->is_boundary is left as true. 
 */
void cell_edge::set_is_boundary(bool new_is_boundary) {
    this->is_boundary = new_is_boundary || this->is_outer_boundary;
}
/** Getter for is_boundary private member variable */
bool cell_edge::get_is_boundary() {
    return this->is_boundary;
}
/** Setter for is_outer_boundary private member variable */
void cell_edge::set_is_outer_boundary(bool new_is_outer_boundary) {
    this->is_outer_boundary = new_is_outer_boundary;
}
    /** Getter for is_outer_boundary private member variable */
bool cell_edge::get_is_outer_boundary() {
    return this->is_outer_boundary;
}

/** Constructor for cell_face. 
 * 
 * \param[in] (id) The id of the cell face. 
 */
cell_face::cell_face(size_t id) : id{id} {
    first_cell = 0;
    second_cell = 0;
    is_boundary = false;
    is_outer_boundary = false;
}

/** Setter for is_boundary private member variable */
void cell_face::set_is_boundary(bool new_is_boundary) {
    this->is_boundary = new_is_boundary || this->is_outer_boundary;
}

/** Getter for is_boundary private member variable */
bool cell_face::get_is_boundary() {
    return this->is_boundary;
}
/** Setter for is_outer_boundary private member variable */
void cell_face::set_is_outer_boundary(bool new_is_outer_boundary) {
    this->is_outer_boundary = new_is_outer_boundary;
}
/** Getter for is_outer_boundary private member variable */
bool cell_face::get_is_outer_boundary() {
    return this->is_outer_boundary;
}

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
 *                          such as that produced by conversion_utils.
 *                          NB: This constructor assumes rows are sorted by cell ID, ascending. 
 * \param[in] (num_parts) The number of parts in the partition (stored as K). 
 * */
geograph3d::geograph3d(const string in_filename, const size_t num_parts) 
    : K{num_parts} {
    part_sizes.resize(K + 1, 0);
    csv_row row;
    std::ifstream in_file(in_filename);
    in_file >> row; // Read first row to get number of cells

    N = stoi(row[0]);
    
    // Edges of cell adjacency graph
    vector<Edge> g_edges;

    // List of vertex indices (by ID) for each cell, indexed by cell ID
    vector<vector<size_t>> vertex_indices;
    // List of face indices (by ID) for each cell, indexed by cell ID
    vector<vector<size_t>> face_indices;

    in_file >> row; // Skip row of column headers

    size_t cell_id;
    // size_t num_neighbors;

    neighbors.reserve(N);
    aug_neighbors.resize(N);

    // Read each row and store info for that cell
    while (in_file >> row) {
        vector<int> cell_neighbors_ints; // Temp list of neighbor IDs as type int
        vector<int> cell_neighbors; // Temp list of neighbor IDs as type size_t
        cell_id = stoi(row[0]);
        // num_neighbors = stoi(row[1]);
        cell_neighbors = delimited_list_to_vector_of_int(row[2], ' ');
        
        neighbors.push_back(cell_neighbors);

        // Add new edges (i.e., to neighbors with greater ID)
        for (int &neighbor_id : cell_neighbors) {
            // Omit negative IDs, which represent walls, from edges
            if ((int) cell_id < neighbor_id) g_edges.push_back({cell_id, (size_t) neighbor_id});
        }

        // Parse and index the (unique) vertices, checking for duplicates
        string vertices_string = (string) row[4];
        vector<size_t> vertex_ids_vec;

        size_t left_pos, right_pos;
        left_pos = vertices_string.find('(');
        while (left_pos < vertices_string.length()) {
            right_pos = vertices_string.find(')', left_pos);
            
            string cell_vertex_tuple = vertices_string.substr(left_pos, right_pos - left_pos + 1);
            vector<double> cell_vertex_coords = tuple_to_vector_of_double(cell_vertex_tuple);
            cell_vertex v(cell_vertices.size(), cell_vertex_coords[0], cell_vertex_coords[1], cell_vertex_coords[2]);
            bool is_new_vertex = true;

            // Check whether this vertex has already been indexed
            for (size_t i = 0; i < cell_vertices.size(); ++i) {
                if (cell_vertices[i] == v) {
                    is_new_vertex = false;
                    vertex_ids_vec.push_back(i);
                    break;
                }
            }
            
            if (is_new_vertex) {
                cell_vertices.push_back(v);
                vertex_ids_vec.push_back(v.id);
            }

            left_pos = vertices_string.find('(', right_pos);
        }

        vertex_indices.push_back(vertex_ids_vec);

        // Parse and index the (unique) faces, checking for duplicates
        string faces_string = string(row[3]);
        vector<size_t> face_ids_vec;

        left_pos = faces_string.find('(');
        // Iterate over face tuples, provided with commas as delimiters
        while (left_pos < faces_string.length()) {
            right_pos = faces_string.find(')', left_pos);
            
            string cell_face_tuple = faces_string.substr(left_pos + 1, right_pos - 1 - left_pos);
            vector<int> face_vertex_indices_ints = delimited_list_to_vector_of_int(cell_face_tuple, ',');
            vector<size_t> face_vertex_indices;
            for (auto &face_vertex_index_int : face_vertex_indices_ints) {
                face_vertex_indices.push_back((size_t) face_vertex_index_int);
            }
            vector<size_t> face_vertex_ids;
            for (auto &vertex_index : face_vertex_indices) {
                face_vertex_ids.push_back(vertex_ids_vec[vertex_index]);
            }
            cell_face f(cell_faces.size());
            std::sort(face_vertex_ids.begin(), face_vertex_ids.end());
            f.vertices = face_vertex_ids;

            bool is_new_face = true;

            // Check whether this face has already been indexed
            for (size_t i = 0; i < cell_faces.size(); ++i) {
                if (cell_faces[i] == f) {
                    is_new_face = false;
                    face_ids_vec.push_back(i);
                    break;
                }
            }

            if (is_new_face) {
                face_ids_vec.push_back(f.id);

                // Link vertices to new face
                for (auto &vertex_id : f.vertices) {
                    cell_vertices[vertex_id].faces.push_back(f.id);
                }

                vector<size_t> face_edge_ids;

                // Construct any new edges
                for (size_t v_index = 0; v_index < face_vertex_indices.size(); ++v_index) {
                    size_t first_v_id = vertex_ids_vec[face_vertex_indices[v_index]];
                    size_t second_v_id = (v_index == face_vertex_indices.size() - 1) ? 
                        vertex_ids_vec[face_vertex_indices[0]] : vertex_ids_vec[face_vertex_indices[v_index + 1]];

                    cell_edge e(cell_edges.size());
                    e.vertices.push_back(first_v_id);
                    e.vertices.push_back(second_v_id);
                    std::sort(e.vertices.begin(), e.vertices.end());

                    bool is_new_edge = true;
                    for (size_t edge_id = 0; edge_id < cell_edges.size(); ++edge_id) {
                        if (cell_edges[edge_id] == e) {
                            is_new_edge = false;
                            face_edge_ids.push_back(edge_id);
                            cell_edges[edge_id].faces.push_back(f.id);
                            break;
                        }
                    }

                    if (is_new_edge) {
                        e.faces.push_back(f.id);
                        cell_edges.push_back(e);
                        face_edge_ids.push_back(e.id);
                        

                        for (auto &vertex_id : e.vertices) {
                            cell_vertices[vertex_id].edges.push_back(e.id);
                        }
                    }
                }

                std::sort(face_edge_ids.begin(), face_edge_ids.end());
                f.edges = face_edge_ids;
                cell_faces.push_back(f);
            }

            left_pos = faces_string.find('(', right_pos);
        }

        face_indices.push_back(face_ids_vec);
    } /** end while(in_file >> row) to read each row of input file */

    // Construct cell adjacency graph
    g = static_graph(g_edges, N);

    // Find all augmented neighbor relationships
    for (size_t i = 0; i < N; ++i) {
        vector<size_t> face_indices_i = face_indices[i];
        std::sort(face_indices_i.begin(), face_indices_i.end());
        vector<size_t> vertex_indices_i = vertex_indices[i];
        std::sort(vertex_indices_i.begin(), vertex_indices_i.end());

        for (size_t j = i + 1; j < N; ++j) {
            vector<size_t> face_indices_j = face_indices[j];
            std::sort(face_indices_j.begin(), face_indices_j.end());
            vector<size_t> vertex_indices_j = vertex_indices[j];
            std::sort(vertex_indices_j.begin(), vertex_indices_j.end());

            // Build sorted vector of common elements
            vector<size_t> shared_vertices;
            std::set_intersection(vertex_indices_i.begin(), vertex_indices_i.end(), 
                                  vertex_indices_j.begin(), vertex_indices_j.end(), 
                                  std::back_inserter(shared_vertices));

            if (shared_vertices.size() > 0) {
                augmented_neighbor aug_neighbor_i(i);
                augmented_neighbor aug_neighbor_j(j);
                aug_neighbor_i.shared_vertices = shared_vertices;
                aug_neighbor_j.shared_vertices = shared_vertices;

                if (shared_vertices.size() == 2) {
                    cell_edge e(0);
                    e.vertices = shared_vertices; // already sorted
                    size_t edge_id; 
                    for (edge_id = 0; edge_id < cell_edges.size(); ++edge_id) {
                        if (cell_edges[edge_id] == e) break;
                    }
                    aug_neighbor_i.shared_edges.push_back(edge_id);
                    aug_neighbor_j.shared_edges.push_back(edge_id);
                } else if (shared_vertices.size() > 2) { // Must share a face
                    cell_face f(0);
                    f.vertices = shared_vertices; // already sorted
                    size_t face_id;
                    for (face_id = 0; face_id < cell_faces.size(); ++face_id) {
                        if (cell_faces[face_id] == f) break;
                    }
                    aug_neighbor_i.shared_faces.push_back(face_id);
                    aug_neighbor_j.shared_faces.push_back(face_id);
                    for (auto &edge_id : cell_faces[face_id].edges) {
                        aug_neighbor_i.shared_edges.push_back(edge_id);
                        aug_neighbor_j.shared_edges.push_back(edge_id);
                    }
                } else {} // Shares a single vertex, do nothing
                
                aug_neighbors[i].push_back(aug_neighbor_j);
                aug_neighbors[j].push_back(aug_neighbor_i);
            }
        }
    }

    // Mark outer boundary elements as boundary
    double min[3] = { __DBL_MAX__, __DBL_MAX__, __DBL_MAX__ };
    double max[3] = { -__DBL_MAX__, -__DBL_MAX__, -__DBL_MAX__ };
    for (auto &cell_vtx : cell_vertices) {
        for (size_t dimension = 0; dimension < 3; ++dimension) {
            if (cell_vtx.pos[dimension] < min[dimension]) min[dimension] = cell_vtx.pos[dimension];
            if (cell_vtx.pos[dimension] > max[dimension]) max[dimension] = cell_vtx.pos[dimension];
        } 
    }

    // Boundary edges and vertices
    for (size_t i = 0; i < cell_edges.size(); ++i) {
        cell_edge e = cell_edges[i];
        cell_vertex v0 = cell_vertices[e.vertices[0]];
        cell_vertex v1 = cell_vertices[e.vertices[1]];
        for (size_t dimension = 0; dimension < 3; ++dimension) {
            if (v0.pos[dimension] == v1.pos[dimension] 
                && (v0.pos[dimension] <= min[dimension] 
                    || v0.pos[dimension] >= max[dimension])) {
                cell_edges[i].set_is_boundary(true);
                cell_edges[i].set_is_outer_boundary(true);
                cell_vertices[v0.id].set_is_boundary(true);
                cell_vertices[v0.id].set_is_outer_boundary(true);
                cell_vertices[v1.id].set_is_boundary(true);
                cell_vertices[v1.id].set_is_outer_boundary(true);
            }
        }
    }

    // Outer boundary faces (all edges are outer boundary)
    for (size_t i = 0; i < cell_faces.size(); ++i) {
        bool is_outer_boundary_face = true;
        for (auto &e_index : cell_faces[i].edges) {
            is_outer_boundary_face = is_outer_boundary_face && cell_edges[e_index].get_is_outer_boundary();
        }
        if (is_outer_boundary_face) {
            cell_faces[i].set_is_outer_boundary(true);
            cell_faces[i].set_is_boundary(true);
        }
    }

    // Populate first_cell and second_cell for all cell faces
    for (size_t cell_id_1 = 0; cell_id_1 < this->N; ++cell_id_1) {
        for (auto &aug_neighbor : aug_neighbors[cell_id_1]) {
            if (aug_neighbor.shared_faces.empty()) continue;
            size_t shared_face_index = aug_neighbor.shared_faces[0];
            if (cell_id_1 < aug_neighbor.id) {
                cell_faces[shared_face_index].first_cell = cell_id_1;
                cell_faces[shared_face_index].second_cell = aug_neighbor.id;
            } else {
                cell_faces[shared_face_index].first_cell = aug_neighbor.id;
                cell_faces[shared_face_index].second_cell = cell_id_1;
            }
        }
    }

    /** Create surface poset graph for each cell */
    for (size_t i = 0; i < N; ++i) {
        vector<string> node_names; 
        vector<Edge> sp_edges;
        size_t num_faces = face_indices[i].size();
        size_t num_vertices = vertex_indices[i].size();
        size_t num_edges = num_vertices + num_faces - 2; // by Euler's formula

        // Add faces, prefixing names with "f"
        for (size_t f_index = 0; f_index < face_indices[i].size(); ++f_index) {
            cell_face f = cell_faces[face_indices[i][f_index]];
            string f_prefix = "f";
            size_t local_f_id = node_names.size();
            node_names.push_back(f_prefix + std::to_string(f.id));
            // Add edges, prefixing names with "e"
            string e_prefix = "e";
            for (auto &edge_index : f.edges) {
                bool is_new_edge;
                cell_edge e = cell_edges[edge_index];
                string e_name = e_prefix + std::to_string(e.id);
                // Check whether the edge has already been added to this surface poset graph
                size_t local_e_id = node_names.size();
                std::vector<string>::iterator it = std::find(node_names.begin(), node_names.end(), e_name);
                is_new_edge = it == node_names.end();

                if (is_new_edge) {
                    node_names.push_back(e_name);
                    sp_edges.push_back({local_f_id, local_e_id});
                } else {
                    local_e_id = std::distance(node_names.begin(), it);
                    sp_edges.push_back({local_e_id, local_f_id});
                    continue;
                }

                // Add vertices, prefixing names with "v"
                string v_prefix = "v";
                for (auto &vertex_index : e.vertices) {
                    bool is_new_vertex;
                    cell_vertex v = cell_vertices[vertex_index];
                    string v_name = v_prefix + std::to_string(v.id);
                    // Check whether the vertex has already been added to this surface poset graph
                    size_t local_v_id = node_names.size();
                    std::vector<string>::iterator it = std::find(node_names.begin(), node_names.end(), v_name);
                    is_new_vertex = it == node_names.end();
                    
                    if (is_new_vertex) {
                        node_names.push_back(v_name);
                        sp_edges.push_back({local_e_id, local_v_id});
                    } else {
                        local_v_id = std::distance(node_names.begin(), it);
                        sp_edges.push_back({local_v_id, local_e_id});
                    }
                }
            }
        }

        // // Debugging degeneracy issues with random Voronoi 20x20x20 instance
        // if (num_faces + num_edges + num_vertices > node_names.size()) {
        //     std::cout << "Node names:\n";
        //     for (auto & node_name : node_names) std::cout << node_name << "\n";
        //     std::cout << "\n";
        //     std::cout << "All faces:\n";
        //     for (auto & face_index : face_indices[i]) std::cout << "f" << face_index << "\n";
        //     std::cout << "\n";
        //     std::cout << "All vertices:\n";
        //     for (auto & vertex_index : vertex_indices[i]) std::cout << "v" << vertex_index << "\n";
        //     std::cout << "\n";
        //     std::cout << "All edges of cell's faces:\n";
        //     for (auto & face_index : face_indices[i]) {
        //         cell_face face = cell_faces[face_index];
        //         for (auto & edge_index : face.edges) std::cout << "e" << edge_index << "\n";
        //     }
        //     std::cout << "\n";
        // }

        static_graph surface_poset_graph = static_graph(sp_edges, num_faces + num_edges + num_vertices, node_names);
        surface_poset_graphs.push_back(surface_poset_graph);
    }

    // Now that all member variables are initialized, create initial assignment
    this->generate_initial_assignment();
}

/**
 * Checks whether flipping a cell to a new part maintains spherical zones (parts). 
 * If the flip is valid, it is made. 
 * 
 * \param[in] (cell_id) The ID of the cell to be flipped
 * \param[in] (new_part) The part to which cell_id will be moved. 
 *                       Enforced to be between 1 and K, inclusive; 
 *                       otherwise, throws an invalid_argument exception.
 */
flip_status geograph3d::attempt_flip(size_t &cell_id, size_t &new_part) {
    // Validate new_part (must be between 1 and K, inclusive)
    if (new_part < 1 || new_part > K) {
        throw std::invalid_argument("New part must be between 1 and K, inclusive.");
    }

    /** Construct S_1, the set of shared faces, edges, and vertices 
     * between the cell to flip and its current (i.e., giving) zone */
    vector<string> S_1;
    std::set<string> S_1_vertices = shared_elements_with_part(cell_id, assignment[cell_id], 0);
    std::set<string> S_1_edges = shared_elements_with_part(cell_id, assignment[cell_id], 1);
    std::set<string> S_1_faces = shared_elements_with_part(cell_id, assignment[cell_id], 2);
    for (auto &S_1_vertex : S_1_vertices) S_1.push_back(S_1_vertex);
    for (auto &S_1_edge : S_1_edges) S_1.push_back(S_1_edge);
    for (auto &S_1_face : S_1_faces) S_1.push_back(S_1_face);
    std::sort(S_1.begin(), S_1.end());

    vector<string> S_1_vertices_sorted{S_1_vertices.begin(), S_1_vertices.end()};
    std::sort(S_1_vertices_sorted.begin(), S_1_vertices_sorted.end());
    vector<string> S_1_edges_sorted{S_1_edges.begin(), S_1_edges.end()};
    std::sort(S_1_edges_sorted.begin(), S_1_edges_sorted.end());
    vector<string> S_1_faces_sorted{S_1_faces.begin(), S_1_faces.end()};
    std::sort(S_1_faces_sorted.begin(), S_1_faces_sorted.end());

    /** Construct X_v, the set of vertices/edges on the surface of cell v 
     * and the surface of v's current (giving) zone */
    vector<string> X_v;
    for (auto &node_name : surface_poset_graphs[cell_id].vertex_name) {
        size_t index = std::stoul(node_name.substr(1, node_name.size()));
        if ((node_name[0] == 'v' && cell_vertices[index].get_is_boundary()) 
            || (node_name[0] == 'e' && cell_edges[index].get_is_boundary())) {
            X_v.push_back(node_name);
        }
    }
    std::sort(X_v.begin(), X_v.end());

    /** Construct Y_v, the set of nodes representing edges and vertices 
     * on the boundary of the shared surface. This function returns a sorted vector. */
    vector<string> Y_v;
    try {
        Y_v = boundary_vertices_and_edges_of_shared_surface(cell_id, S_1);
    } catch (std::invalid_argument &invalid_arg_exception) {
        std::cerr << invalid_arg_exception.what();
        return flip_status::fail_1;
    }

    /** Check all five conditions, 
     * returning flip_status::fail_[i] 
     * if condition (i) fails. */

    /** Check condition (1): S_1 ∩ X_v = Y_v */
    vector<string> S_1_intersect_X_v; 
    std::set_intersection(S_1.begin(), S_1.end(), 
                          X_v.begin(), X_v.end(), 
                          std::back_inserter(S_1_intersect_X_v));


    // Both vectors are sorted, so == operator can be used to check equality
    if (S_1_intersect_X_v != Y_v) {
        return flip_status::fail_1;
    }

    /** Check condition (2):  G_v[S_1] is connected */
    bool satisfies_condition_2 = surface_poset_graphs[cell_id].is_connected_subgraph(S_1);
    if (!satisfies_condition_2) return flip_status::fail_2;

    /** Check condition (3): G_v[S_1^1 ∪ S_1^2] and 
     * G_v[V_v^1 \ S_1^1 ∪ V_v^2 \ S_1^2] are connected */
    vector<string> S_1_edges_and_faces;
    std::set_union(S_1_edges_sorted.begin(), S_1_edges_sorted.end(), 
                   S_1_faces_sorted.begin(), S_1_faces_sorted.end(),
                   std::back_inserter(S_1_edges_and_faces));

    vector<string> all_edges;
    vector<string> all_faces;
    for (auto &node_name : surface_poset_graphs[cell_id].vertex_name) {
        if (node_name[0] == 'e') {
            all_edges.push_back(node_name);
        } else if (node_name[0] == 'f') {
            all_faces.push_back(node_name);
        } else { (void) 0; } // Pass, i.e., ignore vertices
    }
    std::sort(all_edges.begin(), all_edges.end());
    std::sort(all_faces.begin(), all_faces.end());

    vector<string> S_1_edges_complement;
    std::set_difference(all_edges.begin(), all_edges.end(),
                        S_1_edges_sorted.begin(), S_1_edges_sorted.end(),
                        std::back_inserter(S_1_edges_complement));

    vector<string> S_1_faces_complement;
    std::set_difference(all_faces.begin(), all_faces.end(),
                        S_1_faces_sorted.begin(), S_1_faces_sorted.end(),
                        std::back_inserter(S_1_faces_complement));

    vector<string> S_1_edges_and_faces_complement;
    std::set_union(S_1_edges_complement.begin(), S_1_edges_complement.end(), 
                   S_1_faces_complement.begin(), S_1_faces_complement.end(),
                   std::back_inserter(S_1_edges_and_faces_complement));

    bool satisfies_condition_3 = surface_poset_graphs[cell_id].is_connected_subgraph(S_1_edges_and_faces)
                            &&   surface_poset_graphs[cell_id].is_connected_subgraph(S_1_edges_and_faces_complement);
    if (!satisfies_condition_3) return flip_status::fail_3;
    
    /** Construct S_2, the set of shared faces, edges, and vertices 
     * between the cell to flip and its new (i.e., receiving) zone */
    vector<string> S_2;
    std::set<string> S_2_vertices = shared_elements_with_part(cell_id, new_part, 0);
    std::set<string> S_2_edges = shared_elements_with_part(cell_id, new_part, 1);
    std::set<string> S_2_faces = shared_elements_with_part(cell_id, new_part, 2);
    for (auto &S_2_vertex : S_2_vertices) S_2.push_back(S_2_vertex);
    for (auto &S_2_edge : S_2_edges) S_2.push_back(S_2_edge);
    for (auto &S_2_face : S_2_faces) S_2.push_back(S_2_face);

    vector<string> S_2_vertices_sorted{S_2_vertices.begin(), S_2_vertices.end()};
    std::sort(S_2_vertices_sorted.begin(), S_2_vertices_sorted.end());
    vector<string> S_2_edges_sorted{S_2_edges.begin(), S_2_edges.end()};
    std::sort(S_2_edges_sorted.begin(), S_2_edges_sorted.end());
    vector<string> S_2_faces_sorted{S_2_faces.begin(), S_2_faces.end()};
    std::sort(S_2_faces_sorted.begin(), S_2_faces_sorted.end());

    /** Check condition (4): G_v[S_2] is connected */
    bool satisfies_condition_4 = surface_poset_graphs[cell_id].is_connected_subgraph(S_2);
    if (!satisfies_condition_4) return flip_status::fail_4;

    /** Check condition (5): G_v[S_2^1 ∪ S_2^2] and 
     * G_v[V_v^1 \ S_2^1 ∪ V_v^2 \ S_2^2] are connected */
    vector<string> S_2_edges_and_faces;
    std::set_union(S_2_edges_sorted.begin(), S_2_edges_sorted.end(), 
                   S_2_faces_sorted.begin(), S_2_faces_sorted.end(),
                   std::back_inserter(S_2_edges_and_faces));

    vector<string> S_2_edges_complement;
    std::set_difference(all_edges.begin(), all_edges.end(),
                        S_2_edges_sorted.begin(), S_2_edges_sorted.end(),
                        std::back_inserter(S_2_edges_complement));

    vector<string> S_2_faces_complement;
    std::set_difference(all_faces.begin(), all_faces.end(),
                        S_2_faces_sorted.begin(), S_2_faces_sorted.end(),
                        std::back_inserter(S_2_faces_complement));

    vector<string> S_2_edges_and_faces_complement;
    std::set_union(S_2_edges_complement.begin(), S_2_edges_complement.end(), 
                   S_2_faces_complement.begin(), S_2_faces_complement.end(),
                   std::back_inserter(S_2_edges_and_faces_complement));
    bool satisfies_condition_5 = surface_poset_graphs[cell_id].is_connected_subgraph(S_2_edges_and_faces)
                            &&   surface_poset_graphs[cell_id].is_connected_subgraph(S_2_edges_and_faces_complement);
    if (!satisfies_condition_5) return flip_status::fail_5;

    // Flip succeeded, so update assignment and part_sizes
    size_t old_part = assignment[cell_id];
    assignment[cell_id] = new_part;
    this->update_part_sizes(old_part, new_part);

    /** Update is_boundary flags of involved 
     * vertices, edges, and faces. 
     *  - Y_v remains external/boundary. 
     *  - S_1^0 ∪ S_1^1 \ Y_v becomes external/boundary. 
     *  - S_2^0 ∪ S_2^1 \ Y_v becomes internal/not boundary. 
     *  - face(s) between cell cell_id and other cells in new_part become internal/not boundary.
     *  - face(s) between cell cell_id and other cells in old_part become part of zone boundary. 
     */
    // Handle faces first
    this->update_boundary_faces(cell_id, old_part, new_part);

    // Handle edges and vertices together
    vector<string> new_boundary_elements;
    vector<string> S_1_vertices_and_edges;
    std::set_union(S_1_vertices_sorted.begin(), S_1_vertices_sorted.end(), 
                   S_1_edges_sorted.begin(), S_1_edges_sorted.end(),
                   std::back_inserter(S_1_vertices_and_edges));
    std::set_difference(S_1_vertices_and_edges.begin(), S_1_vertices_and_edges.end(), 
                        Y_v.begin(), Y_v.end(), 
                        std::back_inserter(new_boundary_elements));

    for (auto &element_name : new_boundary_elements) {
        size_t index = std::stoul(element_name.substr(1, element_name.size()));
        if (element_name[0] == 'v') {
            cell_vertices[index].set_is_boundary(true);
        } else { // Must be "eXX", by construction
            cell_edges[index].set_is_boundary(true);
        }
    }

    /** Construct the set of nodes representing edges and vertices 
     * on the boundary of the shared surface with the receiving zone (boundary of S_2). 
     * This function returns a sorted vector. */
    vector<string> S_2_boundary;
    try {
        S_2_boundary = boundary_vertices_and_edges_of_shared_surface(cell_id, S_2);
    } catch (std::invalid_argument &invalid_arg_exception) {
        std::cerr << invalid_arg_exception.what();
        return flip_status::fail_5; // Should probably introduce a new failure status for successful flips that fail to update boundaries
    }

    vector<string> new_internal_elements;
    vector<string> S_2_vertices_and_edges;
    std::set_union(S_2_vertices_sorted.begin(), S_2_vertices_sorted.end(),
                   S_2_edges_sorted.begin(), S_2_edges_sorted.end(), 
                   std::back_inserter(S_2_vertices_and_edges));
    std::set_difference(S_2_vertices_and_edges.begin(), S_2_vertices_and_edges.end(), 
                        S_2_boundary.begin(), S_2_boundary.end(), 
                        std::back_inserter(new_internal_elements));

    for (auto &element_name : new_internal_elements) {
        size_t index = std::stoul(element_name.substr(1, element_name.size()));
        if (element_name[0] == 'v') {
            cell_vertices[index].set_is_boundary(false);
        } else { // Must be "eXX", by construction
            cell_edges[index].set_is_boundary(false);
        }
    }

    return flip_status::success;
} /** end attempt_flip */

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
bool geograph3d::attempt_flip_BFS(size_t &cell_id, size_t &new_part) { 
    bool flip_will_succeed = check_flip_BFS(cell_id, new_part);
    
    if (flip_will_succeed) {
        size_t old_part = assignment[cell_id];
        assignment[cell_id] = new_part;
        this->update_boundary_faces(cell_id, old_part, new_part);
        this->update_part_sizes(old_part, new_part);
    }
    
    return flip_will_succeed;
}

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
bool geograph3d::check_flip_BFS(size_t &cell_id, size_t &new_part) {
    // Validate new_part (must be between 1 and K, inclusive)
    if (new_part < 1 || new_part > K) {
        throw std::invalid_argument("New part must be between 1 and K, inclusive.");
    }

    // Make sure cell_id has at least one neighbor in new_part
    bool has_neighbor_in_new_part = false;
    for (auto &neighbor_id : g.adjacency_list[cell_id]) {
        has_neighbor_in_new_part = has_neighbor_in_new_part || (assignment[neighbor_id] == new_part);
    }
    if (!has_neighbor_in_new_part) return false;

    // Create set of unvisited zone neighbors, that is, neighbors of cell cell_id in same zone currently
    std::set<size_t> unvisited_neighbors;
    for (auto &neighbor_id : g.adjacency_list[cell_id]) {
        if (assignment[neighbor_id] == assignment[cell_id]) unvisited_neighbors.insert(neighbor_id);
    }

    // Make sure unvisited_neighbors is not empty (if it is, then flip fails)
    if (unvisited_neighbors.empty()) return false;

    // Indicator vector for whether every cell is in the same zone as cell_id, to be used for BFS
    vector<bool> is_same_zone(this->N, false);
    for (size_t id = 0; id < this->N; ++id) is_same_zone[id] = (assignment[id] == assignment[cell_id]);
    is_same_zone[cell_id] = false;

    // Run BFS from arbitrary zone neighbor to make sure all zone neighbors are reached
    auto first_neighbor_iterator = unvisited_neighbors.begin();
    size_t start = *first_neighbor_iterator;
    unvisited_neighbors.erase(first_neighbor_iterator);

    // Indicator vector for visited cells, indexed by cell id
    vector<bool> is_visited(this->N, false);
    is_visited[start] = true;
    
    // Queue for managing BFS
    std::queue<size_t> frontier;
    frontier.push(start);

    while (!frontier.empty()) {
        size_t current_vertex = frontier.front();
        frontier.pop();

        // Add unvisited neighbors to frontier if in same zone
        for (auto &neighbor_id : g.adjacency_list[current_vertex]) {
            if (is_same_zone[neighbor_id] && !is_visited[neighbor_id]) {
                is_visited[neighbor_id] = true;
                frontier.push(neighbor_id);

                // Remove neighbor_id from unvisited_neighbors, if present
                auto it = unvisited_neighbors.find(neighbor_id);
                if (it != unvisited_neighbors.end()) {
                    unvisited_neighbors.erase(it);
                }
            }

            // Stop early if all zone neighbors have been visited
            if (unvisited_neighbors.empty()) break;
        }
    }

    return unvisited_neighbors.empty(); 
}

/** 
 * Updates is_boundary flags for all vertices and edges. 
 * Inefficient, since it examines every augmented neighbor relationship; 
 * typically used once after a call to set_assignment 
 * for a specific initial partition.
 */
void geograph3d::update_boundary_flags() {
    // Reset all to false
    for (size_t e_index = 0; e_index < cell_edges.size(); ++e_index) cell_edges[e_index].set_is_boundary(false);
    for (size_t v_index = 0; v_index < cell_vertices.size(); ++v_index) cell_vertices[v_index].set_is_boundary(false);
    for (size_t f_index = 0; f_index < cell_faces.size(); ++f_index) cell_faces[f_index].set_is_boundary(false);

    for (size_t id = 0; id < this->N; ++id) {
        for (size_t j = 0; j < aug_neighbors[id].size(); ++j) {
            augmented_neighbor other = aug_neighbors[id][j];
            if (assignment[other.id] != assignment[id]) {
                // Shared vertices, edges, and faces between zones should have is_boundary == true
                for (auto &e_index : other.shared_edges) cell_edges[e_index].set_is_boundary(true);
                for (auto &v_index : other.shared_vertices) cell_vertices[v_index].set_is_boundary(true);
                for (auto &f_index : other.shared_faces) cell_faces[f_index].set_is_boundary(true);
            }
        }
    }
}

} /** end namespace gg3d */

/** Uncomment below to run interactive program allowing user to give flips to attempt */
// int main(int argc, char *argv[]) {
//     string in_filename;
//     int K;
    
//     if (argc >= 3) {
//         in_filename = argv[1];
//         K = std::stoi(argv[2]);
//     } else {
//         string prog_name = argv[0];
//         if (prog_name.empty()) prog_name = "<program name>";
//         cout << "Usage: " << prog_name << " <input_filename> <number_of_parts>\n";
//         return 0;
//     }

//     gg3d::geograph3d geograph = gg3d::geograph3d(in_filename, K);

//     // Summarize geo-graph
//     cout << "The given 3-D geo-graph contains " << geograph.num_cells() << " cells ";
//     cout << "partitioned into " << geograph.num_parts() << " parts.\n";

//     // Summarize cell adjacency graph
//     size_t count_g_edges = 0;
//     for (auto & adj_list : geograph.g.adjacency_list) count_g_edges += adj_list.size();
//     count_g_edges = count_g_edges / 2; // Every edge is double-counted when summing adjacency list sizes
//     cout << "The cell adjacency graph has " << geograph.g.size << " vertices and " << count_g_edges << " edges.\n\n";

//     // Print cell part assignments
//     cout << "Current assignments:\n";
//     vector<size_t> current_assignment = geograph.get_assignment();
//     for (size_t i = 0; i < current_assignment.size(); ++i) {
//         cout << i << " to part " << current_assignment[i] << "\n";
//     }

//     // Print current boundary faces
//     cout << "Zone boundary faces:\n";
//     vector<gg3d::cell_face> cell_faces = geograph.get_cell_faces();
//     for (size_t i = 0; i < cell_faces.size(); ++i) {
//         if (cell_faces[i].get_is_boundary() && !cell_faces[i].get_is_outer_boundary()) {
//             cout << i << ", between cells " << cell_faces[i].first_cell << " and " << cell_faces[i].second_cell << ".\n";
//         }
//     }

//     // Attempt flips
//     size_t cell_id;
//     string response = "Y";
//     while (response[0] == 'Y' || response[0] == 'y') {
//         cout << "Enter the name of the unit/cell to be flipped: ";
//         string name;
//         std::cin >> name;
//         cell_id = stoi(name);
//         if (cell_id >= geograph.g.size) {
//             throw std::runtime_error("Provided cell id is out of bounds.");
//         }
//         current_assignment = geograph.get_assignment();
//         size_t part = current_assignment[cell_id];
//         cout << "\nUnit " << name << " is currently assigned to part " << part << ".\n";
//         cout << "Enter the name of its new part: ";
//         std::cin >> name;
//         size_t new_part = stoi(name);
//         // bool success = geograph.attempt_flip_BFS(cell_id, new_part);
//         gg3d::flip_status result = geograph.attempt_flip(cell_id, new_part);
//         // if (success) { // for attempt_flip_BFS
//         if (result == gg3d::flip_status::success) {
//             cout << "Flip was successful. Unit " << cell_id << " is now assigned to part " << new_part << ".\n";
//         } else {
//             // cout << "Flip failed.\n"; // for attempt_flip_BFS
//             cout << "Flip failed with status " << gg3d::flip_status_string(result) << ".\n";
//         }
        
//         cout << "Try another flip? (Y/N) ";
//         std::cin >> response;
//     }
// } /** end main */
