/** 
 * \file geograph3d.cc
 * \brief Implementation of the 3-D geo-graph. */
#include "geograph3d.hh"
#include "static_graph.hh"
#include "../generate_input/conversion_utils.hh"
#include "../generate_input/vector_utils.hh"
#include <bitset>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
using std::string;
using std::vector;

namespace gg3d {
/** Constructor for augmented_neighbor. 
 * 
 * \param[in] (id) The cell id of the augmented neighbor.
 */
augmented_neighbor::augmented_neighbor(const size_t id) : id{id} {}

/** Constructor for cell_vertex. 
 * 
 * \param[in] (id) The id of the cell vertex. 
 */
cell_vertex::cell_vertex(size_t id, double x, double y, double z) : id{id}, pos{x, y, z} {}

/** Setter for is_boundary private member variable */
void cell_vertex::set_is_boundary(bool new_is_boundary) {
    this->is_boundary = new_is_boundary;
}
/** Getter for is_boundary private member variable */
bool cell_vertex::get_is_boundary() {
    return this->is_boundary;
}

/** Constructor for cell_edge. 
 * 
 * \param[in] (id) The id of the cell edge. 
 */
cell_edge::cell_edge(size_t id) : id{id} {}

/** Setter for is_boundary private member variable */
void cell_edge::set_is_boundary(bool new_is_boundary) {
    this->is_boundary = new_is_boundary;
}
/** Getter for is_boundary private member variable */
bool cell_edge::get_is_boundary() {
    return this->is_boundary;
}

/** Constructor for cell_face. 
 * 
 * \param[in] (id) The id of the cell face. 
 */
cell_face::cell_face(size_t id) : id{id} {}

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
            // Omit negative IDs, which represent walls
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
    } /** end while(in_file >> row) */

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

    // // Check summary measurements of aug neighborhoods and vertices/edges/faces
    // int degree_freq[30];

    // // Count degrees of units with respect to augmented neighborhood
    // for (size_t i = 0; i < 30; ++i) degree_freq[i] = 0;

    // for (size_t i = 0; i < N; ++i) {
    //     size_t degree = aug_neighbors[i].size();
    //     degree_freq[degree] = degree_freq[degree] + 1;
    // }
    // cout << "Frequencies of augmented neighborhood sizes:\n";
    // for (size_t i = 0; i < 30; ++i) cout << "\t" << i << ": " << degree_freq[i] << " units.\n";


    // // Count degrees of vertices
    // for (size_t i = 0; i < 10; ++i) degree_freq[i] = 0;

    // for (size_t i = 0; i < cell_vertices.size(); ++i) {
    //     cell_vertex current_vtx = cell_vertices[i];
    //     int degree = current_vtx.edges.size();
    //     degree_freq[degree] = degree_freq[degree] + 1;
    // }

    // // Count # edges per face
    // for (size_t i = 0; i < 10; ++i) degree_freq[i] = 0;
    
    // for (size_t i = 0; i < cell_faces.size(); ++i) {
    //     cell_face current_face = cell_faces[i];
    //     int degree = current_face.edges.size();
    //     degree_freq[degree] = degree_freq[degree] + 1;
    // }

    // // Count # vertices per face
    // for (size_t i = 0; i < 10; ++i) degree_freq[i] = 0;
    
    // for (size_t i = 0; i < cell_faces.size(); ++i) {
    //     cell_face current_face = cell_faces[i];
    //     int degree = current_face.vertices.size();
    //     degree_freq[degree] = degree_freq[degree] + 1;
    // }

    // // Count # faces per edge
    // for (size_t i = 0; i < 10; ++i) degree_freq[i] = 0;
    
    // for (size_t i = 0; i < cell_edges.size(); ++i) {
    //     cell_edge current_edge = cell_edges[i];
    //     int degree = current_edge.faces.size();
    //     degree_freq[degree] = degree_freq[degree] + 1;
    // }

    // cout << "end of constructor.\n";

    // // Construct surface dual graphs
    // for (int cell_id = 0; cell_id < N; ++cell_id) {
    //     vector<int> cell_neighbors = neighbors[cell_id];
    //     vector<uint64_t> cell_faces = faces[cell_id];

    //     vector<Edge> surface_dual_edges;
    //     vector<string> surface_dual_vertex_names;
        
    //     for (size_t i = 0; i < cell_faces.size(); ++i) {
    //         surface_dual_vertex_names.push_back(std::to_string(cell_neighbors[i]));

    //         for (size_t j = i + 1; j < cell_faces.size(); ++j) {
    //             size_t count_shared_vertices = std::bitset<64>{cell_faces[i] & cell_faces[j]}.count();
    //             // Two faces on the surface share an edge if they share (exactly) two vertices. 
    //             if (count_shared_vertices == 2) {
    //                 Edge new_edge = {(int) i, (int) j};
    //                 surface_dual_edges.push_back(new_edge);
    //             }
    //         }
    //     }

    //     surface_dual.push_back(gg3d::static_graph(surface_dual_edges,
    //                            (int) cell_neighbors.size(),
    //                            surface_dual_vertex_names));
    // }

    // // Construct augmented neighborhood induced subgraphs
    // for (int cell_id = 0; cell_id < N; ++cell_id) {
    //     vector<int> all_aug_neighbors(neighbors[cell_id]);
    //     for (auto & aug_neighbor : aug_neighbors[cell_id]) all_aug_neighbors.push_back(aug_neighbor);

    //     // Convert augmented neighbor IDs to names
    //     vector<string> all_aug_neighbor_names = {};
    //     for (auto & aug_neighbor_id : all_aug_neighbors) {
    //         if (aug_neighbor_id < 0) continue; // Ignore walls (negative IDs)
    //         all_aug_neighbor_names.push_back(g.vertex_name[aug_neighbor_id]);
    //     }

    //     aug_neighbor_graph.push_back(g.induced_subgraph(all_aug_neighbor_names));
    // }

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

        static_graph surface_poset_graph = static_graph(sp_edges, num_faces + num_edges + num_vertices, node_names);
        surface_poset_graphs.push_back(surface_poset_graph);
    }

    // Now that all member variables are initialized, create initial assignment
    generate_initial_assignment();
}

/** TODO: Implement attempt_flip with new conditions */
flip_status geograph3d::attempt_flip(size_t &cell_id, size_t &new_part) {
    // Validate new_part (must be between 1 and K, inclusive)
    if (new_part < 1 || new_part > K) {
        throw std::invalid_argument("New part must be between 1 and K, inclusive.");
    }

    /** Construct S_1, the set of shared faces, edges, and vertices 
     * between the cell to flip and its current (i.e., giving) zone */
    std::set<string> S_1_vertices;
    std::set<string> S_1_edges;
    std::set<string> S_1_faces;
    for (auto &aug_neighbor : aug_neighbors[cell_id]) {
        if (assignment[aug_neighbor.id] != assignment[cell_id]) continue;
        vector<size_t> shared_vertices = aug_neighbor.shared_vertices;
        for (auto &shared_vertex_id : shared_vertices) {
            string vertex_name_str = "v" + std::to_string(shared_vertex_id);
            S_1_vertices.insert(vertex_name_str);
        }
        vector<size_t> shared_edges = aug_neighbor.shared_edges;
        for (auto &shared_edge_id : shared_edges) {
            string edge_name_str = "e" + std::to_string(shared_edge_id);
            S_1_edges.insert(edge_name_str);
        }
        vector<size_t> shared_faces = aug_neighbor.shared_faces;
        for (auto &shared_face_id : shared_faces) {
            string face_name_str = "f" + std::to_string(shared_face_id);
            S_1_faces.insert(face_name_str);
        }
    }

    // Combine shared elements of each dimension to get S_1
    vector<string> S_1;

    for (auto &S_1_vertex : S_1_vertices) S_1.push_back(S_1_vertex);
    for (auto &S_1_edge : S_1_edges) S_1.push_back(S_1_edge);
    for (auto &S_1_face : S_1_faces) S_1.push_back(S_1_face);

    /** Construct X_v, the set of vertices/edges on the surface of cell v 
     * and the surface of v's current (giving) zone */
    vector<string> X_v;
    for (auto &vertex_name : S_1_vertices) {
        size_t vertex_index = std::stoi(vertex_name.substr(1, vertex_name.size()));
        if (cell_vertices[vertex_index].get_is_boundary()) {
            X_v.push_back(vertex_name);
        }
    }
    for (auto &edge_name : S_1_edges) {
        size_t edge_index = std::stoi(edge_name.substr(1, edge_name.size()));
        if (cell_edges[edge_index].get_is_boundary()) {
            X_v.push_back(edge_name);
        }
    }

    /** TODO: Construct Y_v, the set of nodes representing edges and vertices 
     * on the boundary of the shared surface */
    // Not clear how to do this part

    /** TODO: Check all five conditions, 
     * returning flip_status::fail_[i] 
     * if condition (i) fails. 
     */

    /** TODO: Check condition (1): S_1 ∩ X_v = Y_v */


    /** Check condition (2):  G_v[S_1] is connected */
    bool satisfies_condition_2 = surface_poset_graphs[cell_id].is_connected_subgraph(S_1);
    if (!satisfies_condition_2) return flip_status::fail_2;

    /** TODO: Check condition (3):  */
    bool satisfies_condition_3 = false;
    if (!satisfies_condition_3) return flip_status::fail_3;
    
    /** TODO: Check condition (4):  */
    bool satisfies_condition_4 = false;
    if (!satisfies_condition_4) return flip_status::fail_4;

    /** TODO: Check condition (5):  */
    bool satisfies_condition_5 = false;
    if (!satisfies_condition_5) return flip_status::fail_5;

    assignment[cell_id] = new_part; // Flip succeeded, so update assignment
    return flip_status::success;
} /** end attempt_flip */

} /** end namespace gg3d */

// // Quick test of cell_vertex constructor and "==" operator. 
// int main(int argc, char *argv[]) {
//     gg3d::cell_vertex v(0, 1.1, 2.2, 3.3);
//     gg3d::cell_vertex u(0, 1.1000000000000000000001, 2.200000000002, 3.3);
    
//     if (u == v) {
//         cout << "u and v match.\n";
//     } else {
//         cout << "u and v are different.\n";
//     }

//     cout << "Enter any response to close. ";
//     string response;
//     std::cin >> response;
//     cout << response;
// }


/** Uncomment below to run interactive program allowing user to give flips to attempt */
int main(int argc, char *argv[]) {
    string in_filename;
    int K;
    
    if (argc >= 3) {
        in_filename = argv[1];
        K = std::stoi(argv[2]);
    } else {
        string prog_name = argv[0];
        if (prog_name.empty()) prog_name = "<program name>";
        cout << "Usage: " << prog_name << " <input_filename> <number_of_parts>\n";
        return 0;
    }

    gg3d::geograph3d geograph = gg3d::geograph3d(in_filename, K);

    // Summarize geo-graph
    cout << "The given 3-D geo-graph contains " << geograph.num_cells() << " cells ";
    cout << "partitioned into " << geograph.num_parts() << " parts.\n";

    // Summarize cell adjacency graph
    size_t count_g_edges = 0;
    for (auto & adj_list : geograph.g.adjacency_list) count_g_edges += adj_list.size();
    count_g_edges = count_g_edges / 2; // Every edge is double-counted when summing adjacency list sizes
    cout << "The cell adjacency graph has " << geograph.g.size << " vertices and " << count_g_edges << " edges.\n\n";

    // Print cell part assignments
    cout << "Current assignments:\n";
    vector<size_t> current_assignment = geograph.get_assignment();
    for (size_t i = 0; i < current_assignment.size(); ++i) {
        cout << i << " to part " << current_assignment[i] << "\n";
    }

    // Spot-check surface poset graphs (check that of cell 0)
    gg3d::static_graph spg0 = geograph.surface_poset_graphs[0];
    cout << "Surface poset graph of cell 0 (" << spg0.size << " nodes) has edges:\n";
    cout << "\tEdge Src, Edge Dest\n";
    for (size_t i = 0; i < spg0.adjacency_list.size(); ++i) {
        string src_name = spg0.vertex_name[i];
        for (size_t j = 0; j < spg0.adjacency_list[i].size(); ++j) {
            string dest_name = spg0.vertex_name[spg0.adjacency_list[i][j]];
            cout << "\t" << src_name << ",\t" << dest_name << "\n";
        }
        cout << "\n";
    }
    cout << "\n";

    // Attempt flips    
    size_t cell_id;
    string response = "Y";
    while (response[0] == 'Y' || response[0] == 'y') {
        cout << "Enter the name of the unit/cell to be flipped: ";
        string name;
        std::cin >> name;
        cell_id = stoi(name);
        current_assignment = geograph.get_assignment();
        size_t part = current_assignment[cell_id];
        cout << "\nUnit " << name << " is currently assigned to part " << part << ".\n";
        cout << "Enter the name of its new part: ";
        std::cin >> name;
        size_t new_part = stoi(name);
        gg3d::flip_status result = geograph.attempt_flip(cell_id, new_part);
        if (result == gg3d::flip_status::success) {
            cout << "Flip was successful. Unit " << cell_id << " is now assigned to part " << new_part << ".\n";
        } else {
            cout << "Flip failed with status " << gg3d::flip_status_name[result] << ".\n";
        }
        
        cout << "Try another flip? (Y/N) ";
        std::cin >> response;
    }

    /** Uncomment below if attempt flips block above is commented out */
    // string response;
    // cout << "enter any text to close. ";
    // std::cin >> response;
    // cout << response;
} /** end main */
