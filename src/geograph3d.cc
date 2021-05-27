/** 
 * \file geograph3d.cc
 * \brief Implementation of the 3-D geo-graph. */
#include "geograph3d.hh"
#include "staticgraph.hh"
#include "../generate_input/find_augmented_neighbors.hh"
#include "../generate_input/vector_utils.hh"
#include <bitset>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
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
 *  - Builds augmented neighborhood induced subgraphs for each cell. 
 * 
 * \param[in] (in_filename) String name of the input CSV file, 
 *                          such as that produced by find_augmented_neighbors.
 * \param[in] (num_parts) The number of parts in the partition (stored as K). 
 * */
geograph3d::geograph3d(const string in_filename, const int num_parts) 
    : K{num_parts} {
    csv_row row;
    std::ifstream in_file(in_filename);
    in_file >> row; // Read first row to get number of cells

    N = stoi(row[0]);
    assignment.resize(N, -1); // Initialize all assignments to dummy part, -1
    
    // Edges of cell adjacency graph
    vector<Edge> g_edges;

    in_file >> row; // Skip row of column headers

    int id;

    vector<vector<int>> neighbors; // List of lists of neighbor IDs; primary index is cell ID
    vector<vector<int>> aug_neighbors; // List of lists of augmented neighbor (but not neighbor) IDs; primary index is cell ID
    vector<vector<uint64_t>> faces; // List of lists of faces, stored as uint64_t to support bit operations for detecting when two faces share an edge. NB: Breaks if any cell has more than 64 vertices.

    neighbors.reserve(N);
    aug_neighbors.reserve(N);
    faces.reserve(N);

    // Read each row and store info for that cell
    while (in_file >> row) {
        vector<int> cell_neighbors; // Temp list of neighbor IDs
        vector<int> cell_aug_neighbors; // Temp list of augmented neighbor IDs
        vector<uint64_t> cell_faces; // Temp list of bit representations of faces

        id = stoi(row[0]);

        cell_neighbors = delimited_list_to_vector_of_int(row[1], ' ');
        neighbors.push_back(cell_neighbors);

        // Add new edges (i.e., to higher-ID neighbors)
        for (auto &neighbor : cell_neighbors) {
            if (id < neighbor) g_edges.push_back({id, neighbor});
        }

        cell_aug_neighbors = delimited_list_to_vector_of_int(row[2], ' ');
        aug_neighbors.push_back(cell_aug_neighbors);

        vector<string> face_tuples = delimited_tuples_to_vector_of_string(row[3]);
        cell_faces.reserve(face_tuples.size());

        for (auto &face_tuple : face_tuples) {
            string face_tuple_contents = face_tuple.substr(1, face_tuple.length() - 2); // Remove parentheses
            vector<int> face_vertices = delimited_list_to_vector_of_int(face_tuple_contents, ' ');
            
            uint64_t face_bit_repr = 0; // Bit representation of face, where bit i is 1 (0) if vertex i is (not) part of the face
            for (auto &face_vertex : face_vertices) {
                face_bit_repr += 1 << face_vertex; // Bit shift to make bit face_vertex 1
            }
            cell_faces.push_back(face_bit_repr);
        }

        faces.push_back(cell_faces);
    }

    // Construct cell adjacency graph
    g = static_graph(g_edges, N);

    // Construct surface dual graphs
    for (int cell_id = 0; cell_id < N; ++cell_id) {
        vector<int> cell_neighbors = neighbors[cell_id];
        vector<uint64_t> cell_faces = faces[cell_id];

        vector<Edge> surface_dual_edges;
        vector<string> surface_dual_vertex_names;
        
        for (size_t i = 0; i < cell_faces.size(); ++i) {
            surface_dual_vertex_names.push_back(std::to_string(cell_neighbors[i]));

            for (size_t j = i + 1; j < cell_faces.size(); ++j) {
                size_t count_shared_vertices = std::bitset<64>{cell_faces[i] & cell_faces[j]}.count();
                // Two faces on the surface share an edge if they share (exactly) two vertices. 
                if (count_shared_vertices == 2) {
                    Edge new_edge = {(int) i, (int) j};
                    surface_dual_edges.push_back(new_edge);
                }
            }
        }

        surface_dual.push_back(gg3d::static_graph(surface_dual_edges,
                               (int) cell_neighbors.size(),
                               surface_dual_vertex_names));
    }

    // TODO: 
    // Need to remember correspondence between 
    // faces and neighbors for the purpose of determining 
    // vertex subset(s) in surface dual that need to be checked 
    // for connectedness. 
    // Maybe create a surface_dual_graph class that extends static_graph 
    // and adds this information?
    // Or just use static_graph class, but vertex names are the IDs 
    // of the neighboring cells corresponding to the faces? 

    // TODO: 
    // Need to remember IDs of cells corresponding to vertices in induced subgraph
    // of augmented neighborhood. 
    // Could either create an induced subgraph function in static_graph class 
    // or do something here as pre-processing. 
}
}

int main(int argc, char *argv[]) {
    // // Test counting bits of uint64_t
    // uint64_t test = 0x0f0b000000000316;
    // size_t count = std::bitset<64>{test}.count();
    // cout << count << "\n"; // should be 4+3+2+1+2 = 12
    // cout << "all done\n";

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

    // Summarize cell adjacency graph
    size_t count_g_edges = 0;
    for (auto & adj_list : geograph.g.adjacency_list) count_g_edges += adj_list.size();
    count_g_edges = count_g_edges / 2; // Every edge is double-counted when summing adjacency list sizes
    cout << "The cell adjacency graph has " << geograph.g.size << " vertices and " << count_g_edges << " edges.\n";

    // Spot-check surface dual graphs
    int cell_id = 0;
    cout << "The surface dual graph for cell " << cell_id;
    cout << " has " << geograph.surface_dual[cell_id].size << " vertices and ";
    size_t count_edges = 0;
    for (auto & adj_list : geograph.surface_dual[cell_id].adjacency_list) count_edges += adj_list.size();
    count_edges = count_edges / 2; // Every edge is double-counted when summing adjacency list sizes
    cout << count_edges << " edges.\n";
    cout << "The names of the vertices of the surface dual graph for cell " << cell_id << " are";
    for (auto & name : geograph.surface_dual[cell_id].vertex_name) cout << " " << name;
    cout << "\n\n";

    cell_id = 111;
    cout << "The names of the vertices of the surface dual graph for cell " << cell_id << " are";
    for (auto & name : geograph.surface_dual[cell_id].vertex_name) cout << " " << name;
    cout << "\n\n";

    cout << "All done. Type anything and press ENTER to close.\n";
    string temp;
    std::cin >> temp;
    cout << "temp is " << temp << "\n";
}
