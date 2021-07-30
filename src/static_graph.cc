/** 
 * \file staticgraph.cc
 * \brief Implementation of a static, undirected graph. */
#include "static_graph.hh"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <queue>
#include <vector>
using std::string;
using std::vector;

namespace gg3d {

struct Edge;

/** Default constructor for static graph. 
 * Builds an empty graph. 
 */
static_graph::static_graph() : size{0} {}

/** Constructor for static graph. 
 * Builds graph from list of edges and number of vertices. 
 * \param[in] (&edges) pointer to vector of Edge structs, using 0 to n-1 as vertex IDs.
 * \param[in] (number_of_vertices) number of vertices in the graph (a.k.a., 'n')
 * */
static_graph::static_graph(const vector<Edge> &edges, size_t number_of_vertices)
    : size{number_of_vertices} {
    build_adjacency_lists(edges);
    // Assign default vertex names: 0 to size-1
    vertex_name.reserve(number_of_vertices); 
    for (size_t i = 0; i < number_of_vertices; ++i) {
        vertex_name.push_back(std::to_string(i)); 
    }
}

/** Alternate constructor for static graph, 
 * allowing vertex names to be specified (rather than default of 0 to n-1). 
 * \param[in] (&edges) pointer to vector of Edge structs, using 0 to n-1 as vertex IDs.
 * \param[in] (number_of_vertices) number of vertices in the graph (a.k.a., 'n')
 * \param[in] (&vertex_names) pointer to vector of string names of vertices, in order. 
 * */
static_graph::static_graph(const vector<Edge> &edges, size_t number_of_vertices, const vector<string> &vertex_names)
    : size{number_of_vertices}, vertex_name{vertex_names} {
    build_adjacency_lists(edges);
}

/**
 * Check whether the induced subgraph of the given list of vertices 
 * is connected, using a standard breadth-first search (BFS) 
 * from the first vertex in the list. 
 * \param[in] (vertices) pointer to vector of vertex indices (not their string names)
 * */
bool static_graph::is_connected_subgraph(vector<size_t> &vertices) {
    if (vertices.size() <= 0) return false; // For our purposes, the empty graph is not connected

    vector<bool> is_found(this->size, false); // Indicator vector for visited vertices
    vector<bool> is_available(size, false); // Indicator vector for vertices in the induced subgraph
    for (auto &vertex: vertices) is_available[vertex] = true;

    // Run BFS starting from the first vertex
    std::queue<size_t> frontier;
    is_found[vertices[0]] = true;
    frontier.push(vertices[0]);

    while (!frontier.empty()) {
        size_t current_vertex = frontier.front();
        frontier.pop();

        // Add unvisited neighbors to the queue if in the induced subgraph
        for (auto &neighbor_vertex: adjacency_list[current_vertex]) {
            if (is_available[neighbor_vertex] && !is_found[neighbor_vertex]) {
                is_found[neighbor_vertex] = true;
                frontier.push(neighbor_vertex);
            }
        }
    }

    // Did we find them all? 
    bool all_found = true;
    for (size_t i = 0; i < size; ++i) {
        // Need to negate logical XOR(is_available, is_found)
        all_found = all_found && !(is_available[i] != is_found[i]); 
    }

    return all_found;
}

}

string int_vector_to_string(const vector<int> &vector_of_ints) {
    std::stringstream ss;

    ss << "<";
    for (size_t i = 0; i < vector_of_ints.size(); ++i) {
        if (i != 0) {
            ss << ", ";
        }
        ss << vector_of_ints[i];
    }
    ss << ">";

    return ss.str();
}

gg3d::static_graph gg3d::static_graph::induced_subgraph(vector<string> &vertex_names) {
    vector<gg3d::Edge> edges;
    
    vector<size_t> indices;
    for (auto & name : vertex_names) {
        vector<string>::iterator it = find(vertex_name.begin(), vertex_name.end(), name);
        if (it == vertex_name.end()) throw std::invalid_argument("Error in static_graph::induced_subgraph(): Vertex \"" + name + "\" not found.");
        indices.push_back(it - vertex_name.begin());
    }

    // Find edges to add from adjacency lists of vertices in induced subgraph
    for (size_t i = 0; i < indices.size(); ++i) {
        // i is the new index, indices[i] is the old index
        vector<size_t> adj_list_i = adjacency_list[indices[i]];
        for (auto & j_id_old : adj_list_i) {
            // Check whether j_id_old is a vertex in the induced subgraph. 
            // If so, determine its new id, and add the edge from i to j.
            vector<size_t>::iterator it = find(indices.begin(), indices.end(), j_id_old);
            if (it != indices.end()) {
                size_t j_id_new = it - indices.begin();
                gg3d::Edge edge = {i, j_id_new};
                edges.push_back(edge);
            }
        }
    }

    return gg3d::static_graph(edges, (int) vertex_names.size(), vertex_names);
}

// /**
//  * Test static_graph class on dual graph of cube surface
//  * */
// int main() {
//     const vector<gg3d::Edge> edges({
//         {0, 1}, {0, 2}, {0, 3}, {0, 4}, 
//         {1, 5}, {2, 5}, {3, 5}, {4, 5},
//         {1, 2}, {2, 3}, {3, 4}, {1, 4}
//     });

//     int sides_of_cube = 6;
//     std::vector<std::string> face_names = {"bottom", "left", "back", "right", "front", "top"};

//     gg3d::static_graph graph = gg3d::static_graph(edges, sides_of_cube, face_names);

//     for (int i = 0; i < sides_of_cube; ++i) {
//         cout << "i = " << i << "\n";

//         int n_neighbors = graph.adjacency_list[i].size();
//         cout << "\tNo. neighbors: " << n_neighbors << "\n\tNeighbor names: ";

//         for (int j = 0; j < n_neighbors; j++) {
//             cout << graph.vertex_name[graph.adjacency_list[i][j]] << " ";
//         }
//         cout << "\n\tNeighbor indices: ";
        
//         for (int j = 0; j < n_neighbors; j++) {
//             cout << graph.adjacency_list[i][j] << " ";
//         }
//         cout << "\n";
//     }

//     cout << "\n";

//     cout << "CHECK INDUCED SUBGRAPHS\n\n";
//     // Test is_connected_subgraph on several vertex subsets
//     std::vector<std::vector<int>> vertex_subsets = {
//         {0, 1, 2}, {0, 3, 4, 5}, {1, 2, 3, 4}, {0, 1, 2, 3, 4, 5}, // Should all give 'true'
//         {}, {0, 5}, {1, 3}, {2, 4} // Should all give 'false'
//     };

//     for (int subset_index = 0; subset_index < vertex_subsets.size(); subset_index++) {
//         cout << "It is "; 
//         if (graph.is_connected_subgraph(vertex_subsets[subset_index])) {
//             cout << "true";
//         } else {
//             cout << "false";
//         } 
//         cout << " that the subgraph induced by "; 
//         cout << int_vector_to_string(vertex_subsets[subset_index]);
//         cout << " is connected.\n\n";
//     }

//     cout << "CHECK COMPLEMENTS\n\n";

//     // Test whether the complements of vertex_subsets 
//     // induce connected subgraphs. 
//     std::vector<int> vertex_indices(graph.size);
//     for (int i = 0; i < graph.size; ++i) vertex_indices.push_back(i); 

//     for (int subset_index = 0; subset_index < vertex_subsets.size(); subset_index++) {
//         std::vector<int> complement; 
//         for (int i = 0; i < graph.size; ++i) {
//             std::vector<int>::iterator it = find(vertex_subsets[subset_index].begin(), vertex_subsets[subset_index].end(), i);
//             if (it == vertex_subsets[subset_index].end()) complement.push_back(i);
//         }

//         cout << "It is "; 
//         if (graph.is_connected_subgraph(complement)) {
//             cout << "true";
//         } else {
//             cout << "false";
//         } 
//         cout << " that the subgraph induced by "; 
//         cout << int_vector_to_string(complement);
//         cout << " is connected.\n\n";
//     }
// }