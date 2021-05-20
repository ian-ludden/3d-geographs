/** 
 * \file staticgraph.cc
 * \brief Implementation of a static, undirected graph. */
#include "staticgraph.hh"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <queue>
#include <vector>
using namespace std;

namespace gg3d {

    struct Edge {
        // Undirected; by convention, src < dest
        int src, dest;
    };

    class static_graph {
    public:
        vector<vector<int>> adjacency_list; // List of adjacency lists
        int size; // Number of vertices
        vector<string> vertex_name; // List of vertex names, indexed by ID
        
        /** Constructor for static graph. 
         * Builds graph from list of edges and number of vertices. 
         * \param[in] (&edges) pointer to vector of Edge structs, using 0 to n-1 as vertex IDs.
         * \param[in] (number_of_vertices) number of vertices in the graph (a.k.a., 'n')
         * */
        static_graph(const vector<Edge> &edges, int number_of_vertices) {
            // Build strings of default names
            vector<string> vertex_names; 
            vertex_names.reserve(number_of_vertices); 
            for (int i = 0; i < number_of_vertices; ++i) {
                vertex_names.push_back(to_string(i)); 
            }

            static_graph(edges, number_of_vertices, vertex_names);
        };

        /** Alternate constructor for static graph, 
         * allowing vertex names to be specified (rather than default of 0 to n-1). 
         * \param[in] (&edges) pointer to vector of Edge structs, using 0 to n-1 as vertex IDs.
         * \param[in] (number_of_vertices) number of vertices in the graph (a.k.a., 'n')
         * \param[in] (&vertex_names) pointer to vector of string names of vertices, in order. 
         * */
        static_graph(const vector<Edge> &edges, int number_of_vertices, const vector<string> &vertex_names) {
            adjacency_list.resize(number_of_vertices);
            vertex_name.resize(number_of_vertices);
            size = number_of_vertices;

            // For each edge, add src/dest to each other's adj lists
            for (auto &edge: edges) {
                adjacency_list[edge.src].push_back(edge.dest);
                adjacency_list[edge.dest].push_back(edge.src);
            }

            // For each vertex, save its name to vertex_name
            for (int i = 0; i < number_of_vertices; ++i) {
                vertex_name[i] = vertex_names[i];
            }
        };

        /**
         * Check whether the induced subgraph of the given list of vertices 
         * is connected, using a standard breadth-first search (BFS) 
         * from the first vertex in the list. 
         * */
        bool is_connected_subgraph(vector<int> vertices) const {
            if (vertices.size() <= 0) return false;

            vector<bool> is_found(size, false);
            vector<bool> is_available(size, false);
            for (auto &vertex: vertices) is_available[vertex] = true;

            queue<int> frontier;
            is_found[vertices[0]] = true;
            frontier.push(vertices[0]);

            while (!frontier.empty()) {
                int current_vertex = frontier.front();
                frontier.pop();

                for (auto &neighbor_vertex: adjacency_list[current_vertex]) {
                    if (is_available[neighbor_vertex] && !is_found[neighbor_vertex]) {
                        is_found[neighbor_vertex] = true;
                        frontier.push(neighbor_vertex);
                    }
                }
            }

            // Did we find them all? 
            bool all_found = true;
            for (int i = 0; i < size; ++i) {
                // Need to negate XOR(is_available, is_found)
                all_found = all_found && !(is_available[i] != is_found[i]); 
            }

            return all_found;
        };
    };
}

std::string int_vector_to_string(const vector<int> &vector_of_ints) {
    std::stringstream ss;

    ss << "<";
    for (int i = 0; i < vector_of_ints.size(); ++i) {
        if (i != 0) {
            ss << ", ";
        }
        ss << vector_of_ints[i];
    }
    ss << ">";

    return ss.str();
};

/**
 * Test static_graph class on dual graph of cube surface
 * */
int main() {
    const vector<gg3d::Edge> edges({
        {0, 1}, {0, 2}, {0, 3}, {0, 4}, 
        {1, 5}, {2, 5}, {3, 5}, {4, 5},
        {1, 2}, {2, 3}, {3, 4}, {1, 4}
    });

    int sides_of_cube = 6;
    vector<string> face_names = {"bottom", "left", "back", "right", "front", "top"};

    gg3d::static_graph graph = gg3d::static_graph(edges, sides_of_cube, face_names);

    for (int i = 0; i < sides_of_cube; ++i) {
        cout << "i = " << i << "\n";

        int n_neighbors = graph.adjacency_list[i].size();
        cout << "\tNo. neighbors: " << n_neighbors << "\n\tNeighbor names: ";

        for (int j = 0; j < n_neighbors; j++) {
            cout << graph.vertex_name[graph.adjacency_list[i][j]] << " ";
        }
        cout << "\n\tNeighbor indices: ";
        
        for (int j = 0; j < n_neighbors; j++) {
            cout << graph.adjacency_list[i][j] << " ";
        }
        cout << "\n";
    }

    cout << "\n";

    cout << "CHECK INDUCED SUBGRAPHS\n\n";
    // Test is_connected_subgraph on several vertex subsets
    vector<vector<int>> vertex_subsets = {
        {0, 1, 2}, {0, 3, 4, 5}, {1, 2, 3, 4}, {0, 1, 2, 3, 4, 5}, // Should all give 'true'
        {}, {0, 5}, {1, 3}, {2, 4} // Should all give 'false'
    };

    for (int subset_index = 0; subset_index < vertex_subsets.size(); subset_index++) {
        cout << "It is "; 
        if (graph.is_connected_subgraph(vertex_subsets[subset_index])) {
            cout << "true";
        } else {
            cout << "false";
        } 
        cout << " that the subgraph induced by "; 
        cout << int_vector_to_string(vertex_subsets[subset_index]);
        cout << " is connected.\n\n";
    }

    cout << "CHECK COMPLEMENTS\n\n";

    // Test whether the complements of vertex_subsets 
    // induce connected subgraphs. 
    vector<int> vertex_indices(graph.size);
    for (int i = 0; i < graph.size; ++i) vertex_indices.push_back(i); 

    for (int subset_index = 0; subset_index < vertex_subsets.size(); subset_index++) {
        vector<int> complement; 
        for (int i = 0; i < graph.size; ++i) {
            vector<int>::iterator it = find(vertex_subsets[subset_index].begin(), vertex_subsets[subset_index].end(), i);
            if (it == vertex_subsets[subset_index].end()) complement.push_back(i);
        }

        cout << "It is "; 
        if (graph.is_connected_subgraph(complement)) {
            cout << "true";
        } else {
            cout << "false";
        } 
        cout << " that the subgraph induced by "; 
        cout << int_vector_to_string(complement);
        cout << " is connected.\n\n";
    }
}