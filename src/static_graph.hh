/** 
 * \file staticgraph.hh
 * \brief Header file for staticgraph. */

#ifndef STATICGRAPH_HH
#define STATICGRAPH_HH
#include <string>
#include <unordered_map>
#include <vector>
using std::string;
using std::unordered_map;
using std::vector;

namespace gg3d {
typedef unordered_map<string, size_t> name_map;

struct Edge {
    // Undirected; by convention, src < dest
    size_t src, dest;
};

class static_graph {
private:
    /** Map of vertex names (strings) to vertex IDs (local indices) */
    name_map map_names_to_ids;
public:
    /** List of adjacency lists */
    vector<vector<size_t>> adjacency_list;
    /** Number of vertices */
    size_t size;
    /** List of vertex names, indexed by ID */
    vector<string> vertex_name;

    /**
     * Default constructor; creates empty graph (size = 0).
     */
    static_graph();

    /** Constructor for static graph. 
     * Builds graph from list of edges and number of vertices. 
     * \param[in] (&edges) pointer to vector of Edge structs, using 0 to n-1 as vertex IDs.
     * \param[in] (number_of_vertices) number of vertices in the graph (a.k.a., 'n')
     * */
    static_graph(const vector<Edge> &edges, size_t number_of_vertices);
    
    /** Alternate constructor for static graph, 
     * allowing vertex names to be specified (rather than default of 0 to n-1). 
     * \param[in] (&edges) pointer to vector of Edge structs, using 0 to n-1 as vertex IDs.
     * \param[in] (number_of_vertices) number of vertices in the graph (a.k.a., 'n')
     * \param[in] (&vertex_names) pointer to vector of string names of vertices, in order. 
     * */
    static_graph(const vector<Edge> &edges, size_t number_of_vertices, const vector<string> &vertex_names);
    
    /**
     * Check whether the induced subgraph of the given list of vertices 
     * is connected, using a standard breadth-first search (BFS) 
     * from the first vertex in the list. 
     * \param[in] (vertices) pointer to vector of vertex indices (not their string names)
     */
    bool is_connected_subgraph(vector<size_t> &vertices);

    /**
     * Check whether the induced subgraph of the given list of vertex names 
     * is connected. First converts the names to indices, then calls 
     * the is_connected_subgraph function that takes vertex indices. 
     * \param[in] (subgraph_vertex_names) pointer to vector of vertex names (not their IDs/indices)
     */
    bool is_connected_subgraph(vector<string> &subgraph_vertex_names);

    /**
     * Builds and returns a static_graph object representing
     * the induced subgraph of the given set of vertices. 
     * 
     * \param[in] (vertex_names) vector<string> of vertex names
     *                           to include in induced subgraph
     * \return the new static_graph object representing the induced subgraph
     */
    static_graph induced_subgraph(vector<string> &vertex_names);

    /**
     * Converts a vertex name (string) 
     * to a vertex id (local index, of type size_t). 
     * 
     * \param[in] (vertex_name) string name of vertex
     * \return the size_t local index of the vertex with the given name
     */
    size_t get_index_of_name(string &vertex_name);

private:
    /**
     * Builds adjacency list from edges (vector<Edge>). 
     * 
     * \param[in] (edge) a vector<Edge> containing all edges
     */
    void build_adjacency_lists(const vector<Edge> &edges) {
        adjacency_list.resize(size, {}); // Initialize all lists empty

        // For each edge, add src/dest to each other's adj lists
        for (auto &edge: edges) {
            adjacency_list[edge.src].push_back(edge.dest);
            adjacency_list[edge.dest].push_back(edge.src);
        }
    }
};
    
}

#endif /* STATICGRAPH_HH */
