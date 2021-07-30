/** 
 * \file staticgraph.hh
 * \brief Header file for staticgraph. */

#ifndef STATICGRAPH_HH
#define STATICGRAPH_HH
#include <string>
#include <vector>

namespace gg3d {

struct Edge {
    // Undirected; by convention, src < dest
    size_t src, dest;
};

class static_graph {
public:
    /** List of adjacency lists */
    std::vector<std::vector<size_t>> adjacency_list;
    /** Number of vertices */
    size_t size;
    /** List of vertex names, indexed by ID */
    std::vector<std::string> vertex_name;

    /**
     * Default constructor; creates empty graph (size = 0).
     */
    static_graph();

    /** Constructor for static graph. 
     * Builds graph from list of edges and number of vertices. 
     * \param[in] (&edges) pointer to vector of Edge structs, using 0 to n-1 as vertex IDs.
     * \param[in] (number_of_vertices) number of vertices in the graph (a.k.a., 'n')
     * */
    static_graph(const std::vector<Edge> &edges, size_t number_of_vertices);
    
    /** Alternate constructor for static graph, 
     * allowing vertex names to be specified (rather than default of 0 to n-1). 
     * \param[in] (&edges) pointer to vector of Edge structs, using 0 to n-1 as vertex IDs.
     * \param[in] (number_of_vertices) number of vertices in the graph (a.k.a., 'n')
     * \param[in] (&vertex_names) pointer to vector of string names of vertices, in order. 
     * */
    static_graph(const std::vector<Edge> &edges, size_t number_of_vertices, const std::vector<std::string> &vertex_names);
    
    /**
     * Check whether the induced subgraph of the given list of vertices 
     * is connected, using a standard breadth-first search (BFS) 
     * from the first vertex in the list. 
     * \param[in] (vertices) pointer to vector of vertex indices (not their string names)
     * */
    bool is_connected_subgraph(std::vector<size_t> &vertices);

    /**
     * Builds and returns a static_graph object representing
     * the induced subgraph of the given set of vertices. 
     * 
     * \param[in] (vertex_names) vector<string> of vertex names
     *                           to include in induced subgraph
     * \return the new static_graph object representing the induced subgraph
     */
    static_graph induced_subgraph(std::vector<std::string> &vertex_names);

private:
    /**
     * Builds adjacency list from edges (vector<Edge>). 
     * 
     * \param[in] (edge) a vector<Edge> containing all edges
     */
    void build_adjacency_lists(const std::vector<Edge> &edges) {
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
