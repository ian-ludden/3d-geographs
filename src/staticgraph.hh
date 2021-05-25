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
    int src, dest;
};

class static_graph {
public:
    /** List of adjacency lists */
    std::vector<std::vector<int>> adjacency_list;
    /** Number of vertices */
    int size;
    /** List of vertex names, indexed by ID */
    std::vector<std::string> vertex_name;

    static_graph();

    /** Constructor for static graph. 
     * Builds graph from list of edges and number of vertices. 
     * \param[in] (&edges) pointer to vector of Edge structs, using 0 to n-1 as vertex IDs.
     * \param[in] (number_of_vertices) number of vertices in the graph (a.k.a., 'n')
     * */
    static_graph(const std::vector<Edge> &edges, int number_of_vertices);
    
    /** Alternate constructor for static graph, 
     * allowing vertex names to be specified (rather than default of 0 to n-1). 
     * \param[in] (&edges) pointer to vector of Edge structs, using 0 to n-1 as vertex IDs.
     * \param[in] (number_of_vertices) number of vertices in the graph (a.k.a., 'n')
     * \param[in] (&vertex_names) pointer to vector of string names of vertices, in order. 
     * */
    static_graph(const std::vector<Edge> &edges, int number_of_vertices, const std::vector<std::string> &vertex_names);
    
    /**
     * Check whether the induced subgraph of the given list of vertices 
     * is connected, using a standard breadth-first search (BFS) 
     * from the first vertex in the list. 
     * \param[in] (vertices) pointer to vector of vertex indices (not their string names)
     * */
    bool is_connected_subgraph(std::vector<int> &vertices);
};
    
}

#endif /* STATICGRAPH_HH */
