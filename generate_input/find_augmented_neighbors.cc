/**
 * Transform CSV output by Voro++ with cell information 
 * into format needed to construct 3-D geo-graph. 
 */
#include "find_augmented_neighbors.hh"
#include "vector_utils.hh"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
using std::cout;
using std::string;
using std::vector;

/**
 * Tolerance to apply to each coordinate when checking whether
 * two vertices in \reals^3 are equal. 
 */
const double VERTEX_TOLERANCE = 1.0e-10;

/**
 * Determines whether the two given vertices (as vectors of doubles 
 * representing each coordinate) are the same vertex, within a given tolerance.
 * 
 * \param[in] (v1) Vertex 1, given as a vector of `double` coordinates. 
 * \param[in] (v2) Vertex 2, given as a vector of `double` coordinates. 
 * \param[in] (tol) Absolute tolerance for considering two entries equal. 
 * 
 * \return true if the vertices agree (within tolerance tol) in every entry; false otherwise
 */
bool is_same_vertex(vector<double> v1, vector<double> v2, const double tol) {
    if (v1.size() != v2.size()) {
        return false; // Could also consider raising an exception
    }

    for (size_t i = 0; i < v1.size(); ++i) {
        if ((v1[i] > v2[i] + tol) || (v1[i] < v2[i] - tol)) {
            return false;
        }
    }

    return true;
}

void convert_output_csv(std::string in_filename, std::string out_filename) {
    csv_row row;
    int total_particles;

    std::ifstream in_file(in_filename);

    in_file >> row; // Read first line with number of cells
    total_particles = stoi(row[0]);
    
    // Cell IDs (integers from 0 to total_particles - 1, inclusive)
    // mapped to indices for other vectors
    vector<int> id_to_index;
    id_to_index.resize(total_particles);

    // Reverse map of id_to_index
    vector<int> index_to_id;
    index_to_id.reserve(total_particles);

    // Lists of neighbors IDs for each cell. 
    // Neighbors of cell with ID = id should be listed at 
    // neighbors[id_to_index[id]].
    vector<vector<int>> neighbors;
    neighbors.reserve(total_particles);

    // Lists of augmented neighbors (by ID) for each cell, indexed by ID
    vector<vector<int>> aug_neighbors;
    aug_neighbors.reserve(total_particles);

    // List of cell vertices, stored as 3-tuples of doubles
    vector<vector<double>> vertices;

    // List of vertex indices (by ID) for each cell, indexed by ID
    vector<vector<size_t>> vertex_indices;

    // Strings of face tuples, to be used to determine surface dual graphs
    vector<string> faces;
    faces.reserve(total_particles);

    // Vector of cell vertices, to be used to determine augmented neighbors.
    
    in_file >> row; // Skip headers
    int index = 0;

    while (in_file >> row) {
        // Read ID and update index <-> ID maps
        int id = stoi(row[0]);
        id_to_index[id] = index;
        index_to_id.push_back(id);

        // Read neighbors (space-delimited list)
        // int neighbor_count = stoi(row[1]); // No need to read number of neighbors
        string neighbors_line = row[2];

        // Parse neighbor list
        vector<int> neighbor_list = delimited_list_to_vector_of_int(neighbors_line, ' ');

        // Save neighbor list
        neighbors.push_back(neighbor_list);

        // Read and store faces string, but replace commas in tuples with spaces
        string faces_string = row[3];
        std::replace(faces_string.begin(), faces_string.end(), ',', ' ');
        faces.push_back(faces_string);

        // Parse and index the (unique) vertices, checking for duplicates. 
        string vertices_string = (string) row[4];
        vector<size_t> vertex_indices_vec;

        // Parse vertices
        size_t left_pos, right_pos;
        left_pos = vertices_string.find('(');
        while (left_pos < vertices_string.length()) {
            right_pos = vertices_string.find(')', left_pos);
            
            string cell_vertex_tuple = vertices_string.substr(left_pos, right_pos - left_pos + 1);
            vector<double> cell_vertex = tuple_to_vector_of_double(cell_vertex_tuple);

            size_t v_index = vertices.size();
            bool is_new = true;

            // Check whether this vertex has already been indexed
            for (size_t i = 0; i < vertices.size(); ++i) {
                if (is_same_vertex(vertices[i], cell_vertex, VERTEX_TOLERANCE)) {
                    v_index = i;
                    is_new = false;
                    break;
                }
            }
            vertex_indices_vec.push_back(v_index);
            if (is_new) vertices.push_back(cell_vertex);
            left_pos = vertices_string.find('(', right_pos);
        }

        vertex_indices.push_back(vertex_indices_vec);

        // Use vertices to find augmented neighbors.
        // Loop over all other (previous) cells, then 
        // loop over all pairs of vertices from this cell and other.
        vector<int> aug_neighbors_of_current;
        for (int other_index = 0; other_index < index; ++other_index) {
            bool is_aug_neighbor = false;
            for (size_t i = 0; i < vertex_indices[index].size(); ++i) {
                for (size_t j = 0; j < vertex_indices[other_index].size(); ++j) {
                    is_aug_neighbor = is_aug_neighbor || (vertex_indices[index][i] == vertex_indices[other_index][j]);
                }
            }
            // Only add to augmented neighbors list if not already a neighbor
            if (is_aug_neighbor && !(std::count(neighbor_list.begin(), neighbor_list.end(), index_to_id[other_index]))) {
                aug_neighbors_of_current.push_back(index_to_id[other_index]);
                aug_neighbors[other_index].push_back(index_to_id[index]); // Add aug neighbor relationship in both directions
            }
        }
        aug_neighbors.push_back(aug_neighbors_of_current);

        ++index;
    }
    in_file.close();


    // Write output CSV with augmented neighbors included
    std::ofstream out_file(out_filename);

    out_file << total_particles << "\n"; // Number of cells/particles
    out_file << "ID,Neighbors,Aug. Neighbors,Faces\n";

    int id = 0;
    while (id < total_particles) {
        int cell_index = id_to_index[id];
        out_file << id << ",";
        
        vector<int> neighbor_list = neighbors[cell_index];
        vector<int>::size_type neighbor_index = 0;
        // Print first neighbor first to avoid extra space (every cell has at least one neighbor)
        out_file << neighbor_list[neighbor_index++];
        while (neighbor_index < neighbor_list.size()) {
            out_file << " " << neighbor_list[neighbor_index++];
        }
        
        out_file << ",";
        
        // Write list of (additional) augmented neighbors
        vector<int> aug_neighbor_list = aug_neighbors[cell_index];
        vector<int>::size_type aug_neighbor_index = 0;

        if (!aug_neighbor_list.empty()) {
            // Print first aug neighbor before loop to avoid extra space
            out_file << aug_neighbor_list[aug_neighbor_index++];
            while (aug_neighbor_index < aug_neighbor_list.size()) {
                out_file << " " << aug_neighbor_list[aug_neighbor_index++];
            }
        }

        out_file << ",";

        // Write list of faces with commas in tuples replaced by spaces
        out_file << faces[cell_index];

        out_file << "\n";

        ++id;
    }

    out_file.close();
}

int main(int argc, char *argv[]) {
    string in_filename;
    string out_filename;

    if (argc >= 3) {
        in_filename = argv[1];
        out_filename = argv[2];
    } else {
        string prog_name = argv[0];
        if (prog_name.empty()) prog_name = "<program name>";
        cout << "Usage: " << prog_name << " <input_filename> <output_filename>\n";
        return 0;
    }

    convert_output_csv(in_filename, out_filename);
}