/**
 * Transform CSV output by Voro++ with cell information 
 * into format needed to construct 3-D geo-graph. 
 */
#include "find_augmented_neighbors.hh"
#include "vector_utils.hh"
#include <fstream>
#include <iostream>
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

std::istream & operator>>(std::istream & str, csv_row & data) {
    data.read_next_row(str);
    return str;
}

void convert_output_csv(std::string in_filename, std::string out_filename) {
    csv_row row;
    int total_particles;

    std::ifstream in_file(in_filename);

    in_file >> row; // Read first line with number of cells
    total_particles = stoi(row[0]);
    cout << "There are " << total_particles << " cells.\n";
    
    // Cell IDs (integers from 0 to total_particles - 1, inclusive)
    // mapped to indices for other vectors
    vector<int> id_to_index;
    id_to_index.resize(total_particles);

    // Reverse map of id_to_index
    vector<int> index_to_id;
    index_to_id.reserve(total_particles);

    // Number of neighbors for each cell, indexed by ID
    // int neighbor_count[total_particles];

    // Lists of neighbors (by ID) for each cell, indexed by ID
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
        int id = stoi(row[0]);
        id_to_index[id] = index;
        index_to_id.push_back(id);
        int neighbor_count = stoi(row[1]);
        string faces_string = row[3];
        std::replace(faces_string.begin(), faces_string.end(), ',', ' ');
        faces.push_back(faces_string);

        string neighbors_line = row[2];
        size_t pos = 0;
        size_t prev_pos = -1;
        vector<int> neighbor_list;
        neighbor_list.reserve(neighbor_count);
        while ((pos = neighbors_line.find(' ', pos)) != string::npos) {
            int neighbor_id = stoi(neighbors_line.substr(prev_pos + 1, pos));
            if (neighbor_id < 0) neighbor_id = -1; // Collapse all walls to -1
            neighbor_list.push_back(neighbor_id);
            prev_pos = pos;
            ++pos;
        }
        // Catch last neighbor in space-separated list
        int neighbor_id = stoi(neighbors_line.substr(prev_pos + 1, neighbors_line.length() - prev_pos));
        if (neighbor_id < 0) neighbor_id = -1; // Collapse all walls to -1
        neighbor_list.push_back(neighbor_id);
        neighbors.push_back(neighbor_list);
        
        // Index the (unique) vertices, checking for duplicates. 
        // Takes O(V^2) time, added over all iterations of while loop, 
        // where V is the number of vertices. 
        string vertices_string = (string) row[4];

        vector<size_t> vertex_indices_vec;

        size_t left_pos = vertices_string.find('('); 
        size_t right_pos;
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

            if (is_new) {
                vertices.push_back(cell_vertex);
            }

            left_pos = vertices_string.find('(', right_pos);
        }

        vertex_indices.push_back(vertex_indices_vec);

        // Use vertices to find augmented neighbors
        // Loop over all previous cells, then all pairs of vertices from this cell and previous
        vector<int> aug_neighbors_of_current;
        for (int other_index = 0; other_index < index; ++other_index) {
            bool is_aug_neighbor = false;
            for (size_t i = 0; i < vertex_indices[index].size(); ++i) {
                for (size_t j = 0; j < vertex_indices[j].size(); ++j) {
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
        out_file << cell_index << ",";
        
        vector<int> neighbor_list = neighbors[cell_index];
        std::vector<int>::size_type neighbor_index = 0;
        // Print first neighbor first to avoid extra space
        out_file << neighbor_list[neighbor_index++];
        while (neighbor_index < neighbor_list.size()) {
            out_file << " " << neighbor_list[neighbor_index++];
        }
        
        out_file << ",";
        
        // Write list of (additional) augmented neighbors
        vector<int> aug_neighbor_list = aug_neighbors[cell_index];
        std::vector<int>::size_type aug_neighbor_index = 0;

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

int main() {
    vector<double> test_vec = tuple_to_vector_of_double("(10.4, 67.3543, 1935)");

    vector<double> vec1 = {0.333333333333333333, 0.2, 0.4};
    vector<double> vec2 = {0.333333333333333, 0.2, 0.4};
    if (is_same_vertex(vec1, vec2, VERTEX_TOLERANCE)) {
        cout << "vec1 and vec2 are the same\n";
    } else {
        cout << "vec1 and vec2 are different.\n";
    }

    std::string in_filename = "uniform_grid_output_full.csv";
    std::string out_filename = "uniform_grid_output_aug_neighbors.csv";
    
    convert_output_csv(in_filename, out_filename);
}