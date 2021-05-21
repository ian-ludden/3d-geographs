/**
 * Transform CSV output by Voro++ with cell information 
 * into format needed to construct 3-D geo-graph. 
 */
#include "find_augmented_neighbors.hh"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
using std::cout;
using std::string;
using std::vector;

/**
 * Tolerance to apply to each coordinate when checking whether
 * two vertices in \reals^3 are equal. 
 */
const double VERTEX_TOLERANCE = 1.0e-10;

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

    // Strings of face tuples, to be used to determine surface dual graphs
    vector<string> faces;
    faces.reserve(total_particles);

    // Vector of cell vertices, to be used to determine augmented neighbors.
    
    in_file >> row; // Skip headers
    int index = 0;

    while (in_file >> row) {
        id_to_index[stoi(row[0])] = index;
        // neighbor_count[id_to_index[index]] = stoi(row[1]); // Don't really need neighbor counts
        string faces_string = row[3];
        std::replace(faces_string.begin(), faces_string.end(), ',', ' ');
        faces.push_back(faces_string);

        string neighbors_line = row[2];
        size_t pos = 0;
        size_t prev_pos = -1;
        vector<int> neighbor_list;
        neighbor_list.reserve(stoi(row[1]));
        while ((pos = neighbors_line.find(' ', pos)) != string::npos) {
            neighbor_list.push_back(stoi(neighbors_line.substr(prev_pos + 1, pos)));
            prev_pos = pos;
            ++pos;
        }
        neighbors.push_back(neighbor_list);
        
        // Look for vertices
        // Takes O(V^2) time, added over all iterations of while loop, 
        // where V is the number of vertices. 
        string vertices_string = (string) row[4];
        //TODO : Continue here
        // use atof to convert strings to doubles; 
        // probably want to write some utility function to turn a tuple into a vector, 
        // since that keeps coming up. 

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
        // Print first aug neighbor before loop to avoid extra space
        out_file << aug_neighbor_list[aug_neighbor_index++]; // TODO Fix seg fault here
        while (aug_neighbor_index < aug_neighbor_list.size()) {
            out_file << " " << aug_neighbor_list[aug_neighbor_index];
        }

        out_file << ",";

        // Write list of faces with commas in tuples replaced by spaces
        out_file << faces[cell_index];

        out_file << "\n";

        ++index;
    }

    out_file.close();
}

int main() {
    cout << "test print\n";
    cout << "breakpoint on this line\n";
    cout << "past breakpoint\n";

    std::string in_filename = "uniform_grid_output_full.csv";
    std::string out_filename = "uniform_grid_output_aug_neighbors.csv";
    
    convert_output_csv(in_filename, out_filename);
}