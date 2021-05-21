/**
 * Transform CSV output by Voro++ with cell information 
 * into format needed to construct 3-D geo-graph. 
 */
#include "find_augmented_neighbors.hh"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
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
    cout << (string) row[0];
    total_particles = stoi(row[0]);
    
    // Cell IDs (integers from 0 to total_particles - 1, inclusive)
    vector<int> id;
    id.reserve(total_particles);

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
        id[index] = stoi(row[0]);
        // neighbor_count[id[index]] = stoi(row[1]); // Don't really need neighbor counts
        faces[id[index]] = row[4];

        string neighbors_line = row[2];
        size_t pos = 0;
        size_t prev_pos = -1;
        while ((pos = neighbors_line.find(' ', pos)) != string::npos) {
            neighbors[id[index]].push_back(stoi(neighbors_line.substr(prev_pos + 1, pos)));
            prev_pos = pos;
            ++pos;
        }
        
        // Look for vertices
        // Takes O(V^2) time, added over all iterations of while loop, 
        // where V is the number of vertices. 
        string vertices_string = row[5];
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

    index = 0;
    while (index < total_particles) {
        int cell_id = id[index];
        out_file << cell_id << ",";

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