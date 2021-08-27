/**
 * \file local_search.cc 
 * \brief Implementation of random local search for 3-D partitioning. 
 */
#include "geograph3d.hh"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
using std::string;
using std::vector;
using std::cout;

int main(int argc, char *argv[]) {
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

    // Summarize geo-graph
    cout << "The given 3-D geo-graph contains " << geograph.num_cells() << " cells ";
    cout << "partitioned into " << geograph.num_parts() << " parts.\n";

    // Summarize cell adjacency graph
    size_t count_g_edges = 0;
    for (auto & adj_list : geograph.g.adjacency_list) count_g_edges += adj_list.size();
    count_g_edges = count_g_edges / 2; // Every edge is double-counted when summing adjacency list sizes
    cout << "The cell adjacency graph has " << geograph.g.size << " vertices and " << count_g_edges << " edges.\n\n";

    // Print cell part assignments
    cout << "Current assignments:\n";
    vector<size_t> current_assignment = geograph.get_assignment();
    for (size_t i = 0; i < current_assignment.size(); ++i) {
        cout << i << " to part " << current_assignment[i] << "\n";
    }

    // Build set of current boundary faces, printing each
    vector<size_t> boundary_faces;
    vector<size_t> old_boundary_faces;
    cout << "Zone boundary faces:\n";
    vector<gg3d::cell_face> cell_faces = geograph.get_cell_faces();
    for (size_t i = 0; i < cell_faces.size(); ++i) {
        if (cell_faces[i].get_is_boundary() && !cell_faces[i].get_is_outer_boundary()) {
            boundary_faces.push_back(i);
            cout << i << ", between cells " << cell_faces[i].first_cell << " and " << cell_faces[i].second_cell << ".\n";
        }
    }

    size_t num_flip_attempts = 1000;
    size_t current_attempt = 0;

    // Randomly select from boundary_faces, then try a flip from smaller part to larger
    while (current_attempt < num_flip_attempts && !boundary_faces.empty()) {
        current_attempt++;

        auto rand_index = rand() % boundary_faces.size();
        size_t boundary_face_index = boundary_faces[rand_index];

        size_t cell_1 = cell_faces[boundary_face_index].first_cell;
        size_t cell_2 = cell_faces[boundary_face_index].second_cell;
        size_t part_1 = geograph.get_assignment(cell_1);
        size_t part_2 = geograph.get_assignment(cell_2);
        size_t part_1_size = geograph.get_part_size(part_1);
        size_t part_2_size = geograph.get_part_size(part_2);

        // Attempt flip across boundary face, preferring from smaller part to larger part
        gg3d::flip_status result;
        size_t cell_to_flip;
        size_t new_part;
        if (part_1_size < part_2_size) {
            cell_to_flip = cell_2;
            new_part = part_1;
        } else if (part_2_size < part_1_size) {
            cell_to_flip = cell_1;
            new_part = part_2;
        } else { // Equal sizes, so toss a coin
            bool is_heads = (rand() % 2) == 0;
            if (is_heads) {
                cell_to_flip = cell_1;
                new_part = part_2;
            } else {
                cell_to_flip = cell_2;
                new_part = part_1;
            }
        }

        cout << "Attempting flip, cell:," << cell_to_flip << ",part:," << new_part << ",";
        result = geograph.attempt_flip(cell_to_flip, new_part);
        if (result == gg3d::flip_status::success) {
            cout << "Flip was successful, cell," << cell_to_flip << ",is now assigned to part," << new_part << ".\n";
            // Update boundary_faces vector
            old_boundary_faces.clear();
            for (auto &boundary_face_id : boundary_faces) old_boundary_faces.push_back(boundary_face_id);
            boundary_faces.clear();
            cell_faces = geograph.get_cell_faces();
            for (size_t i = 0; i < cell_faces.size(); ++i) {
                if (cell_faces[i].get_is_boundary() 
                    && !cell_faces[i].get_is_outer_boundary()) {
                    boundary_faces.push_back(i);
                }
            }
        } else {
            cout << "Flip failed,status " << gg3d::flip_status_string(result) << ".\n";

            // Try reverse flip
            if (cell_to_flip == cell_1) {
                cell_to_flip = cell_2;
                new_part = part_1;
            } else {
                cell_to_flip = cell_1;
                new_part = part_2;
            }
            cout << "Attempting reverse flip, cell:," << cell_to_flip << ",part:," << new_part << ",";

            result = geograph.attempt_flip(cell_to_flip, new_part);
            if (result == gg3d::flip_status::success) {
                cout << "Reverse flip was successful, cell," << cell_to_flip << ",is now assigned to part," << new_part << ".\n";
                // Update boundary_faces vector
                old_boundary_faces.clear();
            for (auto &boundary_face_id : boundary_faces) old_boundary_faces.push_back(boundary_face_id);
                boundary_faces.clear();
                cell_faces = geograph.get_cell_faces();
                for (size_t i = 0; i < cell_faces.size(); ++i) {
                    if (cell_faces[i].get_is_boundary() 
                        && !cell_faces[i].get_is_outer_boundary()) {
                        boundary_faces.push_back(i);
                    }
                }
            } else {
                cout << "Reverse flip failed,status " << gg3d::flip_status_string(result) << ".\n";
                // Remove failed face from boundary_faces until next successful flip
                boundary_faces.erase(boundary_faces.begin() + rand_index);
            }
        }
    }

    cout << "Terminated after attempting " << current_attempt << " flips.\n";

    vector<vector<size_t>> parts;
    parts.resize(geograph.num_parts(), {});

    cout << "\nNew assignments:\n";
    for (size_t i = 0; i < geograph.num_cells(); ++i) {
        size_t part_id = geograph.get_assignment(i);
        cout << i << "," << part_id << "\n";
        parts[part_id - 1].push_back(i);
    }

    cout << "\nNew parts:\n";
    for (size_t i = 0; i < parts.size(); ++i) {
        if (parts[i].empty()) {
            cout << i+1 << ",empty\n";
            continue;
        }
        cout << "Part " << i+1 << "," << parts[i][0];
        for (size_t j = 1; j < parts[i].size(); ++j) cout << "," << parts[i][j];
        cout << "\n";
    }

    cout << "\nNew part sizes:\n";
    for (size_t i = 1; i <= geograph.num_parts(); ++i) {
        cout << i << "," << geograph.get_part_size(i) << "\n";
    }

    cout << "\nEnter any string to exit: ";
    string response;
    std::cin >> response;
    cout << response;

} /** end main */
