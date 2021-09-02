/**
 * \file local_search_BFS.cc 
 * \brief Implementation of random local search for 3-D partitioning, 
 * using BFS for contiguity verification (as opposed to 3-D geo-graphs 
 * for spherical verification, as in local_search.cc). 
 */
#include "geograph3d.hh"
#include <chrono>
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

    // Build set of current boundary faces
    vector<size_t> boundary_faces;
    vector<size_t> old_boundary_faces;
    // cout << "Zone boundary faces:\n";
    vector<gg3d::cell_face> cell_faces = geograph.get_cell_faces();
    for (size_t i = 0; i < cell_faces.size(); ++i) {
        if (cell_faces[i].get_is_boundary() && !cell_faces[i].get_is_outer_boundary()) {
            boundary_faces.push_back(i);
        }
    }

    size_t num_flip_attempts = 10000;
    size_t num_reverse_flip_attempts = 0;
    size_t current_attempt = 0;

    /** Total time spent on attempt_flip_BFS, in microseconds, separated by failure (0) / success (1) */
    long total_flip_verification_time_us[2] = {0, 0};
    /** Number of failures (0) and successes (1) on attempt_flip_BFS */
    size_t count_flips_with_result[2] = {0, 0};
    /** Total elapsed time during while loop, rounded to nearest second */
    long total_time_seconds;

    auto start_while = std::chrono::high_resolution_clock::now();
    // Randomly select from boundary_faces, then try a flip from smaller part to larger
    while (current_attempt < num_flip_attempts && !boundary_faces.empty()) {
        current_attempt++;

        if (current_attempt % 1000 == 0) cout << "Flip #" << current_attempt << "\n";

        auto rand_index = rand() % boundary_faces.size();
        size_t boundary_face_index = boundary_faces[rand_index];

        size_t cell_1 = cell_faces[boundary_face_index].first_cell;
        size_t cell_2 = cell_faces[boundary_face_index].second_cell;
        size_t part_1 = geograph.get_assignment(cell_1);
        size_t part_2 = geograph.get_assignment(cell_2);
        size_t part_1_size = geograph.get_part_size(part_1);
        size_t part_2_size = geograph.get_part_size(part_2);

        // Attempt flip across boundary face, preferring from smaller part to larger part
        bool success;
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

        auto start = std::chrono::high_resolution_clock::now();
        success = geograph.attempt_flip_BFS(cell_to_flip, new_part);
        auto stop = std::chrono::high_resolution_clock::now();
        total_flip_verification_time_us[success] += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
        count_flips_with_result[success]++;

        if (success) {
            // cout << "Flip was successful, cell," << cell_to_flip << ",is now assigned to part," << new_part << ".\n";
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
            num_reverse_flip_attempts++;
            // Try reverse flip
            if (cell_to_flip == cell_1) {
                cell_to_flip = cell_2;
                new_part = part_1;
            } else {
                cell_to_flip = cell_1;
                new_part = part_2;
            }

            start = std::chrono::high_resolution_clock::now();
            success = geograph.attempt_flip_BFS(cell_to_flip, new_part);
            stop = std::chrono::high_resolution_clock::now();
            total_flip_verification_time_us[success] += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
            count_flips_with_result[success]++;
            
            if (success) {
                // cout << "Reverse flip was successful, cell," << cell_to_flip << ",is now assigned to part," << new_part << ".\n";
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
                // cout << "Reverse flip failed.\n";
                // Remove failed face from boundary_faces until next successful flip
                boundary_faces.erase(boundary_faces.begin() + rand_index);
            }
        }
    }
    auto stop_while = std::chrono::high_resolution_clock::now();
    total_time_seconds = std::chrono::duration_cast<std::chrono::seconds>(stop_while - start_while).count();

    cout << "Terminated after attempting " << current_attempt << " flips.\n";

    // Summarize timing
    cout << "Total elapsed time:\t" << total_time_seconds << " seconds\n";
    size_t total_calls_attempt_flip = current_attempt + num_reverse_flip_attempts;
    cout << "Total calls to attempt_flip_BFS (including reverse flip attempts):\t" << total_calls_attempt_flip << "\n";
    cout << "Average time spent in attempt_flip_BFS:\n"; 
    cout << "\tWith success: " << total_flip_verification_time_us[1] * 1.0 / count_flips_with_result[1] << " microseconds, over " << count_flips_with_result[1] << " samples.\n";
    cout << "\tWith failure: " << total_flip_verification_time_us[0] * 1.0 / count_flips_with_result[0] << " microseconds, over " << count_flips_with_result[0] << " samples.\n";
    cout << "\n";

    float average_time_attempt_flip = (total_flip_verification_time_us[0] + total_flip_verification_time_us[1]) * 1.0 / total_calls_attempt_flip;
    cout << "Average time spent in attempt_flip_BFS:\t" << average_time_attempt_flip << " microseconds.\n";

    cout << "\nNew part sizes:\n";
    for (size_t i = 1; i <= geograph.num_parts(); ++i) {
        cout << i << "," << geograph.get_part_size(i) << "\n";
    }

    cout << "\nEnter any string to exit: ";
    string response;
    std::cin >> response;
    cout << response;

} /** end main */
