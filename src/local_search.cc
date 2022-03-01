/**
 * \file local_search.cc 
 * \brief Implementation of random local search for 3D partitioning. 
 */

#ifndef CELLS_TO_FLIPS_MULTIPLIER
#define CELLS_TO_FLIPS_MULTIPLIER 10
#endif

#ifndef DEBUG 
#define DEBUG 1
#define DEBUG_LOG "local_search_log.csv"
#define RANDOM_SEED 42
#define MILESTONE_COUNT_ATTEMPTS 1000
#endif

// Flag for pausing terminal before exiting (used when testing)
#ifndef PAUSE
#define PAUSE 0
#endif

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
using gg3d::flip_status;

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
    cout << "Current assignments:\nCell ID,Part ID\n";
    vector<size_t> current_assignment = geograph.get_assignment();
    for (size_t i = 0; i < current_assignment.size(); ++i) {
        if (current_assignment[i] != 1) {
            cout << i << "," << current_assignment[i] << "\n";
        }
    }

    // Build set of current boundary faces
    vector<size_t> boundary_faces;
    vector<size_t> old_boundary_faces;
    vector<gg3d::cell_face> cell_faces = geograph.get_cell_faces();
    for (size_t i = 0; i < cell_faces.size(); ++i) {
        if (cell_faces[i].get_is_boundary() && !cell_faces[i].get_is_outer_boundary()) {
            boundary_faces.push_back(i);
        }
    }

    // Print augmented neighborhood sizes
    cout << "\n";
    vector<size_t> aug_neighborhood_size_freq; 
    aug_neighborhood_size_freq.resize(100, 0);
    for (size_t i = 0; i < geograph.num_cells(); ++i) {
        size_t aug_neighborhood_size = geograph.count_aug_neighbors(i);
        aug_neighborhood_size_freq[aug_neighborhood_size] += 1;
    }
    cout << "Augmented neighborhood size frequencies\nSize\tCount\n";
    for (size_t i = 0; i < aug_neighborhood_size_freq.size(); ++i) {
        if (aug_neighborhood_size_freq[i] > 0) {
            cout << i << "\t" << aug_neighborhood_size_freq[i] << "\n";
        }
    }
    cout << "\n";

    /* Attempt a number of flips proportional to number of cells */
    size_t num_flip_attempts = CELLS_TO_FLIPS_MULTIPLIER * geograph.num_cells();
    size_t current_attempt = 0;

    /** Total time spent on attempt_flip, in microseconds, separated based on result */
    long total_gg3d_flip_time_us[6] = {0, 0, 0, 0, 0, 0}; // There are six statuses in gg3d::flip_status enum
    /** Number of flips with each possible result */
    size_t count_gg3d_flips_with_result[6] = {0, 0, 0, 0, 0, 0};
    /** Total elapsed time during while loop, rounded to nearest second */
    long total_time_seconds;

    /** Total time spent on check_flip_BFS, in microseconds, separated by failure (0) / success (1) */
    long total_BFS_flip_time_us[2] = {0, 0};
    /** Number of failures (0) and successes (1) on check_flip_BFS */
    size_t count_BFS_flips_with_result[2] = {0, 0};

    #if DEBUG
    srand(RANDOM_SEED);
    std::ofstream log_file;
    log_file.open(DEBUG_LOG);
    log_file << "Attempt Index,Cell ID,Giving Zone,Receiving Zone,BFS Flip Status,gg3d Flip Status\n";
    #endif

    /* TODO : Consider changing initial assignment, 
    or running for several flips aiming for balanced zones 
    before running timing analysis for random flips. */

    auto start_while = std::chrono::high_resolution_clock::now();
    // Randomly select from boundary_faces, then try a flip from smaller part to larger
    while (current_attempt < num_flip_attempts && !boundary_faces.empty()) {
        current_attempt++;

        #if DEBUG
        if (current_attempt % MILESTONE_COUNT_ATTEMPTS == 0) cout << "Flip #" << current_attempt << "\n";
        #endif

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

        #if DEBUG
        log_file << current_attempt << "," << cell_to_flip << "," << geograph.get_assignment(cell_to_flip) << "," << new_part << ",";
        // cout << current_attempt << "," << cell_to_flip << "," << geograph.get_assignment(cell_to_flip) << "," << new_part << ",\n";
        
        // if (cell_to_flip == 207) {
        // cout << "boundary face at index " << boundary_face_index << ": f" << cell_faces[boundary_face_index].id << " with vertices \n";
        // for (auto &vertex_id : cell_faces[boundary_face_index].vertices) {
        //     cout << "v" << vertex_id << " -- (" << cell_vertices[vertex_id].pos[0] << "," << cell_vertices[vertex_id].pos[1] << "," << cell_vertices[vertex_id].pos[2] << ")\n";
        // }
        // cout << "\n";
        // }
        
        #endif

        /* Check whether flip maintains contiguity with BFS, 
        but only make the flip if it maintains spherical zones 
        (i.e., attempt_flip succeeds). 
        Allows for side-by=side comparison of 3DGG and BFS, 
        both in terms of time and success/failure. */

        auto start = std::chrono::high_resolution_clock::now();
        bool bfs_result = geograph.check_flip_BFS(cell_to_flip, new_part);
        auto stop = std::chrono::high_resolution_clock::now();
        #if DEBUG
        log_file << (bfs_result ? "success" : "failure") << ",";
        // cout << (bfs_result ? "success" : "failure") << ",";
        #endif
        total_BFS_flip_time_us[bfs_result] += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
        count_BFS_flips_with_result[bfs_result]++;

        start = std::chrono::high_resolution_clock::now();
        result = geograph.attempt_flip(cell_to_flip, new_part);
        stop = std::chrono::high_resolution_clock::now();
        #if DEBUG
        log_file << gg3d::flip_status_string(result) << "\n";
        // cout << gg3d::flip_status_string(result) << "\n";
        #endif
        total_gg3d_flip_time_us[static_cast<size_t>(result)] += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
        count_gg3d_flips_with_result[static_cast<size_t>(result)]++;
        if (result == gg3d::flip_status::success) {
            if (!bfs_result) {
                cout << "unexpected scenario: BFS flip fails, gg3d flip succeeds\n";
            }

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
            // Remove failed face from boundary_faces until next successful flip
            boundary_faces.erase(boundary_faces.begin() + rand_index);
        }
        
    }
    auto stop_while = std::chrono::high_resolution_clock::now();
    total_time_seconds = std::chrono::duration_cast<std::chrono::seconds>(stop_while - start_while).count();

    cout << "Terminated after attempting " << current_attempt << " flips.\n\n";
    
    #if DEBUG
    log_file.close();
    #endif

    // Summarize timing
    cout << "Total elapsed time," << total_time_seconds << ",seconds\n";
    size_t total_calls_attempt_flip = current_attempt;
    cout << "Total calls to attempt_flip," << total_calls_attempt_flip << "\n";
    cout << "Average time spent in attempt_flip\n"; 
    for (int status_val = flip_status::fail_1; status_val <= flip_status::success; ++status_val) {
        flip_status status = static_cast<flip_status>(status_val);
        cout << "\tWith result " << gg3d::flip_status_string(status) << "," << total_gg3d_flip_time_us[status_val] * 1.0 / count_gg3d_flips_with_result[status_val] << ",microseconds over," << count_gg3d_flips_with_result[status_val] << ",samples.\n";
    }
    cout << "\n";

    long sum_total_gg3d_flip_times = 0;
    for (int i = 0; i < 6; ++i) sum_total_gg3d_flip_times += total_gg3d_flip_time_us[i];
    float average_time_gg3d_flip = sum_total_gg3d_flip_times * 1.0 / total_calls_attempt_flip;
    cout << "Average time spent in attempt_flip," << average_time_gg3d_flip << ",microseconds\n";

    size_t total_calls_check_flip_BFS = current_attempt;
    cout << "Total calls to check_flip_BFS," << total_calls_check_flip_BFS << "\n";
    cout << "Average time spent in check_flip_BFS\n"; 
    cout << "\tWith success," << total_BFS_flip_time_us[1] * 1.0 / count_BFS_flips_with_result[1] << ",microseconds over," << count_BFS_flips_with_result[1] << ",samples.\n";
    cout << "\tWith failure," << total_BFS_flip_time_us[0] * 1.0 / count_BFS_flips_with_result[0] << ",microseconds over," << count_BFS_flips_with_result[0] << ",samples.\n";
    cout << "\n";

    float average_time_check_flip_BFS = (total_BFS_flip_time_us[0] + total_BFS_flip_time_us[1]) * 1.0 / total_calls_check_flip_BFS;
    cout << "Average time spent in check_flip_BFS," << average_time_check_flip_BFS << ",microseconds.\n";

    // Summarize new parts (zones)
    cout << "\nNew part sizes\n";
    for (size_t i = 1; i <= geograph.num_parts(); ++i) {
        cout << i << "," << geograph.get_part_size(i) << "\n";
    }

    #if PAUSE
    cout << "\nEnter any string to exit: ";
    string response;
    std::cin >> response;
    cout << response;
    #endif

} /** end main */
