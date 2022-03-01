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

    #if DEBUG
    srand(RANDOM_SEED);
    std::ofstream log_file;
    log_file.open(DEBUG_LOG);

    // Summarize geo-graph
    log_file << "The given 3-D geo-graph contains " << geograph.num_cells() << " cells ";
    log_file << "partitioned into " << geograph.num_parts() << " parts.\n";

    // Summarize cell adjacency graph
    size_t count_g_edges = 0;
    for (auto & adj_list : geograph.g.adjacency_list) count_g_edges += adj_list.size();
    count_g_edges = count_g_edges / 2; // Every edge is double-counted when summing adjacency list sizes
    log_file << "The cell adjacency graph has " << geograph.g.size << " vertices and " << count_g_edges << " edges.\n\n";

    // Print cell part assignments
    log_file << "Current assignments:\nCell ID,Part ID\n";
    vector<size_t> current_assignment = geograph.get_assignment();
    for (size_t i = 0; i < current_assignment.size(); ++i) {
        if (current_assignment[i] != 1) {
            log_file << i << "," << current_assignment[i] << "\n";
        }
    }
    #endif

    // Build set of current boundary faces
    vector<size_t> boundary_faces;
    vector<size_t> old_boundary_faces;
    vector<gg3d::cell_face> cell_faces = geograph.get_cell_faces();
    for (size_t i = 0; i < cell_faces.size(); ++i) {
        if (cell_faces[i].get_is_boundary() && !cell_faces[i].get_is_outer_boundary()) {
            boundary_faces.push_back(i);
        }
    }

    #if DEBUG
    // Print augmented neighborhood sizes
    log_file << "\n";
    vector<size_t> aug_neighborhood_size_freq; 
    aug_neighborhood_size_freq.resize(100, 0);
    for (size_t i = 0; i < geograph.num_cells(); ++i) {
        size_t aug_neighborhood_size = geograph.count_aug_neighbors(i);
        aug_neighborhood_size_freq[aug_neighborhood_size] += 1;
    }
    log_file << "Augmented neighborhood size frequencies\nSize\tCount\n";
    for (size_t i = 0; i < aug_neighborhood_size_freq.size(); ++i) {
        if (aug_neighborhood_size_freq[i] > 0) {
            log_file << i << "\t" << aug_neighborhood_size_freq[i] << "\n";
        }
    }
    log_file << "\n";
    #endif

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
        if (current_attempt % MILESTONE_COUNT_ATTEMPTS == 0) log_file << "Flip #" << current_attempt << "\n";
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
        #endif
        total_BFS_flip_time_us[bfs_result] += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
        count_BFS_flips_with_result[bfs_result]++;

        start = std::chrono::high_resolution_clock::now();
        result = geograph.attempt_flip(cell_to_flip, new_part);
        stop = std::chrono::high_resolution_clock::now();
        #if DEBUG
        log_file << gg3d::flip_status_string(result) << "\n";
        #endif
        total_gg3d_flip_time_us[static_cast<size_t>(result)] += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
        count_gg3d_flips_with_result[static_cast<size_t>(result)]++;
        if (result == gg3d::flip_status::success) {
            if (!bfs_result) {
                cout << "Unexpected scenario: BFS flip fails, gg3d flip succeeds\n";
                exit(-1);
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

    #if DEBUG
    log_file << "Terminated after attempting " << current_attempt << " flips.\n\n";
    
    // Summarize timing
    log_file << "Total elapsed time," << total_time_seconds << ",seconds\n";
    size_t total_calls_attempt_flip = current_attempt;
    log_file << "Total calls to attempt_flip," << total_calls_attempt_flip << "\n";
    log_file << "Average time spent in attempt_flip\n"; 
    for (int status_val = flip_status::fail_1; status_val <= flip_status::success; ++status_val) {
        flip_status status = static_cast<flip_status>(status_val);
        log_file << "\tWith result " << gg3d::flip_status_string(status) << "," << total_gg3d_flip_time_us[status_val] * 1.0 / count_gg3d_flips_with_result[status_val] << ",microseconds over," << count_gg3d_flips_with_result[status_val] << ",samples.\n";
    }
    log_file << "\n";

    long sum_total_gg3d_flip_times = 0;
    for (int i = 0; i < 6; ++i) sum_total_gg3d_flip_times += total_gg3d_flip_time_us[i];
    float average_time_gg3d_flip = sum_total_gg3d_flip_times * 1.0 / total_calls_attempt_flip;
    log_file << "Average time spent in attempt_flip," << average_time_gg3d_flip << ",microseconds\n";

    size_t total_calls_check_flip_BFS = current_attempt;
    log_file << "Total calls to check_flip_BFS," << total_calls_check_flip_BFS << "\n";
    log_file << "Average time spent in check_flip_BFS\n"; 
    log_file << "\tWith success," << total_BFS_flip_time_us[1] * 1.0 / count_BFS_flips_with_result[1] << ",microseconds over," << count_BFS_flips_with_result[1] << ",samples.\n";
    log_file << "\tWith failure," << total_BFS_flip_time_us[0] * 1.0 / count_BFS_flips_with_result[0] << ",microseconds over," << count_BFS_flips_with_result[0] << ",samples.\n";
    log_file << "\n";

    float average_time_check_flip_BFS = (total_BFS_flip_time_us[0] + total_BFS_flip_time_us[1]) * 1.0 / total_calls_check_flip_BFS;
    log_file << "Average time spent in check_flip_BFS," << average_time_check_flip_BFS << ",microseconds.\n";

    // Summarize new parts (zones)
    log_file << "\nNew part sizes\n";
    for (size_t i = 1; i <= geograph.num_parts(); ++i) {
        log_file << i << "," << geograph.get_part_size(i) << "\n";
    }

    log_file.close();
    #endif


    // Write timing summary to cout
    cout << "number of cells,flip type,status,average microseconds,number of samples\n";

    for (int status_val = flip_status::fail_1; status_val <= flip_status::success; ++status_val) {
        flip_status status = static_cast<flip_status>(status_val);
        cout << geograph.num_cells() << ",geograph3d," << gg3d::flip_status_string(status) << "," << total_gg3d_flip_time_us[status_val] * 1.0 / count_gg3d_flips_with_result[status_val] << "," << count_gg3d_flips_with_result[status_val] << "\n";
    }

    cout << geograph.num_cells() << ",bfs,failure," << total_BFS_flip_time_us[0] * 1.0 / count_BFS_flips_with_result[0] << "," << count_BFS_flips_with_result[0] << "\n";
    cout << geograph.num_cells() << ",bfs,success," << total_BFS_flip_time_us[1] * 1.0 / count_BFS_flips_with_result[1] << "," << count_BFS_flips_with_result[1] << "\n";


    #if PAUSE
    cout << "\nEnter any string to exit: ";
    string response;
    std::cin >> response;
    cout << response;
    #endif

} /** end main */
