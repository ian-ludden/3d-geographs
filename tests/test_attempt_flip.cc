/**
 * Tests the attempt_flip function of the geograph3d class. 
 */
#include "test_attempt_flip.hh"
#include "../src/geograph3d.hh"
#include "../src/static_graph.hh"
#include "../generate_input/conversion_utils.hh"
#include "../generate_input/vector_utils.hh"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
using std::string;
using std::vector;

int main(int argc, char *argv[]) {
    /** Filename of converted voro++ output into format for constructing geograph3d object */
    string in_filename;
    
    if (argc >= 2) {
        in_filename = argv[1];
    } else {
        cout << "Missing/incorrect parameter(s).\n";
        string prog_name = argv[0];
        if (prog_name.empty()) prog_name = "<program name>";
        cout << "Usage: " << prog_name << " <input_filename>\n";
        return 0;
    }

    cout << in_filename << "\n";

    // Read test cases from CSV file
    csv_row row;
    string test_cases_filename = "test_attempt_flip_cases.csv";
    std::ifstream test_cases_file(test_cases_filename);

    test_cases_file >> row; // Read first line (number of test cases)
    int num_test_cases = stoi(row[0]);
    test_cases_file >> row; // Read (and ignore) column headers

    vector<bool> results; // true for passing test, false for failing test
    results.reserve(num_test_cases);

    int curr_test_case = 0;
    while(curr_test_case < num_test_cases) {
        test_cases_file >> row;
        string test_id = row[0];
        size_t num_parts = stoi(row[1]);
        size_t cell_id = stoi(row[2]);
        size_t new_part = stoi(row[3]);
        string expected_result_str = row[4]; // Given as "TRUE" or "FALSE"
        bool expected_result = (expected_result_str.compare("TRUE") == 0);
        
        vector<size_t> init_assignment;
        init_assignment.resize(27, num_parts); // 3x3x3 uniform grid used for all test cases

        for (size_t part = 1; part < num_parts; ++part) { // Don't need last part, can leave as default
            vector<int> part_cell_ids = delimited_list_to_vector_of_int(row[4 + part], ' ');
            for (auto & id : part_cell_ids) {
                init_assignment[id] = part;
            }
        }

        gg3d::geograph3d geograph = gg3d::geograph3d(in_filename, num_parts);
        geograph.set_assignment(init_assignment);

        bool result = geograph.attempt_flip(cell_id, new_part);
        results.push_back(result == expected_result);

        if (!results[curr_test_case]) {
            cout << "\tFailed test ID " << test_id << ".\n";
        }

        curr_test_case++;
    }
    
    int success_count = 0;
    for (size_t i = 0; i < results.size(); ++i) {
        if (results[i]) ++success_count;
    }

    cout << "All test cases complete. Passed: " << success_count << "/" << num_test_cases << "\n";
    string exit_entry;
    cout << "Enter any string to exit. ";
    std::cin >> exit_entry;
    cout << exit_entry;
}
