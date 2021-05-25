/**
 * Utility functions for creating/manipulating vectors.
 */
#include "vector_utils.hh"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
using std::cout;

// TODO: Consider seeing if C++ supports parameterized types 
// so the tuple_to_vector_of_<type> functions can be combined. 

vector<int> tuple_to_vector_of_int(const string tuple) {
    vector<int> vec;
    
    size_t start = tuple.find('(');
    size_t end = tuple.find(')');
    size_t left_pos = start + 1;
    size_t right_pos;

    while (left_pos < end) {
        right_pos = tuple.find(',', left_pos); 
        if (right_pos > end) {
            right_pos = end;
        }

        // Construct entry from left_pos (inclusive) to right_pos (exclusive)
        int entry = std::stoi(tuple.substr(left_pos, right_pos - left_pos));
        vec.push_back(entry);
        left_pos = right_pos + 1;
    }

    return vec;
}

vector<double> tuple_to_vector_of_double(const string tuple) {
    vector<double> vec;

    size_t start = tuple.find('(');
    size_t end = tuple.find(')');

    if (start > tuple.length()) {
        std::cerr << "Tuple must contain '('.\n";
        throw std::invalid_argument("Invalid tuple string given to tuple_to_vector_of_double: Tuple must contain '('.\n");
    }
    if (end > tuple.length()) {
        std::cerr << "Tuple must contain ')'.\n";
        throw std::invalid_argument("Invalid tuple string given to tuple_to_vector_of_double: Tuple must contain ')'.\n");
    }

    size_t left_pos = start + 1;
    size_t right_pos = tuple.find(',', left_pos);

    while (left_pos < end && right_pos < end) {
        // Construct entry from left_pos (inclusive) to right_pos (exclusive)
        double entry = std::stod(tuple.substr(left_pos, right_pos - left_pos));
        vec.push_back(entry);
        left_pos = right_pos + 1;
        right_pos = tuple.find(',', left_pos); 
    }

    // Catch last entry (right_pos == end)
    double entry = std::stod(tuple.substr(left_pos, end - left_pos));
    vec.push_back(entry);

    return vec;
}