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

vector<int> delimited_list_to_vector_of_int(const string list, const char delim) {
    int num_entries = 1;
    for (size_t i = 0; i < list.length(); ++i) {
        if (delim == list[i]) num_entries++;
    }
    
    vector<int> vec;
    vec.reserve(num_entries);

    size_t left_pos = 0;
    size_t right_pos;

    for (int i = 0; i < num_entries; ++i) {
        right_pos = list.find(delim, left_pos);
        if (right_pos > list.length()) right_pos = list.length(); // Handle last entry in list

        int entry = stoi(list.substr(left_pos, right_pos - left_pos)); // Substring excludes character at index right_pos
        vec.push_back(entry);
        
        left_pos = right_pos + 1;
    }

    return vec;
}

vector<string> delimited_tuples_to_vector_of_string(const string list) {
    // Determine number of entries (tuples) by counting right parentheses
    int num_entries = 0;
    for (size_t i = 0; i < list.length(); ++i) {
        if (')' == list[i]) num_entries++;
    }

    vector<string> vec;
    vec.reserve(num_entries);

    size_t left_pos = 0;
    size_t right_pos;

    for (int i = 0; i < num_entries; ++i) {
        right_pos = list.find(')', left_pos);
        if (right_pos > list.length()) right_pos = list.length() - 1;
        
        string entry = list.substr(left_pos, right_pos - left_pos + 1); // Substring *includes* character at index right_pos
        vec.push_back(entry);

        left_pos = right_pos + 2; // Move from ')' to next '(', skipping delimiter (space character) in between
    }

    return vec;
}