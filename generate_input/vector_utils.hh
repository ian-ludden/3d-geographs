/**
 * Header file for vector_utils.cc. 
 */
#ifndef VECTOR_UTILS_HH
#define VECTOR_UTILS_HH
#include <iostream>
#include <string>
#include <vector>
using std::string;
using std::vector;

/**
 * Converts a tuple of integers from a string to a vector. 
 * 
 * \param[in] (tuple) a string object 
 *            containing an n-tuple, i.e., 
 *            (x_1, x_2, x_3, ..., x_n), 
 *            each entry of which is integral. 
 *            Whitespace before/after each entry is optional. 
 * 
 * \return A vector<int> object of length n containing 
 *         the entries of the n-tuple, in order. 
 */
vector<int> tuple_to_vector_of_int(const string tuple);

/**
 * Converts a tuple of double from a string to a vector. 
 * 
 * \param[in] (tuple) a string object 
 *            containing an n-tuple, i.e., 
 *            (x_1, x_2, x_3, ..., x_n), 
 *            each entry of which is a real value. 
 * 
 * \return A vector<double> object of length n containing 
 *         the entries of the n-tuple, in order. 
 */
vector<double> tuple_to_vector_of_double(const string tuple);

/**
 * Converts a space-delimited list of integers to a vector. 
 * 
 * \param[in] (list) a string object containing 
 *                   delimited integers. 
 *                   There should be no leading or trailing spaces 
 *                   and only one character between each integer. 
 * \param[in] (delim) the delimiter character (char) used between integers
 */
vector<int> delimited_list_to_vector_of_int(const string list, const char delim);

/**
 * Converts a single-character-delimited list of tuples, 
 * each surrounded by parentheses and having no nested parentheses, 
 * to a vector of strings of tuples, keeping the parentheses. 
 * 
 * \param[in] (list) a string object containing tuples delimited by some character. 
 *                   There should be no leading or trailing whitespace 
 *                   and only one delimiter character between each tuple. 
 */
vector<string> delimited_tuples_to_vector_of_string(const string list);

#endif /** VECTOR_UTILS_HH */