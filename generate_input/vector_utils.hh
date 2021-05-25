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
 * \param[in] (&tuple) pointer to a string object 
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
 * \param[in] (&tuple) pointer to a string object 
 *            containing an n-tuple, i.e., 
 *            (x_1, x_2, x_3, ..., x_n), 
 *            each entry of which is a real value. 
 * 
 * \return A vector<double> object of length n containing 
 *         the entries of the n-tuple, in order. 
 */
vector<double> tuple_to_vector_of_double(const string tuple);

#endif /** VECTOR_UTILS_HH */