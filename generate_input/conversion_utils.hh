/**
 * Header file for conversion_utils.cc. 
 */
#ifndef CONVERSION_UTILS_HH
#define CONVERSION_UTILS_HH
#include <iostream>
#include <string>
#include <vector>
using std::cout;

/**
 * Tolerance to apply to each coordinate when checking whether
 * two vertices in \reals^3 are equal. 
 */
const double VERTEX_TOLERANCE = 1.0e-10;

/**
 * Object representation of CSV row, based on StackOverflow answer: 
 * https://stackoverflow.com/a/1120224 (CC-by-SA 4.0)
 * 
 * Adapted to avoid separating coordinate tuples, that is, 
 * any commas inside of parentheses are not treated as delimeters. 
 * 
 * Not guaranteed to handle nested parentheses properly. 
 */ 
class csv_row {
private:
    std::string m_line;
    std::vector<int> m_data;

public:
    std::string operator[](std::size_t index) const {
        return std::string(&m_line[m_data[index] + 1], m_data[index + 1] - (m_data[index] + 1));
    };

    int size() const {
        return m_data.size();
    }

    void read_next_row(std::istream & str) {
        std::getline(str, m_line);

        m_data.clear();
        m_data.push_back(-1);
        
        std::string::size_type pos = 0;

        while ((pos = m_line.find(',', pos)) != std::string::npos) {
            m_data.push_back(pos);
            ++pos;

            // Don't separate tuples of coordinates; 
            // check for parentheses, and skip over any commas inside them
            if (m_line[pos] == '(') {                
                // Skip to end of tuple (right parenthesis)
                pos = m_line.find(')', pos);

                // Skip other parentheses, too
                std::string::size_type next_paren_pos = m_line.find('(', pos);
                std::string::size_type next_comma_pos = m_line.find(',', pos);

                do {
                    pos = m_line.find(')', next_paren_pos);
                    next_paren_pos = m_line.find('(', pos);
                    next_comma_pos = m_line.find(',', pos);
                } while (next_paren_pos < next_comma_pos 
                         && next_paren_pos > 0);
                
                if (next_comma_pos > 0
                    && next_comma_pos < m_line.size()) {
                        pos = next_comma_pos;
                }
            }
        }

        // Check for trailing comma with no data
        pos = m_line.size();
        m_data.push_back(pos);
    };
};

std::istream & operator>>(std::istream & str, csv_row & data);

/**
 * Determines whether the two given vertices (as vectors of doubles 
 * representing each coordinate) are the same vertex, within a given tolerance.
 * 
 * \param[in] (v1) Vertex 1, given as a vector of `double` coordinates. 
 * \param[in] (v2) Vertex 2, given as a vector of `double` coordinates. 
 * \param[in] (tol) Absolute tolerance for considering two entries equal. 
 * 
 * \return true if the vertices agree (within tolerance tol) in every entry; false otherwise
 */
bool is_same_vertex(std::vector<double> v1, std::vector<double> v2, const double tol);

/**
 * Converts the output csv file (from random_points_box.cc, uniform_grid.cc, 
 * or another input generator with the same output format) into a version for 
 * building a 3-D geo-graph (geograph3d object). 
 * 
 * One of the primary subtasks is finding the augmented neighbors of each cell. 
 */
void convert_output_csv(std::string in_filename, std::string out_filename);

#endif /** CONVERSION_UTILS_HH */