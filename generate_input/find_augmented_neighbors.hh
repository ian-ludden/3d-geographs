/**
 * Header file for find_augmented_neighbors.cc. 
 */
#ifndef FIND_AUGMENTED_NEIGHBORS_HH
#define FIND_AUGMENTED_NEIGHBORS_HH
#include <iostream>
#include <string>
#include <vector>
using std::cout;

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
                } while (next_paren_pos < next_comma_pos);
                
                pos = next_comma_pos;
            }
        }

        // Check for trailing comma with no data
        pos = m_line.size();
        m_data.push_back(pos);
    };
};

#endif /** FIND_AUGMENTED_NEIGHBORS_HH */