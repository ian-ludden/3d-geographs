// Create Voronoi cells for randomly generated seed points in a rectangular prism (box)
#include "voro++.hh"
#include <chrono>
#include <iostream>
#include <cmath>
#include <stdexcept>
using namespace voro;

const bool DEBUG = false;

// Define bounding box
const int x_min = 0;
const int x_max = 200;
const int y_min = 0; 
const int y_max = 200;
const int z_min = 0; 
const int z_max = 200; 

const int box_vol = (x_max - x_min) * (y_max - y_min) * (z_max - z_min);

// Define number of seed points (particles)
const int num_particles = 8000;

// Define minimum squared distance between particles
const double MIN_DISTANCE_SQUARED = 0.25;

// Define max number of failures allowed 
// (attempts to create a particle that is too close to another)
const int MAX_FAILURES = 1000000;

// Define parameters for container blocks
const int blocks_x = 5, blocks_y = 5, blocks_z = 5;
const int particles_per_block = ((box_vol / blocks_x) / blocks_y) / blocks_z;

// Define uniform(a, b) pseudorandom number generator
double uniform(double a, double b) {
    return a + (b - a) * (double(rand()) / RAND_MAX);
}

int main() {
    int index;
    int count_failures;
    double x[num_particles], y[num_particles], z[num_particles];

    std::cout << "Min distance squared: " << MIN_DISTANCE_SQUARED << "\n";
    std::cout << "Min distance: " << std::sqrt(MIN_DISTANCE_SQUARED) << "\n";

    // Create a Voro++ container object, non-periodic in each dimension, 
    // with the specified dimensions and block sizes. 
    // Pre-allocate space for avg_count particles per block. 
    container con(x_min, x_max, y_min, y_max, z_min, z_max, blocks_x, 
            blocks_y, blocks_z, false, false, false, particles_per_block);

    // Add particle at random coordinates
    index = 0;
    count_failures = 0;

    while (index < num_particles) {
        x[index] = uniform(x_min, x_max);
        y[index] = uniform(y_min, y_max);
        z[index] = uniform(z_min, z_max);

        // Make sure new point is not too close to any existing points
        bool is_too_close = false;
        for (size_t i = 0; i < index; ++i) {
            double dist_sq = (x[i] - x[index])*(x[i] - x[index]) + (y[i] - y[index])*(y[i] - y[index]) + (z[i] - z[index])*(z[i] - z[index]);
            if (dist_sq < MIN_DISTANCE_SQUARED) {
                is_too_close = true;
                std::cout << "Tried to create two points within a distance of " << std::sqrt(dist_sq) << ".\n";
                break;
            }
        }
        
        if (!is_too_close) {
            con.put(index, x[index], y[index], z[index]);
            if (DEBUG) printf("Point %d: (%g, %g, %g)\n", index, x, y, z);
            index++;
        } else {
            count_failures++;
        }

        if (count_failures > MAX_FAILURES) {
            throw std::runtime_error("Failed too many times to create a Voronoi seed far enough away from existing seeds.\n");
        }
    }

    // Determine minimum distance
    double min_dist_sq = box_vol;
    for (size_t i = 0; i < num_particles; ++i) {
        for (size_t j = i+1; j < num_particles; ++j) {
            double dist_sq = (x[i] - x[j])*(x[i] - x[j]) + (y[i] - y[j])*(y[i] - y[j]) + (z[i] - z[j])*(z[i] - z[j]);
            if (dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
            }
        }
    }
    std::cout << "Min distance squared (actual): " << min_dist_sq << "\n";
    std::cout << "Min distance (actual): " << std::sqrt(min_dist_sq) << "\n";

	// Save particles in POV-Ray format
	// con.draw_particles_pov("random_points_box_p.pov");

	// Save Voronoi cells in POV-Ray format
	// con.draw_cells_pov("random_points_box_v.pov");

    // Write CSV output file
    const char * output_fname = "random_voronoi_8000.csv";
    std::cout << "Writing Voronoi cells to " << output_fname << ".\n";
    auto start = std::chrono::high_resolution_clock::now();
    FILE * fp = safe_fopen(output_fname, "w");
    fprintf(fp, "%d\n", index); // Print number of cells on a line by itself
    fprintf(fp, "ID,No. Neighbors,Neighbors,Faces,Vertices\n");
    // Output Voronoi cells with neighbor info. 
    // See http://math.lbl.gov/voro++/doc/voro++_overview.pdf, 
    // Section 6, to interpret format string. 
    con.print_custom("%i,%s,%n,%t,%P", fp);
    fclose(fp);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Elapsed time: " << duration.count() << " milliseconds.\n";

    std::string response;
    std::cout << "Enter any string to exit: ";
    std::cin >> response;
    std::cout << response;
}