// Create Voronoi cells for randomly generated seed points in a rectangular prism (box)
#include "voro++.hh"
#include <chrono>
#include <iostream>
using namespace voro;

const bool DEBUG = false;

// Define bounding box
const int x_min = 0;
const int x_max = 300;
const int y_min = 0; 
const int y_max = 300;
const int z_min = 0; 
const int z_max = 300; 

const int box_vol = (x_max - x_min) * (y_max - y_min) * (z_max - z_min);

// Define number of seed points (particles)
const int num_particles = 27000;

// Define parameters for container blocks
const int blocks_x = 5, blocks_y = 5, blocks_z = 5;
const int particles_per_block = ((box_vol / blocks_x) / blocks_y) / blocks_z;

// Define uniform(a, b) pseudorandom number generator
double uniform(double a, double b) {
    return a + (b - a) * (double(rand()) / RAND_MAX);
}

int main() {
    int index;
    double x, y, z;

    // Create a Voro++ container object, non-periodic in each dimension, 
    // with the specified dimensions and block sizes. 
    // Pre-allocate space for avg_count particles per block. 
    container con(x_min, x_max, y_min, y_max, z_min, z_max, blocks_x, 
            blocks_y, blocks_z, false, false, false, particles_per_block);

    // Add particle at random coordinates
    index = 0;

    while (index < num_particles) {
        x = uniform(x_min, x_max);
        y = uniform(y_min, y_max);
        z = uniform(z_min, z_max);

        con.put(index, x, y, z);
        if (DEBUG) printf("Point %d: (%g, %g, %g)\n", index, x, y, z);

        index++;
    }

	// Save particles in POV-Ray format
	// con.draw_particles_pov("random_points_box_p.pov");

	// Save Voronoi cells in POV-Ray format
	// con.draw_cells_pov("random_points_box_v.pov");

    // Write CSV output file
    const char * output_fname = "random_voronoi_27000.csv";
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