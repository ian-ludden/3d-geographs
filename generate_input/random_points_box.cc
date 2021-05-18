// Create Voronoi cells for randomly generated seed points in a rectangular prism (box)
#include "voro++.hh"
using namespace voro;

const bool DEBUG = false;

// Define bounding box
const int x_min = 0;
const int x_max = 10;
const int y_min = 0; 
const int y_max = 20;
const int z_min = 0; 
const int z_max = 10; 

const int box_vol = (x_max - x_min) * (y_max - y_min) * (z_max - z_min);

// Define number of seed points (particles)
const int num_particles = 20;

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
	con.draw_particles_pov("random_points_box_p.pov");

	// Save Voronoi cells in POV-Ray format
	con.draw_cells_pov("random_points_box_v.pov");

    // Print column headers of CSV output file
    FILE * fp = safe_fopen("random_points_box_output.csv", "w");
    fprintf(fp, "ID,Coordinates,No. Neighbors,Neighbors,Face Vertices\n");
    // Output Voronoi cells with neighbor info. 
    // See http://math.lbl.gov/voro++/doc/voro++_overview.pdf, 
    // Section 6, to interpret format string. 
    con.print_custom("%i,%q,%s,%n,%t", fp);
    fclose(fp);
}