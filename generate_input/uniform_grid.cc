// Create Voronoi cells for a 3-D grid in a rectangular prism (box)
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

// Define parameters for container blocks
const int blocks_x = 5, blocks_y = 5, blocks_z = 5;
const int particles_per_block = ((box_vol / blocks_x) / blocks_y) / blocks_z;

int main() {
    int index;
    int i, j, k;
    double x, y, z;

    // Create a Voro++ container object, non-periodic in each dimension, 
    // with the specified dimensions and block sizes. 
    // Pre-allocate space for avg_count particles per block. 
    container con(x_min, x_max, y_min, y_max, z_min, z_max, blocks_x, 
            blocks_y, blocks_z, false, false, false, particles_per_block);

    // Add particle at center of each unit cube
    index = 0;
    for (i = x_min; i < x_max; i++) {
        x = i + 0.5;
        for (j = y_min; j < y_max; j++) {
            y = j + 0.5;
            for (k = z_min; k < z_max; k++) {
                z = k + 0.5;
                con.put(index, x, y, z);
                if (DEBUG) printf("Point %d: (%g, %g, %g)\n", index, x, y, z);
                index++;
            }
        }
    }

	// Save particles in POV-Ray format
	con.draw_particles_pov("uniform_grid_p.pov");

	// Save Voronoi cells in POV-Ray format
	con.draw_cells_pov("uniform_grid_v.pov");

    // Print column headers of CSV output file
    const char * full_output_fname = "uniform_grid_output_full.csv";
    FILE * fp = safe_fopen(full_output_fname, "w");
    fprintf(fp, "%d\n", index); // Print number of cells on a line by itself
    fprintf(fp, "ID,No. Neighbors,Neighbors,Faces,Vertices\n");
    // Output Voronoi cells with neighbor info. 
    // See http://math.lbl.gov/voro++/doc/voro++_overview.pdf, 
    // Section 6, to interpret format string. 
    con.print_custom("%i,%s,%n,%t,%P", fp);
    fclose(fp);
}