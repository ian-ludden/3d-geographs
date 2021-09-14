// Create Voronoi cells for a body-centered cubic (bcc) lattice in a bounding cube
#include "voro++.hh"
#include <chrono>
#include <iostream>
#include <string>
using namespace voro;

// Define bounding box
const int x_min = 0;
const int x_max = 20;
const int y_min = 0; 
const int y_max = 20;
const int z_min = 0; 
const int z_max = 20; 

const int box_vol = (x_max - x_min) * (y_max - y_min) * (z_max - z_min);

// Define parameters for container blocks (used to pre-allocate memory for container object)
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

    // Add particles for body-centered cubic (bcc) lattice
    index = 0;

    // Start with unit cube centers
    for (i = x_min; i < x_max; i++) {
        x = i + 0.5;
        for (j = y_min; j < y_max; j++) {
            y = j + 0.5;
            for (k = z_min; k < z_max; k++) {
                z = k + 0.5;
                con.put(index, x, y, z);
                index++;
            }
        }
    }

    // Now, add unit cube vertices, excluding outer boundary
    for (i = x_min + 1; i < x_max; i++) {
        x = i;
        for (j = y_min + 1; j < y_max; j++) {
            y = j;
            for (k = z_min + 1; k < z_max; k++) {
                z = k;
                con.put(index, x, y, z);
                index++;
            }
        }
    }

	// Save particles in POV-Ray format (only feasible for small instances)
	// con.draw_particles_pov("bcc_lattice_p.pov");

	// Save Voronoi cells in POV-Ray format (only feasible for small instances)
	// con.draw_cells_pov("bcc_lattice_v.pov");

    // Print column headers of CSV output file
    const char * output_fname = "bcc_lattice_20x20x20.csv";
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