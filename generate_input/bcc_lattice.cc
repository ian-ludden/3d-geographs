// Create Voronoi cells for a body-centered cubic (bcc) lattice in a bounding cube
#include "voro++.hh"
#include <chrono>
#include <iostream>
#include <sstream>
#include <string>
using namespace voro;

// Define bounding box
const int x_min = 0;
const int y_min = 0; 
const int z_min = 0; 

// Define parameters for container blocks (used to pre-allocate memory for container object)
const int blocks_x = 1, blocks_y = 1, blocks_z = 1;

int main(int argc, char *argv[]) {
    int index;
    int i, j, k;
    int x_max, y_max, z_max;
    int box_vol;
    double x, y, z;
    bool interactive = false;

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << "--sidelength S --output[optional] OUTFILENAME --interactive[optional]" << std::endl;
        std::cerr << "Alias for --sidelength is -s, for --output is -o, --interactive is -i.";
        return 1;
    }
    
    int side_length = -1;
    std::string out_fname;
    for (i = 1; i < argc; ++i) {
        if ("--sidelength" == std::string(argv[i]) 
            || "-s" == std::string(argv[i])) {
                if (i + 1 < argc) { side_length = atoi(argv[++i]); }
                else { std::cerr << "--sidelength (-s) option requires one argument." << std::endl; }
        }
        if ("--output" == std::string(argv[i]) 
            || "-o" == std::string(argv[i])) {
                if (i + 1 < argc) { out_fname = argv[++i]; }
                else { std::cerr << "--output (-o) option requires one argument." << std::endl; }
        }
        if ("--interactive" == std::string(argv[i]) 
            || "-i" == std::string(argv[i])) {
            interactive = true;
        }
    }

    if (side_length < 1) {
        std::cerr << "--sidelength is required and must be a positive integer.";
        return 1;
    }

    if (out_fname.empty()) {
        std::ostringstream oss;
        oss << "bch_" << side_length << ".csv";
        out_fname = oss.str();
    }

    // Set Voro++ parameters
    x_max = side_length;
    y_max = side_length;
    z_max = side_length;
    box_vol = (x_max - x_min) * (y_max - y_min) * (z_max - z_min);

    // Create a Voro++ container object, non-periodic in each dimension, 
    // with the specified dimensions and block sizes. 
    // Pre-allocate space for avg_count particles per block. 
    container con(x_min, x_max, y_min, y_max, z_min, z_max, blocks_x, 
            blocks_y, blocks_z, false, false, false, box_vol);

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
	con.draw_particles_pov("bcc_lattice_p.pov");

	// Save Voronoi cells in POV-Ray format (only feasible for small instances)
	con.draw_cells_pov("bcc_lattice_v.pov");

    // Print column headers of CSV output file
    // const char * output_fname = "bcc_lattice_5x5x5.csv";
    std::cout << "Writing Voronoi cells to " << out_fname << ".\n";
    auto start = std::chrono::high_resolution_clock::now();
    FILE * fp = safe_fopen(out_fname.c_str(), "w");
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

    if (interactive) {
        std::string response;
        std::cout << "Enter any string to exit: ";
        std::cin >> response;
        std::cout << response;
    }
}