# 3d-geographs
Implementation of three-dimensional geo-graphs, a dynamic data structure for efficient contiguity verification in 3D graph partitioning. 

## Dependencies
- [Voro++](http://math.lbl.gov/voro++/), version 0.4.6. See [http://math.lbl.gov/voro++/download/](http://math.lbl.gov/voro++/download/) for download and installation instructions. Used in "Generate input files" step below. 
- (Optional) [Jupyter](https://jupyter.org/) for "Analyze experimental results" step below. The Jupyter Notebook includes several Python module dependencies. 

## Replication steps
(WIP) The following steps replicate the experiments in (paper link when available), starting from a clean clone of this repository. 

### Generate input files
The `generate_input/` subfolder contains C++ programs for creating cubic honeycomb (CH) and bitruncated cubic honeycomb (BCH) Voronoi tesselations using Voro++. Commands to compile and run these programs:
```
cd generate_input
make 
./generate_all_inputs.sh
```
The script will create `grid` and `bcc` folders in `generate_input/` containing the CH and BCH instances, respectively. For each type of honeycomb, the bounding box side length varies from 10 to 40 in increments of 5. The entire process takes a few hours on a standard workstation. 

### Run random local search
The `src/` folder contains the 3D geo-graph source files and an implementation of random local search for testing the 3D geo-graph. To compile the 3D geo-graph and local search implementation, use 
```
cd src
make
```
To then run all the local search experiments, comparing the 3D geo-graph against breadth-first search (BFS), call `./run_all_local_search.sh`. This will create output files `results/data_[bcc or grid].csv` and `results/partition_[bcc or grid]_K[value of K].out` with flip timing, status details, and final partitions. 

### Analyze experimental results
The `results/Experiments_Analysis.ipynb` Jupyter Notebook loads the experimental results from `results/data_[bcc or grid].csv`. The notebook then generates plots and tables for publication. 

### Visualizing final partitions
The following steps produce PNG files for each final partition and the individual zones. 
1. Run `generate_input/create_all_semicolon_versions.sh` to produce semicolon-separated versions of the honeycomb input files. 
2. From the `results` directory, run `python produce_povrays_of_partitions.py`. This creates a POV file for each partition and zone from the semicolon-separated honeycomb files (from Step 1) and the final partition CSV files (from "Run random local search" above). 
3. (Requires POV-Ray installation) Run `results/run_all_povray.sh` to invoke the POV-Ray engine for each POV file. This produces a PNG file from each POV and moves both to the `results/img/` directory. 
