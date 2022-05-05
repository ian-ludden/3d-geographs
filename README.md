# 3d-geographs
Implementation of three-dimensional geo-graphs, a dynamic data structure for efficient contiguity verification in 3D graph partitioning. 

## Dependencies
- [Voro++](http://math.lbl.gov/voro++/), version 0.4.6. See [http://math.lbl.gov/voro++/download/](http://math.lbl.gov/voro++/download/) for download and installation instructions. Used in "Generate input files" step below. 
- (Optional) [Jupyter](https://jupyter.org/) for "Analyze experimental results" step below. The Jupyter Notebook includes several Python module dependencies. 

## Replication steps
The following steps 

### Generate input files
The `generate_input/` subfolder contains C++ programs for creating cubic honeycomb (CH) and bitruncated cubic honeycomb (BCH) Voronoi tesselations using Voro++. Commands to compile and run these programs:
```
cd generate_input
make 
./generate_all_inputs.sh
```
The script will create `grid` and `bcc` folders in `generate_input/` containing the CH and BCH instances, respectively. For each type of honeycomb, the bounding box side length varies from 10 to 40 in increments of 5. 

### Run random local search
The `src/` folder contains the 3D geo-graph source files and an implementation of random local search for testing the 3D geo-graph. To compile the 3D geo-graph and local search implementation, use 
```
cd src
make
```
To then run all the local search experiments, comparing the 3D geo-graph against breadth-first search (BFS), call `./run_all_local_search.sh`. This will create output files `results/local_search_grid.out` and `results/local_search_bcc.out` with flip timing and status details. 

### Analyze experimental results
The `results/Experiments_Analysis.ipynb` Jupyter Notebook loads the experimental results from `results/local_search_[bcc/grid].out` and generates plots for publication. 