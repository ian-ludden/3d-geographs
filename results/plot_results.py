import csv
from math import isnan

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

grid_results_fname = 'local_search_grid.out'
bcc_results_fname = 'local_search_bcc.out'
NUM_BOX_SIZES = 5
BOX_SIZES = [10, 15, 20, 25, 30]
CUBE_SIZES = [box_size**3 for box_size in BOX_SIZES]
GEOGRAPH_STATUS_TO_INDEX = {'fail_1': 0, 'fail_2': 1, 'fail_3': 2, 'fail_4': 3, 'fail_5': 4, 'success': 5}
BFS_STATUS_TO_INDEX = {'failure': 0, 'success': 1}

# Plotting parameters
MARKER_SIZE = 12

"""
Parses run_all_local_search.sh output file, 
either for grid (uniform cubic honeycomb) or 
BCC (body-centered cubic lattice) results, 
to extract data for geograph3d and 
breadth-first search (BFS) flip experiments. 

Parameters:
===============
    fname - name of the output file to be parsed

Returns:
===============
    geograph_microseconds   - NumPy array of geograph3d flip times, in microseconds, 
                              for different bounding box sizes and flip statuses
    geograph_samples        - NumPy array of geograph3d sample counts, 
                              for different bounding box sizes and flip statuses
    bfs_microseconds        - NumPy array of BFS flip times, in microseconds, 
                              for different bounding box sizes and flip statuses
    bfs_samples             - NumPy array of BFS sample counts, 
                              for different bounding box sizes and flip statuses
"""
def parse_local_search_output_file(fname):
    geograph_microseconds = np.zeros((NUM_BOX_SIZES, 6))
    geograph_samples = np.zeros((NUM_BOX_SIZES, 6))
    bfs_microseconds = np.zeros((NUM_BOX_SIZES, 2))
    bfs_samples = np.zeros((NUM_BOX_SIZES, 2))

    box_size_index = -1
    with open(fname, 'r') as f:
        f_reader = csv.reader(f)
        current_box_size = 0
        for row in f_reader:
            if not row or len(row[0]) < 4: continue
            if row[0][:4] == '../g':
                current_box_size = row[0].split('x')[1]
                print(current_box_size)
                box_size_index += 1
                continue
            
            if "With result" in row[0]: # Geograph results
                avg_microsec = float(row[1])
                if isnan(avg_microsec): avg_microsec = 0
                num_samples = int(row[3])
                print(avg_microsec, num_samples)
                status_index = GEOGRAPH_STATUS_TO_INDEX[row[0].split(' ')[-1]]
                geograph_microseconds[box_size_index, status_index] = avg_microsec
                geograph_samples[box_size_index, status_index] = num_samples
                continue
            
            if "With " in row[0]: # BFS results
                avg_microsec = float(row[1])
                if isnan(avg_microsec): avg_microsec = 0
                num_samples = int(row[3])
                print(avg_microsec, num_samples)
                status_index = BFS_STATUS_TO_INDEX[row[0].split(' ')[-1]]
                bfs_microseconds[box_size_index, status_index] = avg_microsec
                bfs_samples[box_size_index, status_index] = num_samples
                continue

    return geograph_microseconds, geograph_samples, bfs_microseconds, bfs_samples


if __name__ == '__main__':
    grid_geograph_us, grid_geograph_samples, grid_bfs_us, grid_bfs_samples = \
        parse_local_search_output_file(grid_results_fname)

    bcc_geograph_us, bcc_geograph_samples, bcc_bfs_us, bcc_bfs_samples = \
        parse_local_search_output_file(bcc_results_fname)
    
    print(grid_geograph_samples)
    print(grid_bfs_us)
    print(bcc_geograph_us)

    point_markers = ["s", "P", "^", "x", "+", "D", ">", "*", "v", "o"]

    plt.rcParams['font.size'] = '18'
    fig, ax = plt.subplots(figsize=(12, 6))
    for i in range(6):
        if i == 1: continue # Skip 3DGG fail_2, since never happened
        ax.plot(CUBE_SIZES, grid_geograph_us[:, i], '--', marker=point_markers[i], markersize=MARKER_SIZE)
    
    for i in range(2):
        ax.plot(CUBE_SIZES, grid_bfs_us[:, i], '--', marker=point_markers[i + 6], markersize=MARKER_SIZE)
    
    plt.legend(['3DGG fail_1', '3DGG fail_3', '3DGG fail_4', '3DGG fail_5', '3DGG success', 'BFS failure', 'BFS success'], bbox_to_anchor=(1.05, 1))
    plt.tight_layout()
    plt.show()