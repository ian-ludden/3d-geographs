import csv
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

grid_results_fname = 'local_search_grid.out'
bcc_results_fname = 'local_search_bcc.out'

if __name__ == '__main__':
    print('TODO: plot results')

    with open(grid_results_fname, 'r') as grid_f:
        grid_f_reader = csv.reader(grid_f)
        current_box_size = 0
        for row in grid_f_reader:
            if not row: continue
            if row[0][:4] == '../g':
                current_box_size = row[0].split('x')[1]
                print(current_box_size)
                continue
            
            if "With result" in row[0]:
                avg_microsec = float(row[1])
                num_samples = int(row[3])
                print(avg_microsec, num_samples)
                continue

    with open(bcc_results_fname, 'r') as grid_f:
        grid_f_reader = csv.reader(grid_f)
        current_box_size = 0
        for row in grid_f_reader:
            if not row: continue
            if row[0][:4] == '../g':
                current_box_size = row[0].split('x')[1]
                print(current_box_size)
                continue
            
            if "With result" in row[0]:
                avg_microsec = float(row[1])
                num_samples = int(row[3])
                print(avg_microsec, num_samples)
                continue


