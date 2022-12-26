#! /usr/bin/env python

from ovito.io import import_file, export_file
import ovito.modifiers as om
from ovito.pipeline import FileSource
from ovito.data import DataTable
import numpy as np
import argparse

fnn = {"bcc": np.sqrt(3) * 0.5, "fcc": 1.0 / np.sqrt(2)}

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] structure a0', description = 'Performs a cluster analysis on the simulation snapshots')
# parser.add_argument('structure', choices = ['bcc', 'fcc'], help = "Base crystal structure")
# parser.add_argument('a0', type = float, help = "The lattice parameter of the system")
# parser.add_argument('-m', '--min', type = int, default = 3, help = "The minimum cluster size to report")
args = parser.parse_args()

# cutoff = fnn[args.structure] * args.a0
cutoff = 4.4 # Taken from https://www.princeton.edu/~maelabs/mae324/glos324/xenon.htm
pipeline = import_file("*.dump") # load all dump files in the current directory
n_types = len(set(pipeline.compute().particles.particle_type))
n_frames = pipeline.source.num_frames
for i in range(2, n_types + 1): # the len(...) command gets the number of unique atom types. The +1 allows the range to include that number
    pipeline.modifiers.append(om.SelectTypeModifier(types = {i})) # selects type i atoms
    pipeline.modifiers.append(om.ClusterAnalysisModifier(cutoff = cutoff, sort_by_size = True, compute_com = True, only_selected = True))

export_file(pipeline, f'cluster_results_type2_d{cutoff}_initial.txt', 'txt/table', key = 'clusters', frame = 0)
export_file(pipeline, f'cluster_results_type2_d{cutoff}_final.txt', 'txt/table', key = 'clusters', frame = n_frames)
if n_types > 2:
    for i in range(3, n_types + 1):
        export_file(pipeline, f'cluster_results_type{i}_initial.txt', 'txt/table', key = f'clusters.{i-1}', frame = 0)
        export_file(pipeline, f'cluster_results_type{i}_final.txt', 'txt/table', key = f'clusters.{i-1}', frame = n_frames)
