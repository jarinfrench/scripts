#! /usr/bin/env python3

from sys import argv
import argparse
import statistics
import scipy.stats as sp
import numpy as np
from tabulate import tabulate
import natsort

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] file [file file file ...]', description = "Calculate statistics for a set of data.")
parser.add_argument('file', nargs = '+', help = "The file(s) to process")
parser.add_argument('-l','--line-start', type = int, default = 1, help = "The line number to start data analysis on")
parser.add_argument('-c','--column', type = int, help = "The column number (starting from 1) of the data to process.")
args = parser.parse_args()

args.file = natsort.natsorted(args.file)

for file in args.file:
    data = []
    print(file)
    with open(file) as f:
        for _ in range(args.line_start - 1):
            next(f);
        for line in f:
            if line.startswith("#") or line[0].isalpha():
                print("Possible header: use to specify which data set to process.")
                print(line)
            else:
                sub_data = []
                for i in line.split():
                    sub_data.append(float(i))
                data.append(sub_data)

    data = np.array(data)
    if not args.column:
        while True:
            dataset_to_process = int(input("Please specify the number (1 - {}) of the dataset to process: ".format(len(data[0]))))
            if dataset_to_process > len(data[0]) or dataset_to_process < 1:
                continue
            else:
                break
    else:
        dataset_to_process = args.column

    data_to_process = data[:,dataset_to_process - 1]

    try:
        geometric_mean = sp.gmean(data_to_process)
    except {ValueError, RuntimeWarning}:
        print("Unable to calculate geometric mean")
        geometric_mean = 0

    try:
        harmonic_mean = sp.hmean(data_to_process)
    except ValueError:
        print("Unable to calculate harmonic mean")
        harmonic_mean = 0

    arithmetic_mean = statistics.mean(data_to_process)
    median = statistics.mean(data_to_process)
    maximum = max(data_to_process)
    minimum = min(data_to_process)
    standard_deviation = statistics.stdev(data_to_process)
    variance = statistics.variance(data_to_process)

    print(tabulate([["G_Mean", geometric_mean], ["H_Mean", harmonic_mean], \
        ["A_Mean", arithmetic_mean], ["Median", median], ["Max", maximum], \
        ["Min", minimum],["SD", standard_deviation], ["Variance", variance]], headers = \
        ["Statistic", "Value"]))
    print("\n")
# print("{0:8s} {1:8s} {2:8s} {3:8s} {4:8s} {5:8s} {6:8s}".format("G_Mean", "H_Mean", "A_mean", "Median", "Range", "SD", "Variance", width = 10))
# print("{0:6.2f} {1:6.2f} {2:6.2f} {3:6.2f} {4:6.2f} {5:6.2f} {6:6.2f}".format(geometric_mean, harmonic_mean, arithmetic_mean, median, range, standard_deviation, variance, width = 10))
