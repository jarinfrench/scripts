#! /usr/bin/env python3

from sys import argv
import argparse
import statistics
import scipy.stats as sp
import numpy as np
from tabulate import tabulate
import natsort
import warnings

np.seterr(all='warn')
warnings.filterwarnings('error')

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] file [file file file ...]', description = "Calculate statistics for a set of data.")
parser.add_argument('file', nargs = '+', help = "The file(s) to process")
parser.add_argument('-l','--line-start', type = int, default = 1, help = "The line number to start data analysis on")
parser.add_argument('-c','--column', type = int, help = "The column number (starting from 1) of the data to process.")
parser.add_argument('--summary', action = 'store_true', help = "Only show the summary of the data")
args = parser.parse_args()

args.file = natsort.natsorted(args.file)
g_mean_summary = 0
h_mean_summary = 0
a_mean_summary = 0
q1_summary = 0
median_summary = 0
q3_summary = 0
max_summary = 0
min_summary = 0
sd_summary = 0
variance_summary = 0
range_summary = 0
mode_summary = 0


for file in args.file:
    data = []
    if not args.summary:
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
    except (ValueError, RuntimeWarning):
        print("Unable to calculate geometric mean")
        geometric_mean = 0

    try:
        harmonic_mean = sp.hmean(data_to_process)
    except (ValueError, RuntimeWarning):
        print("Unable to calculate harmonic mean")
        harmonic_mean = 0

    arithmetic_mean = statistics.mean(data_to_process)
    quartile_25, median, quartile_75 = np.percentile(data_to_process, [25,50,75])
    maximum = max(data_to_process)
    minimum = min(data_to_process)
    range_val = maximum - minimum
    standard_deviation = statistics.stdev(data_to_process)
    variance = statistics.variance(data_to_process)

    try:
        mode = statistics.mode(data_to_process)
    except:
        mode = statistics.mode([round(i,-1) for i in data_to_process])

    print("segregation energy: {}".format(mode - minimum))
    if not args.summary:
        print(tabulate([["G_Mean", geometric_mean], ["H_Mean", harmonic_mean], \
            ["A_Mean", arithmetic_mean], ["Min", minimum], ["Quartile1", quartile_25], \
            ["Median", median], ["Quartile3", quartile_75], ["Max", maximum],  ["Mode", mode], \
            ["Range", range_val], ["SD", standard_deviation], ["Variance", variance]], headers = \
            ["Statistic", "Value"]))
        print("\n")

    g_mean_summary += geometric_mean
    h_mean_summary += harmonic_mean
    a_mean_summary += arithmetic_mean
    q1_summary += quartile_25
    median_summary += median
    q3_summary += quartile_75
    max_summary += maximum
    min_summary += minimum
    range_summary += range_val
    sd_summary += standard_deviation
    variance_summary += variance
    mode_summary += mode

g_mean_summary /= len(args.file)
h_mean_summary /= len(args.file)
a_mean_summary /= len(args.file)
q1_summary /= len(args.file)
median_summary /= len(args.file)
q3_summary /= len(args.file)
max_summary /= len(args.file)
min_summary /= len(args.file)
range_summary /= len(args.file)
sd_summary /= len(args.file)
variance_summary /= len(args.file)
mode_summary /= len(args.file)

print("Statistics Summary:")
print(tabulate([["G_Mean", g_mean_summary], ["H_Mean", h_mean_summary], \
    ["A_Mean", a_mean_summary], ["Min", min_summary], ["Quartile1", q1_summary], \
    ["Median", median_summary], ["Quartile3", q3_summary], ["Max", max_summary], ["Mode", mode_summary], \
    ["Range", range_summary], ["SD", sd_summary], ["Variance", variance_summary]], headers = \
    ["Statistic", "Value"]))
print("\n")
