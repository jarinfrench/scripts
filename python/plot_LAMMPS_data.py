#! /usr/bin/env python3

from sys import argv, exit
from collections import defaultdict
from tabulate import tabulate
import matplotlib.pyplot as plt
import numpy as np
import csv, argparse
from myModules import *

# counts the number of times each item appears in a list
def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return ((key,locs) for key,locs in tally.items() if len(locs) > 1)

# simple way to return either the length of a list, or the "length" of a single number
def depth(data):
    try:
        return len(data)
    except:
        return 1

def determine_label(label_set, label_dict, unit_dict):
    l = []
    for i in range(depth(label_set)):
        if label_set[i] in label_dict:
            l.append(unit_dict.get(label_dict.get(label_set[i])))
    if not l or l.count(l[0]) == len(l):
        if depth(l) == 1:
            return l
        else:
            return l[0]
    else:
        return "Unknown Units"

def print_data_summary(data_set_counter, data_all, labels_all):
    print("There are {} data sets:".format(data_set_counter))
    print(tabulate([[i+1, len(data_all[i][0]), len(labels_all[i])] for i in range(data_set_counter)],
                    headers = ["Data set number", "Number of Entries", "Number of Labels"],
                    numalign = 'left'))

def parsePlottedInput(string):
    # we need to make sure that we only have brackets [], commas ',' and numbers in our string
    string = ''.join(str(string).replace(' ','')) # remove all spaces
    remaining = [i for i in list(string) if i not in ['[',']', ','] and not i.isdigit()] # anything left after removing brackets, commas, and numbers
    if all(i == ':' for i in remaining): # we handle slicing notation here
        colon_indices = [i for i,j in enumerate(list(string)) if j == ':'] # identify the indices of each colon
        for j in colon_indices: # for each index
            left_number_start_index = max([string.rfind('[', 0, j), string.rfind(',', 0, j)]) # finds the occurrence of [ or , before the number previous to the colon
            tmp = [string.find(']', j), string.find(',', j)] # finds the occurrence of , or ] after the number after the colon
            right_number_end_index = min(tmp) if min(tmp) != -1 else max(tmp) # handles the case where min returns an error code of -1
            first_num = int(string[left_number_start_index + 1 : j]) # the first number is everything from the left number start index to the colon (not inclusive)
            end_num = int(string[j + 1 : right_number_end_index]) # similarly for the end number
            tmp = [k for k in range(first_num, end_num)] # create a list
            tmp = str(tmp)[1:-1].replace(' ','') # make it a string, and remove any extra spaces)
            string = string.replace(string[left_number_start_index + 1:right_number_end_index], tmp) # replace the slice notation in the original with the list
            return eval(string) # we return a nested list of indices to plot
    if not remaining:
        return eval(string)

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] parsed_output.txt', description = "Script that extracts datasets from a parsed LAMMPS output file (using the parse_lammps_output script) and plots user-specified subsets.")
parser.add_argument('file', metavar = 'parsed_output', help = 'The parsed output file')
parser.add_argument('--one-dataset', action = 'store_true', help = "Combine all subsets into one dataset")
parser.add_argument('--no-combination', action = 'store_true', help = "Do not combine any datasets (default: prompt for combinations)")
parser.add_argument('--no-average', action = 'store_true', help = "Do not average any datasets (default: prompt for averaging)")
parser.add_argument('--print-used-data', action = 'store_true', help = 'Option to print the data used in plotting in individual text files')
args = parser.parse_args()

if args.print_used_data:
    output_filename = "plotted_data.txt"
fopen = open(args.file, 'r')
reader = csv.reader(fopen)

# Now get the information
# Note that the first 4 lines are expected to give the LAMMPS version, the box size, the number of atoms, and the unit style
# in that order.
# LAMMPS version
line = next(reader)
if line[0] != "LAMMPS version:":
    print("Unable to determine LAMMPS version.")
    version = "Unknown"
else:
    version = line[1]

# Get the box size
line = next(reader)
if line[0] != "Lx Ly Lz =":
    print("Unable to determine box size")
    Lx, Ly, Lz = [-1,-1,-1]
else:
    Lx, Ly, Lz = line[1:4]

# Number of atoms
line = next(reader)
if line[0] != "N =":
    print("Unable to determine the number of atoms.")
    N = -1
else:
    N = int(line[1])

# Unit style
line = next(reader)
if line[0] != "Unit style:":
    print("Unable to determine unit style")
    unit_style = "Unknown"
else:
    unit_style = line[1]

# Create the relevant unit labels dictionary
if unit_style == "metal":
    unit_labels = {"mass": "g/mol", "distance": r"$\AA$", "time": "ps",
                   "energy": "eV", "velocity": r"$\AA$/ps", "force": r"eV/$\AA$",
                   "torque": "eV", "temperature": "K", "pressure": "bars",
                   "dynamic viscosity": "Poise", "charge": r"Multiples of e$^-1$ charge",
                   "dipole": r"charge*$\AA$", "electric field": r"V/$\AA$",
                   "volume": r"$\AA^3$",
                   "density": r"g/cm$^3$", "none": " ", "cpu": "s"}
elif unit_style == "si":
    unit_labels = {"mass": "kg", "distance": "m", "time": "s",
                   "energy":"J", "velocity": "m/s", "force": "N",
                   "torque": "N*m", "temperature": "K", "pressure": r"N/m$^2$",
                   "dynamic viscosity": r"N*s/m$^2$", "charge": "C",
                   "dipole": "C*m", "electric field": "V/m",
                   "volume": r"m$^3$",
                   "density": r"kg/m$^3$", "none": " ", "cpu": "s"}
else:
    print("Please see the LAMMPS manual for relevant units.")
    unit_labels = "Not implemented yet"

# Dictionary from labels to units
lammps_thermo = {"Step": "none", "Elapsed": "time", "Elaplong": "time",
                 "Dt": "time", "Time": "time", "CPU": "cpu", "Atoms": "none",
                 "Temp": "temperature", "Press": "pressure", "PotEng": "energy",
                 "KinEng": "energy", "TotEng": "energy", "Enthalpy": "none",
                 "Pxx": "pressure", "Pyy": "pressure", "Pzz": "pressure",
                 "Pxy": "pressure", "Pxz": "pressure", "Pyz": "pressure",
                 "Fmax": "force", "Fnorm": "force", "E_per_atom": "energy",
                 "Lx": "distance", "Ly": "distance", "Lz": "distance",
                 "Volume": "volume", "PE_per_atom": "energy",
                 "Enthalpy_per_atom": "energy", "Volume_per_atom": "volume"}

# From this point on, we have options.  The next line will either contain the time step,
# or it will contain the labels list for the 1st data set
line = next(reader)
if line[0] == "This data is (assumed) to be from a minimization":
    time_step = [None]
    labels = next(reader)
elif line[0] == "Time step:":
    time_step = [float(line[1])]
    labels = next(reader)
else:
    print("Unable to determine time step")
    time_step = [None]
    if all(i in lammps_thermo.keys() for i in line):
        labels = line
    else:
        print("Unable to determine labels")
        print("Is there something wrong with the output file?")
        exit(1)

data = [] # initialize our data set container
for i in range(len(labels)):
    data.append([])

# We need a way to store all of the data, so we create a "super container" for our
# data, and our labels
data_all = []
labels_all = []
data_set_counter = 0 # counter for the number of data sets

# Read the remainder of the file
for row in reader:
    # two possibilities here: blank line, or more data
    # a blank line indicates a separation between data sets
    if not row:
        # new data set
        data_all.append(data[:]) #store the last data set
        labels_all.append(labels[:])
        data = [] # reset the current data set
        data_set_counter += 1
        # we know that the next line will be either data (using the same labels
        # as the previous data set), or a new timestep, and/or a new label set

        row = next(reader)
        if row[0] == "Time step:": # We have a new timestep to take into account
            time_step.append(float(row[1]))
            row = next(reader)
            if all(i in lammps_thermo.keys() for i in row): # if the next line has a list of labels
                labels = row # reassign the current labels list
                for i in range(len(labels)):
                    data.append([])
            else: # we found data, now we need to process it
                for i in range(len(labels)):
                    try:
                        data.append([])
                        data[i].append(float(row[i]))
                    except:
                        print("Error appending data to set {}".format(data_set_counter + 1)) # we use the +1 because we are indexing from 0
                        print("Error occured after time step found, and labels found")

        elif all(i in lammps_thermo.keys() for i in row): # new set of labels, same timestep
            time_step.append(time_step[data_set_counter - 1])
            labels = row
        else: # the next line is assumed to have data
            time_step.append(time_step[data_set_counter - 1])
            for i in range(len(labels)):
                try:
                    data.append([])
                    data[i].append(float(row[i]))
                except:
                    print("Error appending data to set {}".format(data_set_counter + 1)) # we use the +1 because we are indexing from 0
                    print("New labels found")

    else:
        for i in range(len(labels)):
            try:
                data[i].append(float(row[i]))
            except:
                print("Error appending data to set {}".format(data_set_counter + 1)) # we use the +1 because we are indexing from 0
                print(data)

# having finished reading the file, we need to put the last dataset into our containers
data_all.append(data)
labels_all.append(labels)
data_set_counter += 1

# Now we do some extra calculations to calculate per atom values of certain properties
for i in range(data_set_counter):
    # Total energy per atom
    if "TotEng" in labels_all[i] and N > 1:
        data_all[i].append([])
        totEngIndex = labels_all[i].index("TotEng")
        for j in range(len(data_all[i][totEngIndex])):
            data_all[i][len(labels_all[i])].append(data_all[i][totEngIndex][i] / N)
        labels_all[i].append("E_per_atom")

    # Potential Energy per atom
    if "PotEng" in labels_all[i] and N > 1:
        data_all[i].append([])
        potEngIndex = labels_all[i].index("PotEng")
        for j in range(len(data_all[i][potEngIndex])):
            data_all[i][len(labels_all[i])].append(data_all[i][potEngIndex][i] / N)
        labels_all[i].append("PE_per_atom")

    # Enthalpy per atom
    if "Enthalpy" in labels_all[i] and N > 1:
        data_all[i].append([])
        enthalpyIndex = labels_all[i].index("Enthalpy")
        for j in range(len(data_all[i][enthalpyIndex])):
            data_all[i][len(labels_all[i])].append(data_all[i][enthalpyIndex][i] / N)
        labels_all[i].append("Enthalpy_per_atom")

    # Volume per atom
    if "Volume" in labels_all[i] and N > 1:
        data_all[i].append([])
        volumeIndex = labels_all[i].index("Volume")
        for j in range(len(data_all[i][volumeIndex])):
            data_all[i][len(labels_all[i])].append(data_all[i][volumeIndex][i] / N)
        labels_all[i].append("Volume_per_atom")

    # Calculate the times if not already printed
    # NOTE: This may need to be edited to properly account for the time difference
    # when the time step changes
    if "Step" in labels_all[i] and not "Time" in labels_all[i] and time_step[i] is not None:
        data_all[i].append([])
        stepIndex = labels_all[i].index("Step")
        for j in range(len(data_all[i][stepIndex])):
            data_all[i][len(labels_all[i])].append(data_all[i][stepIndex][j] * time_step[i])
        labels_all[i].append("Time")

print_data_summary(data_set_counter, data_all, labels_all)

while not args.no_average:
    average = input("Would you like to average any subsets? (y|n): ")
    if len(average) == 1 and average in ['y', 'n']:
        break
    else:
        print("Please enter y or n")
if not args.no_average:
    if average == 'y':
        sets_to_average = input("Please enter the subsets to average using a space separated list or slicing notation: ").split()
        for i in range(len(sets_to_average)): # parse the input
            if ':' in sets_to_average[i]:
                tmp = sets_to_average[i].split(':')
                tmp = [j for j in range(int(tmp[0]) - 1, int(tmp[1]))]
                sets_to_average[i:i+1] = tmp
            else:
                sets_to_average[i] = int(sets_to_average[i]) - 1
        sets_to_average.sort()
        data = []
        labels = []
        for j in range(len(data_all[sets_to_average[0]])):
            data.append([])
            labels.append([])
        avg_data_set_index = data_set_counter + 1

        for i in sets_to_average:
            if i >= data_set_counter: # Quietly ignore out of bounds input
                break
            for j in range(len(data_all[i])):
                labels[j] = (labels_all[i][j] + '_avg')
                if labels_all[i][j] in ['Step', 'Time']: # we don't average step numbers or Time
                    data[j].append(str(data_all[i][j][0]) + "--" + str(data_all[i][j][-1]))
                else:
                    data[j].append(np.mean(data_all[i][j]))

        lammps_thermo_avg = dict()
        for key in lammps_thermo.keys():
            lammps_thermo_avg[key + "_avg"] = lammps_thermo[key]
        lammps_thermo.update(lammps_thermo_avg)
        data_all.append(data)
        labels_all.append(labels)

print("Added average dataset")

while not args.no_combination:
    combine = input("Would you like to combine any subsets? (y|n): ")
    if len(combine) == 1 and combine in ['y', 'n']:
        break
    else:
        print("Please enter y or n)")
if not args.no_combination:
    if combine == 'y':
        sets_to_combine = input("Please enter the subsets to combine using a space separated list or slicing notation: ").split()
        for i in range(len(sets_to_combine)):
            if ':' in sets_to_combine[i]:
                tmp = sets_to_combine[i].split(':')
                tmp = [j for j in range(int(tmp[0]) - 1, int(tmp[1]))]
                sets_to_combine[i:i+1] = tmp
            else:
                sets_to_combine[i] = int(sets_to_combine[i]) - 1
        sets_to_combine.sort()
        if avg_data_set_index in sets_to_combine:
            print("Removing averaged data set from combinining command")
            sets_to_combine.remove(avg_data_set_index)

        # sanity check: make sure we don't go past the maximum index
        if max(sets_to_combine) > data_set_counter:
            print("Removing unaffiliated indices")
            sets_to_combine = [i for i in sets_to_combine if i <= data_set_counter]
        for i in sets_to_combine[1:]:
            for j in range(len(data_all[i])):
                data_all[sets_to_combine[0]][j] += data_all[i][j]
            data_set_counter -= 1

    # this removes duplicate step results from each subset of data
    for i in range(data_set_counter): # for each subset of data
        for dup in sorted(list_duplicates(data_all[i][0]), reverse = True):
            for j in range(len(labels_all[i])):
                del data_all[i][j][dup[1][1]]

plot_again = True
num_new_plots = 0
while plot_again:
    print("Please enter the dataset(s) you would like to plot.")
    print("All numbers must be within the range 1-{}".format(data_set_counter))
    if not args.no_average:
        print("The averaged data set has index number {}".format(avg_data_set_index))
    data_sets_to_plot = [-1];
    while not set(data_sets_to_plot).issubset(list(range(1,data_set_counter+1)) + [avg_data_set_index]):
        data_sets_to_plot = input("Enter as a space separated list of integers, or using slicing (i.e. \'1 2 3:10\'): ").split()

        for i in range(len(data_sets_to_plot)):
            if ':' in data_sets_to_plot[i]:
                tmp = data_sets_to_plot[i].split(':')
                tmp = [j for j in range(int(tmp[0]), int(tmp[1]) + 1)]
                data_sets_to_plot[i:i+1] = tmp
            else:
                data_sets_to_plot[i] = int(data_sets_to_plot[i])
        data_sets_to_plot.sort()

    # The total number of options is going to be the sum of each of the number of
    # labels of each set specified
    max_num_labels = max(len(labels_all[i - 1]) for i in data_sets_to_plot)
    num_options = sum(len(labels_all[i - 1]) for i in data_sets_to_plot)
    option_num = 1

    # create the table
    table = []
    data_to_use = []
    labels_to_use = []
    for i in range(len(data_sets_to_plot)):
        table.append([])
        tmp = []
        for j in range(len(labels_all[data_sets_to_plot[i] - 1])):
            curr_label = labels_all[data_sets_to_plot[i] - 1][j]
            if curr_label in ['Step_avg', 'Time_avg']:
                tmp.append("() - {label}".format(option_num = option_num, label = curr_label))
            else:
                tmp.append("{option_num} - {label}".format(option_num = option_num, label = curr_label))
                data_to_use.append(data_all[data_sets_to_plot[i] - 1][j])
                labels_to_use.append(curr_label)
                option_num += 1
        if len(tmp) < max_num_labels:
            for j in range(max_num_labels - len(tmp)):
                tmp.append("--------")
        table[i] = tmp

    if len(data_sets_to_plot) > 5:
        print(bcolors.WARNING + "Warning: unable to show table of index options" + bcolors.ENDC)
        print(bcolors.WARNING + "To determine what indices to use, identify the number of labels in each data set starting\nwith 1, and then add one for each label in each consecutive set." + bcolors.ENDC)
    else:
        # Now transpose everything so it's labeled correctly
        transposed_list = list(list(zip(*table)))

        print(tabulate(transposed_list, \
            headers = ["Data set {}".format(i) for i in data_sets_to_plot]))

    print("Please enter the numbers of the datasets you would like to plot")
    print("Note: to plot multiple x or y values, simply include them in the a nested list in the desired location")
    print("i.e. to plot data sets 2, 3, and 4 against 1 on the same plot, you would type [[1,[2,3,4]]] or [[[2,3,4],1]]")
    print("Slicing is allowed")
    plotted = eval(input("Use tuple format (ex: [[x_data,y_data], [x_data2,y_data2]]): "))

    plotted = parsePlottedInput(plotted)

    data_len = []
    for i in range(len(plotted)): # find the lengths of each data set being plotted
        data_len.append([i, [depth(plotted[i][0]), depth(plotted[i][1])]])

    for i in range(len(plotted)):
        try:
            x = [j - 1 for j in plotted[i][0]] # index of x data and label
        except:
            x = plotted[i][0] - 1

        try:
            y = [j - 1 for j in plotted[i][1]] # index of y data and label
        except:
            y = plotted[i][1] - 1

        if not depth(y) == 1: # if there are multiple y axis data
            title_label_main = "["
            for j in range(len(y)):
                title_label_main += labels_to_use[y[j]]
                if not j + 1 == len(y):
                    title_label_main += ", "
                else:
                    title_label_main += "] vs "
        else:
            title_label_main = labels_to_use[y] + " vs "

        if not depth(x) == 1: # same thing with the x axis
            title_label_main += "["
            for j in range(len(x)):
                title_label_main += labels_to_use[x[j]]
                if not j + 1 == len(x):
                    title_label_main ++ ", "
                else:
                    title_label_main += "]"
        else:
            title_label_main += labels_to_use[x]

        title_label_right = "LAMMPS version: " + version + "; N = " + str(N)

        if depth(x) > 1: # determine the units for the x axis
            labels_x = []
            for j in range(len(x)):
                labels_x.append(labels_to_use[x[j]])
            x_label = determine_label(labels_x, lammps_thermo, unit_labels)
        else:
            x_label = labels_to_use[x] + " (" + determine_label([labels_to_use[x]], lammps_thermo, unit_labels)[0] + ")"

        if depth(y) > 1:
            labels_y = []
            for j in range(len(y)):
                labels_y.append(labels_to_use[y[j]])
            y_label = determine_label(labels_y, lammps_thermo, unit_labels)
        else:
            y_label = labels_to_use[y] + " (" + determine_label([labels_to_use[y]], lammps_thermo, unit_labels)[0] + ")"

        plt.figure(i + 1)

        # now we actually plot the data
        if depth(x) > 1 and depth(y) > 1:
            for j in range(depth(x)):
                for m in range(depth(y)):
                    plt.plot(data_to_use[x[j]], data_to_use[y[m]], label = labels_to_use[y[m]])
        elif depth(x) > 1 and depth(y) == 1:
            for j in range(depth(x)):
                plt.plot(data_to_use[x[j]], data_to_use[y], label = labels_to_use[y])
        elif depth(x) == 1 and depth(y) > 1:
            for j in range(depth(y)):
                plt.plot(data_to_use[x], data_to_use[y[j]], label = labels_to_use[y[j]])
        else:
            plt.plot(data_to_use[x], data_to_use[y], label=labels_to_use[y])

        # label the plot
        plt.title(title_label_main + "\n")
        plt.suptitle("\n\n" + title_label_right, fontsize=8) # makes a sub title
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        if depth(x) > 1 or depth(y) > 1: # add a legend if multiple x/y values
            plt.legend(loc='best')

    if args.print_used_data:
        unique_indices = list(set([item for sublist in plotted for item in sublist]))
        unique_indices.sort()
        fopen = open(output_filename, 'w')
        for i in unique_indices:
            fopen.write(labels_to_use[i - 1] + ' ')
        fopen.write('\n')
        for i in range(len(data_to_use[0])):
            for j in unique_indices:
                try:
                    fopen.write(repr(data_to_use[j - 1][i]) + ' ')
                except:
                    continue
            fopen.write('\n')
    plt.show()

    plot_again = promptForContinue()
    if plot_again:
        num_new_plots += 1
        output_filename = "plotted_data_{}.txt".format(num_new_plots)
