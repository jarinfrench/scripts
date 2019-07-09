#! /usr/bin/env python3

from sys import argv, exit
import numpy as np
import argparse, random, datetime

def range_type(astr, _min = 0, _max = 1):
    value = float(astr)
    if _min <= value <= _max:
        return value
    else:
        raise argparse.ArgumentTypeError('Value not in range {}-{}'.format(_min, _max))


parser = argparse.ArgumentParser(usage = '%(prog)s [-h] num_grains',
    description = "Script that generates the txt file for use in anisotropic grain growth in MOOSE")
parser.add_argument('num_grains', type = int, help = "The number of grains that will be simulated")
parser.add_argument('-e','--iso_energy', default = 1.0, type = float, help = "The isotropic grain boundary energy value (J/m^2)")
parser.add_argument('-m','--iso_mobility', default = 1.0e-6, type = float, help = "The isotropic grain boundary mobility value (m^4/Js)")
parser.add_argument('-q','--iso_q', default = 1.0, type = float, help = "The isotropic grain boundary mobility activation energy value (eV)")
parser.add_argument('-s','--seed', default = 0, type = int, help = "The random number generator seed")
parser.add_argument('-f','--fraction', metavar = '0<=f<=1', default = 0, type = range_type, help = "The fraction of boundaries that will have a different value given")
parser.add_argument('-a','--aniso1', metavar = ('[\'e\',\'m\',\'q\']','value'), nargs = 2, help = "The (1) property that will be made anisotropic")
parser.add_argument('--aniso2', metavar = ('[\'em\',\'eq\',\'mq\']','value','value'), nargs = 3, help = "The (2) properties that will be made anisotropic")
parser.add_argument('--aniso3', metavar = 'value', nargs = 3, help = "All three properties will be made anisotropic")
parser.add_argument('-o','--output', default = '<num_grains>_grain.txt', help = "The name of the output file")
parser.add_argument('--only-energy', action = 'store_true', help = 'Flag to write only the energies to the file \'energies.txt\' (incompatible with --only-mobility and --only-q)')
parser.add_argument('--only-mobility', action = 'store_true', help = 'Flag to write only the mobilities to the file \'mobilities.txt\' (incompatible with --only-energy and --only-q)')
parser.add_argument('--only-q', action = 'store_true', help = 'Flag to write only the activatio energies to the file \'activation_energies.txt\' (incompatible with --only-energy and --only-mobility)')
parser.add_argument('-l','--logfile', default = 'aniso_grains.log', help = "The name of the log file containing the list of anisotropic boundaries")

args = parser.parse_args()

if (args.only_energy and args.only_mobility) or (args.only_energy and args.only_q) or (args.only_mobility and args.only_q):
    print("Error: Only one of --only-energy, --only-mobility, and --only-q can be specified.")
    exit(2)

if args.output == '<num_grains>_grains.txt':
    args.output = '{}_grains.txt'.format(args.num_grains)

date = datetime.datetime.now().strftime("%d %B %Y %H:%M")

energies = np.full((args.num_grains, args.num_grains), args.iso_energy)
mobilities = np.full((args.num_grains, args.num_grains), args.iso_mobility)
activation_energies = np.full((args.num_grains, args.num_grains), args.iso_q)

aniso_string = '<None>'
num_changed = 0
if args.fraction > 0:
    aniso_string = ''
    num_changed = int(args.num_grains * args.fraction)
    if num_changed < 1:
        print("No values to change - not enough grains (num_changed = {}).".format(num_changed))
        exit(1)
    elif abs(args.fraction - 1.0) < 1e-8: #100 percent change
        if args.aniso3 is not None:
            energies = np.full((args.num_grains, args.num_grains), args.aniso3[0])
            mobilities = np.full((args.num_grains, args.num_grains), args.aniso3[1])
            activation_energies = np.full((args.num_grains, args.num_grains), args.aniso3[2])
        elif args.aniso2 is not None:
            if 'e' in args.aniso2[0]:
                energies = np.full((args.num_grains, args.num_grains), args.aniso2[1])
            if 'm' in args.aniso2[0]:
                if 'm' == args.aniso2[0][1]:
                    mobilities = np.full((args.num_grains, args.num_grains), args.aniso2[2])
                else: #'m' == args.aniso2[0][0] --> 'mq'
                    mobilities = np.full((args.num_grains, args.num_grains), args.aniso2[1])
            if 'q' in args.aniso2[0]:
                activation_energies = np.full((args.num_grains, args.num_grains), args.aniso2[2])
        elif args.aniso1 is not None:
            if 'e' == args.aniso1[0]:
                energies = np.full((args.num_grains, args.num_grains), args.aniso1[1])
            elif 'm' == args.aniso1[0]:
                mobilities = np.full((args.num_grains, args.num_grains), args.aniso1[1])
            else: # 'q' == args.aniso1[0]
                activation_energies = np.full((args.num_grains, args.num_grains), args.aniso1[1])
    else:
        # generates the indices for a num_grains by num_grains upper triangular matrix (shifted right by 1)
        # indices[0] will refer to the row index, indices[1] will refer to the column index
        random.seed(args.seed) # seed the random number generator for reproducibility
        indices = np.triu_indices(args.num_grains,1)
        index_tuples = [(i,j) for i,j in zip(indices[0], indices[1])] # generates the tuples
        indices_to_change = random.sample(index_tuples, num_changed)

        flog = open(args.logfile, 'a')
        flog.write('Outfile: {output}\nDate: {date}'.format(output = args.output, date = date))
        flog.write('Anisotropic boundary properties between the following grains:\n')
        for i,j in indices_to_change:
            flog.write('gr{} gr{}\n'.format(i,j))
        flog.close()
        indices_to_change += [(j,i) for i,j in indices_to_change] # mirror the changes across the diagonal
        array_indices_to_change = (np.array([i for i,j in indices_to_change]),
                                   np.array([j for i,j in indices_to_change]))

        if args.aniso3 is not None:
            energies[array_indices_to_change] = args.aniso3[0]
            mobilities[array_indices_to_change] = args.aniso3[1]
            activation_energies[array_indices_to_change] = args.aniso3[2]
            aniso_string += 'energy, mobility, and activation energy'
        elif args.aniso2 is not None: # i.e. we have ['eq', <value_for_e>, <value_for_q>]
            if 'e' in args.aniso2[0]: # always 'em' or 'eq', so e is specified first
                energies[array_indices_to_change] = args.aniso2[1]
                aniso_string += 'energy'
            if 'm' in args.aniso2[0]:
                if 'm' == args.aniso2[0][1]: # for the case of 'em'
                    mobilities[array_indices_to_change] = args.aniso2[2]
                    aniso_string += ' and mobility'
                else: # 'm' == args.aniso2[0][0]: # for the case of 'mq'
                    mobilities[array_indices_to_change] = args.aniso2[1]
                    aniso_string += 'mobility'
            if 'q' in args.aniso2[0]: # always 'eq' or 'mq', so q is specified last
                activation_energies[array_indices_to_change] = args.aniso2[2]
                aniso_string += ' and activation energy'
        elif args.aniso1 is not None: # i.e. we have ['e', <value>]
            if 'e' == args.aniso1[0]:
                energies[array_indices_to_change] = args.aniso1[1]
                aniso_string = 'energy'
            elif 'q' == args.aniso1[0]:
                activation_energies[array_indices_to_change] = args.aniso1[1]
                aniso_string = 'activation energy'
            else: # if 'm' == args.aniso1[0]:
                mobilities[array_indices_to_change] = args.aniso1[1]
                aniso_string = 'mobility'
        else:
            print("Error: specifying --fraction requires additionally specifying one of either --aniso1, --aniso2, or --aniso3.")
            exit(2)
values_array = np.vstack((energies,mobilities,activation_energies))

if args.only_energy or args.only_mobility or args.only_q:
    if args.only_energy:
        fout = open('energies.txt', 'w')
        data = energies
    elif args.only_mobility:
        fout = open('mobilities.txt', 'w')
        data = mobilities
    elif args.only_q:
        fout = open('activation_energies.txt', 'w')
        data = activation_energies
    for i in range(len(data)):
        for j in range(len(data[i])):
            fout.write('{} '.format(data[i][j]))
        fout.write('\n')
    fout.close()
    exit(0)

fout = open(args.output, 'w')

fout.write('Generated using generate_polycrystal_data.py with N = {num_grains}, iso_energy = {iso_energy}, '
           'iso_mobility = {iso_mobility}, iso_q = {iso_q}\n'.format(num_grains = args.num_grains,
           iso_energy = args.iso_energy, iso_mobility = args.iso_mobility, iso_q = args.iso_q))
fout.write('seed = {seed}, fraction = {frac} ({num_changed} values changed), and generating anisotropic values in {aniso}\n'
           .format(seed = args.seed, frac = args.fraction, num_changed = num_changed, aniso = aniso_string))

for i in range(len(values_array)):
    for j in range(len(values_array[i])):
        fout.write('{} '.format(values_array[i][j]))
    fout.write('\n')
fout.close()
