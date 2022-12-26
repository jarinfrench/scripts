#! /usr/bin/env python3.9

import argparse
from natsort import natsorted

def positiveInt(num):
    try:
        num = int(num)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid int value {num}")
    if num <= 0:
        raise argparse.ArgumentTypeError(f"Invalid positive int value {num}")
    else:
        return num

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] file(s) [options]',
    description = "Calculates the center of mass of a collection of grains based on their grain ID",
    epilog = "This script reads output formatted after the find_grains executable output")
parser.add_argument('files', nargs = '+', help = "The file(s) to process")
parser.add_argument('-g', '--show-grain', nargs = '+', type = positiveInt, default = [0], help = "The grain(s) to show results for (default: [0] --> all grains)")
parser.add_argument('-y', '--types', nargs = '+', type = positiveInt, default = [0], help = "The atom types to process (default: [0] --> all types evaluated)")
parser.add_argument('-N','--grains', type = positiveInt, default = 2, help = "The number of grains in the simulation")
parser.add_argument('-B', '--bounds', nargs = 6, type = float, help = "The bounds of the simulation. If given, the center of the grain will also be shown as a percent of the system bounds. Given as xlo xhi ylo yhi zlo zhi")
parser.add_argument('-p', action = 'store_true', help = "Flag to show the location of the grain center(s) in percent only. Requires -B")
parser.add_argument('-t', action = 'store_true', help = "Flag to show space separated values of the grain center. If -B is specified, values shown are percent")

args = parser.parse_args()

if args.p and not args.bounds:
    raise argparse.ArgumentTypeError(f"The -p flag must be specified with the -B flag")

for file in natsorted(args.files):
    centers = [[0.0, 0.0, 0.0]] * args.grains
    counts = [0] * args.grains
    if not args.t:
        print(f"Processing file {file}...")
    with open(file) as f:
        var_names = f.readline().split(",")
        x_idx = [i for i,val in enumerate(var_names) if "X" in val and not "Xu" in val][0]
        y_idx = [i for i,val in enumerate(var_names) if "Y" in val and not "Yu" in val][0]
        z_idx = [i for i,val in enumerate(var_names) if "Z" in val and not "Zu" in val][0]
        type_idx = [i for i,val in enumerate(var_names) if "Atom Type" in val][0]
        grain_num_idx = [i for i,val in enumerate(var_names) if "Grain Number" in val][0]

        for line in f:
            if line.startswith("#"):
                continue
            data = line.split()[0:grain_num_idx + 1]

            type = int(data[type_idx])
            if not args.types[0] == 0 and type not in args.types:
                continue
            grain_num = int(data[grain_num_idx])
            x = float(data[x_idx])
            y = float(data[y_idx])
            z = float(data[z_idx])
            counts[grain_num - 1] += 1
            centers[grain_num - 1] = [i+j for i,j in zip(centers[grain_num - 1], [x,y,z])]

        for i, tup in enumerate(zip(centers, counts)):
            grain = tup[0]
            c = tup[1]
            if c == 0:
                if args.t:
                    print("0.5 0.5 0.5")
                else:
                    print(f"  Grain {i + 1} no longer exists")

            if not args.show_grain[0] == [0] and i + 1 not in args.show_grain:
                continue
            grain = [j / c for j in grain]
            if args.bounds:
                bounds = [args.bounds[1] - args.bounds[0],
                          args.bounds[3] - args.bounds[2],
                          args.bounds[5] - args.bounds[4]]
                grain_percent = [j/k for j,k in zip(grain,bounds)]

                if args.p:
                    if args.t:
                        print(f"{' '.join(f'{gp:.4f}' for gp in grain_percent)}")
                    else:
                        print(f"  Center of grain {i + 1} ({c} atoms): {'[' + ','.join(f'{gp:.2%}' for gp in grain_percent) + ']'}")
                else:
                    if args.t:
                        print(f"{' '.join(f'{gg:.4f}' for gg in grain)}")
                    else:
                        print(f"  Center of grain {i + 1} ({c} atoms): {'[' + ','.join(f'{gg:.4f}' for gg in grain) + ']'} ({'[' + ','.join(f'{gp:.2%}' for gp in grain_percent) + ']'})")
            else:
                if args.t:
                    print(f"{' '.join(f'{gg:.4f}' for gg in grain)}")
                else:
                    print(f"  Center of grain {i + 1} ({c} atoms): {'[' + ','.join(f'{gg:.4f}' for gg in grain) + f']'}")
