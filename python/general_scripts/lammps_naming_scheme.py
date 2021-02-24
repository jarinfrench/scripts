#! /usr/bin/env python3

import regex as re
import argparse
from colorama import Fore, Back, Style
from sys import exit
from os import rename

compound_regex = r"_(?P<compound>(?:[A-Z](?:[a-z])?(?:[1-9])*)+)"
axis_regex = r"_(?P<axis>[0-9][0-9][0-9])"
size_regex = r"_(?P<size>[1-9][0-9]*x[1-9][0-9]*x[1-9][0-9]*)"
potential_regex = r"_p(?P<potential>(?:[A-Z][a-z]*)+)"
angle_regex = r"(?:_(?P<angle>[0-9]+(?:\.[0-9]+)?)degree)"
radius_regex = r"(?:_r(?P<radius>[1-9][0-9]*)A)"
concentration_regex = r"(?:_c(?P<concentration>0\.[0-9]+)(?P<type>wt|at)%(?P<solute>[A-Z][a-z]))"
extra_regex = r"(?:_(?P<extra>[a-zA-Z0-9_.-]*))?\.(?P<extension>dat|lmp)$"
compound_re = re.compile(compound_regex)
axis_re = re.compile(axis_regex)
size_re = re.compile(size_regex)
angle_re = re.compile(angle_regex)
r_re = re.compile(radius_regex)
pot_re = re.compile(potential_regex)
u_re = re.compile(concentration_regex)
extra_re = re.compile(extra_regex)
filename_regex = r"LAMMPS" + compound_regex + axis_regex + size_regex + potential_regex + angle_regex + "?" + radius_regex + "?" + concentration_regex + "?" + extra_regex
filename_re = re.compile(filename_regex)

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] [options]', description = "Validates the format of the LAMMPS data file name given to it. Offers suggestions if it fails.")
parser.add_argument('file', help = "The file name to validate")

args = parser.parse_args()

fancy_checkmark = Fore.GREEN + Style.BRIGHT + u'\u2713' + Style.RESET_ALL
fancy_x = Fore.RED + Style.BRIGHT + u'\u2717' + Style.RESET_ALL

if __name__ == "__main__":

    if filename_re.findall(args.file):
        print(f"{fancy_checkmark} {args.file}")
        exit(0)

    print(f"{fancy_x} {args.file}")
    suggestion = []

    if len(filename_re.findall(args.file)) == 0:
        result = args.file
        suggestion += ["LAMMPS"]
        if not args.file.startswith("LAMMPS"):
            print(f"{Fore.YELLOW}Start the filename with \"LAMMPS\"{Style.RESET_ALL}")
            args.file = "LAMMPS_" + args.file
        result = result.replace("LAMMPS", '')
        # independent pieces - required
        compound_res = compound_re.findall(args.file)
        result = result.replace(''.join(compound_res), '')
        axis_res = axis_re.findall(args.file)
        result = result.replace(''.join(axis_res), '')
        size_res = size_re.findall(args.file)
        result = result.replace(''.join(size_res), '')
        pot_res = pot_re.findall(args.file)
        result = result.replace(''.join(["p"] + pot_res), '')

        # if one is here, they both should be
        angle_res = angle_re.findall(args.file)
        result = result.replace(''.join(angle_res + ["degree"]), '')
        r_res = r_re.findall(args.file)
        result = result.replace(''.join(["r"] + r_res + ["A"]), '')

        # optional stuff
        u_res = u_re.findall(args.file)
        for match in u_res:
            result = result.replace(''.join(["c"] + [match[0]] + [match[1]] + ["%"] + [match[2]]), '')
        extra_res = extra_re.search(args.file)
        result = result.replace(''.join(extra_res.group("extension")), '')
        result = result.lstrip('_')
        result = result.rstrip('.')
        result = [result.rstrip('_')]

        # test the independent pieces first
        if not compound_res:
            suggestion += ["<compound>"]
            print(f"{Fore.YELLOW}Put one compound in the filename{Style.RESET_ALL}")
        else:
            suggestion += compound_res

        if not axis_res:
            suggestion += ["<axis>"]
            print(f"{Fore.YELLOW}Put the orientation in the filename (i.e. which crystal direction is aligned with the z direction of the simulation cell?){Style.RESET_ALL}")
        else:
            suggestion += axis_res

        if not size_res:
            suggestion += ["<size>"]
            print(f"{Fore.YELLOW}Put the size of the system in the filename{Style.RESET_ALL}")
        else:
            suggestion += size_res

        if not pot_res:
            suggestion += ["<potential>"]
            print(f"{Fore.YELLOW}Put the potential in the filename{Style.RESET_ALL}")
        else:
            suggestion += ["p" + ''.join(pot_res)]

        if angle_res or r_res:
            if angle_res and r_res:
                suggestion += [''.join(angle_res) + "degree"]
                suggestion += ["r" + ''.join(r_res) + "A"]
            else:
                suggestion += ["<angle>"]
                suggestion += ["<radius>"]
                print(f"{Fore.YELLOW}Both the angle and the radius should be in the filename{Style.RESET_ALL}")

        if u_res:
            for conc, type, el in u_res:
                suggestion += ["c" + ''.join(conc) + ''.join(type) + "%" + ''.join(el)]

        if not extra_res:
            print(f"{Fore.YELLOW}No file extension found: use either dat or lmp{Style.RESET_ALL}")
        else:
            suggestion += [''.join(result) + "." + extra_res.group("extension")]

        suggestion = "_".join(suggestion)
        print(suggestion)

        if not any([c in suggestion for c in ['<', '>']]):
            change_name = input(f"{Fore.GREEN}Would you like to rename the file to the suggestion shown?{Style.RESET_ALL} (Y|n): ")
            if any([c in change_name for c in ['y', 'Y']]):
                rename(args.file, suggestion)
