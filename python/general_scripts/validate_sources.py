#! /usr/bin/env python3

import os
from natsort import natsorted
from colorama import Fore, Back, Style
import argparse

#TODO: Include globbing and regex to allow for parsing the * and () in the Sources.txt

parser = argparse.ArgumentParser(usage = 'validate_sources.py', description = "Checks ${HOME}/Pictures/Research_images/Outside_sources to see if all images are cited in that directory's Sources.txt")
parser.add_argument('-n', '--not-cited', action = 'store_true', help = 'Flag to only show the uncited images')
parser.add_argument('-i', '--ignore', nargs = '+', help = "Subdirectories to ignore")
parser.add_argument('-r', '--reference', default = os.getenv('HOME')+"/Pictures/Research_images/Outside_sources", help = "The reference directory to check (default: ${HOME}/Pictures/Research_images/Outside_sources)")
parser.add_argument('-s', '--sources', default = "Sources.txt", help = "The file containing the image sources (default: Sources.txt)")
args = parser.parse_args()
if args.ignore is None:
    args.ignore = []

checkmark = Fore.GREEN + Style.BRIGHT + u'\u2713'
fancy_x = Fore.RED + Style.BRIGHT + u'\u2717'
fancy_question = Fore.YELLOW + Style.BRIGHT + u'\u2753'

os.chdir(args.reference)

files = []
for (dirpath, dirnames, filenames) in os.walk("."):
    files.extend(os.path.join(dirpath, filename) for filename in filenames)

files = ["/".join(file.split("/")[1:]) for file in files] # remove the leading ./
files.remove(args.sources) # remove Sources.txt from the list
files = natsorted(files) # sort the list
files = [i for i in files if i.split("/")[0] not in args.ignore] # Remove files in directories specified by args.ignore

source_entries = []

with open(args.sources) as f:
    for line in f:
        if (line.split()[0]).split("/")[0] not in args.ignore: # skip the same directories in the sources file
            source_entries.extend([line.split()[0]])

source_entries = natsorted(source_entries)

num_not_cited = 0
title_printed = False
for file in files:
    if file not in source_entries:
        if args.not_cited and not title_printed:
            print("Files without a source")
            title_printed = True
        print(f"{fancy_x} - {file}{Style.RESET_ALL}")
        num_not_cited += 1
    else:
        if not args.not_cited:
            print(f"{checkmark} - {file}{Style.RESET_ALL}")
print(f"{num_not_cited} images not cited\n")

num_cited_but_not_found = 0
title_printed = False
for source in source_entries:
    if source not in files:
        if not title_printed:
            print("Sources specified, but no corresponding file")
            title_printed = True
        print(f"{fancy_question} - {source}{Style.RESET_ALL}")
        num_cited_but_not_found += 1

print(f"{num_cited_but_not_found} extraneous sources")
