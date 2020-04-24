#! /usr/bin/env python3

import os
from sys import exit

cwd = os.getcwd().split("/")

shorten_individual_dirs = True if os.getenv('PS1_SHORTEN_INDIVIDUAL') is not None and int(os.getenv('PS1_SHORTEN_INDIVIDUAL')) > 0 else False

if len(cwd) < 3:
    result = ""
else:
    if "/".join(cwd[0:3]) == "/home/jarinf":
        result = "~"
        cwd = cwd[3:] # truncate the beginning part of the path
    elif "/".join(cwd[0:3]) == "/media/jarinf":
        result = "m~"
        if len(cwd[3].split()) > 1: # More than one word in this directory
            result += cwd[3].split()[0][0:3].upper()
        else: # just one word
            for i,char in enumerate(cwd[3]):
                if i >= 3:
                    break
                else:
                    result += char.upper()
        cwd = cwd[4:] # truncate the beginning part of the path
    else:
        result = ""

if shorten_individual_dirs:
    if len("/".join([result] + cwd)) < 50:
        print("/".join([result] + cwd))
        exit(0)
    for directory in cwd:
        if len(directory) > 8:
            result += "/{}...".format(directory[0:5])
        else:
            result += "/{}".format(directory)
    print(result)
else:
    if len(cwd) > 4:
        cwd = [cwd[0], cwd[1], "...", cwd[-2], cwd[-1]]
    print("/".join([result] + cwd))
