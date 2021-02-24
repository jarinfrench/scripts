#! /usr/bin/env python3

import os
import math
from sys import argv, exit

cwd = os.getcwd().split("/")
shorten_individual_dirs = False if os.getenv('PROMPT_DIRTRIM') is not None and int(os.getenv('PROMPT_DIRTRIM')) > 0 else True
user = os.getenv('USER')
sysname = os.getenv('SYSNAME')
if len(argv) == 2:
    num_dirs = int(argv[1])
    if num_dirs == 0: # ignoring the shortening effect
        num_dirs = 100000 # we set the number of directories to be a number that no one has business having directories this deep

else:
    if os.getenv('PROMPT_DIRTRIM') is not None:
        num_dirs = int(os.getenv('PROMPT_DIRTRIM'))
        if num_dirs <= 3:
            num_dirs = 4
    else:
        num_dirs = 4

if len(cwd) < 3:
    result = ""
else:
    if "/".join(cwd[0:3]) == f"/home/{user}":
        result = "~"
        cwd = cwd[3:] # truncate the beginning part of the path
    elif "/".join(cwd[0:3]) == f"/media/{user}":
        result = "m~"
        if len(cwd) > 3:

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
            cwd = [""]
    elif len(cwd) > 4 and "/".join(cwd[0:4]) == f"/work/{sysname}/{user}":
        result = "w~"
        cwd = cwd[4:]
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
    if len(cwd) > num_dirs:
        cwd = cwd[0:math.floor(num_dirs / 2)] + ["..."] + cwd[-math.floor(num_dirs / 2):]
    print("/".join([result] + cwd))
