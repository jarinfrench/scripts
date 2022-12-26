#! /usr/bin/env python3

import numpy as np

def getSubStr(mystr, char1, char2):
    return mystr[mystr.find(char1) + 1 : mystr.find(char2)]

Ts = [i for i in range(1050, 1401, 50)]

for T in Ts:
    filenames = [f"T{T}/large_r/dir_final_{i}/diffusion_results.txt" for i in [1, 2, 3]]
    diffusions = {}

    for file in filenames:
        with open(file, 'r') as f:
            for line in f:
                l = line.split()
                if not line.split():
                    continue
                if l[0] == "MSD":
                    key1 = getSubStr(l[2], '_', '.')
                if l[0] == "GB" or l[0] == "Bulk":
                    key2 = l[0] # GB or Bulk
                    key3 = l[3] # x or y or z or xy or xz or yz or xyz
                    val = float(l[5])
                    if key1 in diffusions.keys():
                        if key2 in diffusions[key1].keys():
                            if key3 in diffusions[key1][key2].keys():
                                diffusions[key1][key2][key3].append(val)
                            else:
                                diffusions[key1][key2][key3] = [val]
                        else:
                            diffusions[key1][key2] = {}
                    else:
                        diffusions[key1] = {}

    for atom in diffusions.keys():
        for diff_type in diffusions[atom].keys():
            for direction in diffusions[atom][diff_type].keys():
                avg = np.mean(diffusions[atom][diff_type][direction])
                count = len(diffusions[atom][diff_type][direction])

                fout = open(f"{atom}_{diff_type}_{direction}_diffusion.txt", 'a')
                if T == 1050: # will not write this line if T1050 doesn't exist for whatever reason
                    fout.write("#T avg_diffusion (n_vals)\n")
                fout.write(f"{T} {avg} ({count})\n")
