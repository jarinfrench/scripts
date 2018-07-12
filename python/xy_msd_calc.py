#! /usr/bin/env python3

import os

dir_path = "/media/jarinf/Seagate Expansion Drive/UO2/grain_growth/MSD_data/20degree/U/average_data/"
filenames = [f for f in os.listdir("/media/jarinf/Seagate Expansion Drive/UO2/grain_growth/MSD_data/20degree/U/average_data") if "with_xy" not in f] # <-- this doesn't quite work. It doesn't give me the absolute file path.
new_filenames = [dir_path + os.path.splitext(base)[0] + "_with_xy.dat" for base in filenames]

for f1,f2 in zip(filenames, new_filenames):
    line_num = 1
    timesteps = []
    msd_xys = []
    with open(dir_path + f1, 'r') as fin, open(f2, 'w') as fout:
        for line in fin:
            if line_num == 1:
                var_list = line.split()
                var_list.append("\"MSD_xy\"")
                fout.write(' '.join(var_list) + "\n")
            elif line_num == 2:
                fout.write(line)
            else:
                data = line.split()
                msd_xy = float(data[1]) + float(data[2])
                data.append(str("%8.6f"%msd_xy))
                fout.write(' '.join(data) + "\n")
            line_num += 1
