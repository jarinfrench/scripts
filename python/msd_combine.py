import numpy as np
import natsort as ns
import glob
import csv

roots = ["UO2_100_20degree_Basak",
        "UO2_100_45degree_Basak",
        "UO2_110_20degree_Basak",
        "UO2_110_45degree_Basak",
        "UO2_111_20degree_Basak",
        "UO2_111_45degree_Basak",
        "UO2_111_sigma7_Basak",
        "UO2_112_45degree_Basak"
        ]
for r in roots:
    files_1 = ns.natsorted(glob.glob(f"{r}_*dir_1_MSD_U.dat"))
    files_2 = ns.natsorted(glob.glob(f"{r}_*dir_2_MSD_U.dat"))
    files_3 = ns.natsorted(glob.glob(f"{r}_*dir_3_MSD_U.dat"))

    for idx,f_list in enumerate((files_1, files_2, files_3)):
        data = [np.loadtxt(i, dtype={'names': ('Timestep', 'MSD X', 'MSD Y', 'MSD Z', 'MSD'), 'formats': ('int', 'float', 'float', 'float', 'float')}) for i in f_list]
        temps = [i.split('_')[4] for i in f_list]
        m = max([i.size for i in data])
        for jdx, _ in enumerate(data):
            diff = m - data[jdx].size
            if diff == 0:
                data[0]['Timestep'] = data[jdx]['Timestep']
                continue
            data[jdx] = np.pad(data[jdx], (0,diff), 'constant', constant_values = (-1))
        result_list = np.empty((m, 4*len(data) + 1))
        result_list[:,0] = data[0]['Timestep']
        col_idx = 1
        for a in data:
            result_list[:,col_idx] = a['MSD X']
            col_idx += 1
            result_list[:,col_idx] = a['MSD Y']
            col_idx += 1
            result_list[:,col_idx] = a['MSD Z']
            col_idx += 1
            result_list[:,col_idx] = a['MSD']
            col_idx += 1
        header = ["#", "timestep"] + ns.natsorted(temps * 4)
        with open(f"{r}_dir_{idx+1}.txt", 'w') as f:
            writer = csv.writer(f, delimiter = ' ')
            writer.writerow(header)
            for row in result_list:
                writer.writerow(row)
