#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
from io import StringIO
import warnings

kB = 8.617e-5 # eV/K
def TtoArr(x):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return 1 / (kB * x)

def ArrtoT(x):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return 1 / (x * kB)

files = [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith("diffusion.txt")]

ddtypes = {"names": ("T", "D", "nfiles"), "formats": ("int", "float", "int")}
all_data = {}

for f in files:
    atom_type = f.split('_')[0] # U, Mo, or Xe
    diff_type = f.split('_')[1] # Bulk or GB
    diff_direction = f.split('_')[2] # x, y, z, xy, xz, yz, or xyz
    if not atom_type in all_data.keys():
        all_data[atom_type] = {diff_type: {diff_direction: None}}
    if not diff_type in all_data[atom_type].keys():
        all_data[atom_type][diff_type] = {diff_direction: None}
    if atom_type not in ["Mo", "U", "Xe"]:
        print(f"Unrecognized atom type {atom_type}")
        exit()
    if diff_type not in ["Bulk", "GB"]:
        print(f"Unrecognized diffusion type {diff_type}")
        exit()
    # data = np.loadtxt(f, dtype = ddtypes, usecols = (0,1))
    with open(f, 'r') as tmp:
        dtmp = tmp.read().replace('(', '')
        dtmp = dtmp.replace(')','')
    data = np.genfromtxt(StringIO(dtmp), dtype = ddtypes)
    all_data[atom_type][diff_type][diff_direction] = data

for atom_type in all_data.keys():
    #if all(["GB" in all_data[atom_type].keys(), "Bulk" in all_data[atom_type].keys()]):
        for diff_dir in all_data[atom_type]["GB"].keys():
            data_GB = all_data[atom_type]["GB"][diff_dir]
            fig, ax = plt.subplots()
            xs = 1 / (kB * data_GB["T"])
            ys = np.log(data_GB["D"])
            ax.plot(xs, data_GB["D"], '.', color = '#1b9e77', label = "GB")
            neg_Q_gb, ln_D0_gb = np.polyfit(xs, ys, 1)
            ax.plot(xs, np.exp(ln_D0_gb + neg_Q_gb * xs), '--', color = '#1b9e77', alpha = 0.3)
            annotation = f"$Q^{{GB}}$ = {-neg_Q_gb:.2f} eV, $D_0^{{GB}}$ = {np.exp(ln_D0_gb):.2E}"
            if "Bulk" in all_data[atom_type].keys() and diff_dir in all_data[atom_type]["Bulk"].keys():
                data_Bulk = all_data[atom_type]["Bulk"][diff_dir]
                xs = 1 / (kB * data_Bulk["T"])
                ys = np.log(data_Bulk["D"])
                ax.plot(xs, data_Bulk["D"], '.', color = '#d95f02', label = "Bulk")
                if len(data_Bulk) > 2:
                    neg_Q_bulk, ln_D0_bulk = np.polyfit(xs, ys, 1)
                    ax.plot(xs, np.exp(ln_D0_bulk + neg_Q_bulk * xs), '--', color = '#d95f02', alpha = 0.3)
                    annotation += f"\n$Q^{{Bulk}}$ = {-neg_Q_bulk:.2f} eV, $D_0^{{Bulk}}$ = {np.exp(ln_D0_bulk):.2E}"
            ax.text(0.5, 0.7, annotation, transform = ax.transAxes)
            ax.set_xlabel(r"1/$k_BT$ ($eV^{{-1}}$)")
            ax.set_ylabel(r"D ($m^2/s$)")
            ax2 = ax.secondary_xaxis('top', functions = (ArrtoT, TtoArr))
            ax.tick_params(which = 'both', bottom = True, left = True, right = True, direction = "in")
            ax2.tick_params(direction = "in")
            ax2.set_xlabel('Temperature (K)')
            ax.set_yscale('log')
            ax.set_title(f"{atom_type} {diff_dir.upper()} Diffusion")
            plt.legend(loc = 'best')
            plt.savefig(f"{atom_type}_{diff_dir}_diffusion_results.png")
            plt.close()
        # Also plot the xy and z diffusion behavior on the same plot
        xy_gb_data = all_data[atom_type]["GB"]["xy"]
        z_gb_data = all_data[atom_type]["GB"]["z"]
        if "Bulk" in all_data[atom_type].keys():
            fig, ax = plt.subplots(1, 2, figsize = (12,4))
            plt.subplots_adjust(wspace = 0.02)
            xy_bulk_data = all_data[atom_type]["Bulk"]["xy"]
            z_bulk_data = all_data[atom_type]["Bulk"]["z"]
            ax[1].set_title("Bulk Diffusion")
            ax[1].plot(1 / (kB * xy_bulk_data["T"]), xy_bulk_data["D"], '.', color = '#1b9e77', label = "XY")
            ax[1].plot(1 / (kB * z_bulk_data["T"]), z_bulk_data["D"], '.', color = '#d95f02', label = "Z")
            ax[1].set_xlabel(r"1/$k_BT$ ($eV^{{-1}}$)")
            ax[1].set_yscale('log')
            ax[1].tick_params(which = 'both', bottom = True, left = True, right = True, labelleft = False, direction = "in")
            ax2 = ax[1].secondary_xaxis('top', functions = (ArrtoT, TtoArr))
            ax2.set_xlabel('Temperature (K)')
            ax2.tick_params(direction = "in")

            ax[0].set_title("GB Diffusion")
            ax[0].plot(1 / (kB * xy_gb_data["T"]), xy_gb_data["D"], '.', color = '#1b9e77', label = "XY")
            ax[0].plot(1 / (kB * z_gb_data["T"]), z_gb_data["D"], '.', color = '#d95f02', label = "Z")
            ax[0].set_xlabel(r"1/$k_BT$ ($eV^{{-1}}$)")
            ax[0].set_ylabel(r"D ($m^2/s$)")
            ax[0].set_yscale('log')
            ax[0].tick_params(which = 'both', bottom = True, left = True, right = True, direction = "in")
            ax2 = ax[0].secondary_xaxis('top', functions = (ArrtoT, TtoArr))
            ax2.set_xlabel('Temperature (K)')
            ax2.tick_params(direction = "in")
            ax[0].legend(loc = 'best')
            ax[1].set_ylim(ax[0].get_ylim())
        else:
            fig, ax = plt.subplots()
            ax.set_title("GB Diffusion")
            ax.plot(1 / (kB * xy_gb_data["T"]), xy_gb_data["D"], '.', color = '#1b9e77', label = "XY")
            ax.plot(1 / (kB * z_gb_data["T"]), z_gb_data["D"], '.', color = '#d95f02', label = "Z")
            ax.set_xlabel(r"1/$k_BT$ ($eV^{{-1}}$)")
            ax.set_ylabel(r"D ($m^2/s$)")
            ax.set_yscale('log')
            ax.tick_params(which = 'both', bottom = True, left = True, right = True, direction = "in")
            ax2 = ax.secondary_xaxis('top', functions = (ArrtoT, TtoArr))
            ax2.set_xlabel('Temperature (K)')
            ax2.tick_params(direction = "in")
            ax.legend(loc = 'best')

        plt.savefig(f"{atom_type}_xy_z_diff_comparison.png", bbox_inches = 'tight')
