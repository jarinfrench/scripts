#! /usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd

kB = 8.617e-5

data = pd.read_csv('diffusion_results.txt', sep=' ', skiprows=1, names=(
    'concentration', 'temperature', 'run', 'msd_xyz', 'std_err'),
    index_col=['concentration', 'temperature'])

group_count = data.groupby(['concentration', 'temperature']).count()
group_count.columns = ['Count', 'i1', 'i2']
group_count.drop(['i1', 'i2'], axis=1, inplace=True)
summed_vals = data.groupby(['concentration', 'temperature']).sum()
summed_vals.columns = ['Run Sum', 'MSD Sum', 'StdErr Sum']
summed_vals.drop(['Run Sum', 'StdErr Sum'], axis=1, inplace=True)
sqrt_sum_sq = data.pow(2).groupby(
    ['concentration', 'temperature']).sum().pow(0.5)
sqrt_sum_sq.columns = ['i1', 'i2', 'StdErr Calc']
sqrt_sum_sq.drop(['i1', 'i2'], axis=1, inplace=True)
msd_std_dev = data.groupby(['concentration', 'temperature']).std()
msd_std_dev.columns = ['i1', 'StDev', 'i2']
msd_std_dev.drop(['i1', 'i2'], axis=1, inplace=True)
result = group_count.join(summed_vals).join(
    sqrt_sum_sq).join(msd_std_dev).reset_index()
result['MSD Average'] = result['MSD Sum'] / result['Count']
result['StdErr Result'] = result['StdErr Calc'] / result['Count']
result['StdErr Alt'] = result['StDev'] / result['Count'].pow(0.5)

for T in data.reset_index().temperature.unique():
    rtmp = result.query("temperature == @T")
    x = rtmp.concentration.to_list()
    y = rtmp['MSD Average'].to_list()
    yerr = rtmp['StdErr Alt'].to_list()
    plt.errorbar(x, y, yerr, fmt='.')
    plt.yscale('log')
    plt.savefig(f"T{T}_defect_diffusion_vs_composition.png")
    plt.clf()

for c in data.reset_index().concentration.unique():
    rtmp = result.query("concentration == @c")
    x = rtmp.temperature.to_list()
    x_Arrhenius = [1 / (kB * i) for i in x]
    y = rtmp['MSD Average'].to_list()
    yerr = rtmp['StdErr Alt'].to_list()
    plt.errorbar(x_Arrhenius, y, yerr, fmt='.', label=f"U{c}Mo")
plt.yscale('log')
plt.legend()
plt.savefig("defect_diffusion_vs_temperature.png")
