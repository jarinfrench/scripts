#! /usr/bin/env python3

import argparse
import colorama as clr
import matplotlib.pyplot as plt
import matplotlib.text as mtext
import numpy as np
import sys, pwlf

from scipy.stats import iqr
from GPyOpt.methods import BayesianOptimization
from sklearn.linear_model import ElasticNetCV

def addInfo(file):
    with open(args.output, 'a') as f:
        f.write(f"No results calculated for {file}\n")

def press(event):
    if event.key == 'y':
        ax.set_title("")
        plt.savefig("mobility_plot.png", bbox_inches = 'tight')
        plot_text = f"{args.file} dA/dt = {avg_params[0]:.2f} error = {avg_slope_err:.2f} r_sq = {avg_rsq:.4f} fit to points {index_start + 1}:{index_end + 1}\n"
        with open(args.output, 'a') as f:
            f.write(plot_text)
        sys.exit(0)
    elif event.key == 'n':
        print(clr.Fore.RED  + clr.Style.BRIGHT + "Image not saved" + clr.Style.RESET_ALL)
        sys.exit(0)
    else:
        return # do nothing

def calculateFit(x, y):
    if len(x) < 3:
        raise TypeError("Not enough data points to calculate error")
    # params contains the slope and the intercept, in that order
    params, sum_sq_res, *_ = np.polyfit(x, y, 1, full = True)
    _, cov_mat = np.polyfit(x, y, 1, cov = True)
    avg = np.mean(y)
    sum_sq = np.sum((y - avg)**2)
    r_sq = 1 - sum_sq_res[0] / sum_sq
    slope_err = np.sqrt(cov_mat[0][0])
    return params, r_sq, slope_err

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] file', description = "Creates a plot of the area data for grain growth, and calculates the slope of the area vs time plot.")
parser.add_argument('file', help = "The data file containing the area vs time data.")
parser.add_argument('-o', '--output', default = "slope_calc.txt", help = "File name for the data file containing the fitting results. Default: slope_calc.txt")
parser.add_argument('-m', '--image-name', default = "mobility_plot.png", help = "File name for the image file. Default: mobility_plot.png")
parser.add_argument('-n', '--segments', type = int, default = 6, help = "The initial number of line segments to fit to the data. Default = 6")
parser.add_argument('-p', '--percent', type = float, default = 0.25, help = "The minimum percent growth for calculations to occur. Default = 0.25")
parser.add_argument('-f', '--force', action = "store_true", help = "Force the usage of the specified number of segments")
parser.add_argument('-F', '--force-fit', action = "store_true", help = "Force fitting the data, even if not enough growth")
parser.add_argument('-r', '--fit-range', nargs = 2, type = int, default = [-1, -2], help = "First and last breakpoint to include in overall fit. The '0' breakpoint is the start of the data, and the '-1' breakpoint is the end of the data. (Default: determined by script)")
parser.add_argument('-g', '--growth', action = "store_true", help = "Flag to only show the percent growth")
parser.add_argument('-i', '--interactive', action = "store_true", help = "Flag to do interactive calculations")

args = parser.parse_args()
clr.init()

# First, we need to read the data
t = []; y1 = []; y2 = [] # initialize the data arrays

with open(args.file) as f:
    for line in f:
        if line.startswith("#"):
            continue # ignore commented lines
        data = line.split()
        t.append(float(data[0]))
        y1.append(float(data[1]))
        y2.append(float(data[2]))

t = np.array(t)
y1 = np.array(y1)
y2 = np.array(y2)
# get the initial area
initial_grain_size = [y1[0],y2[0]] # we don't know which data set is the important one yet
end_grain_size = [y1[-1],y2[-1]]

if initial_grain_size[0] - end_grain_size[0] < 0: # end size is bigger than initial size
    initial_grain_size = y2[0]
    end_grain_size = y2[-1]
    y = y2
else:
    initial_grain_size = y1[0]
    end_grain_size = y1[-1]
    y = y1

# Check if there was enough growth overall
percent_growth = 1 - end_grain_size / initial_grain_size

if percent_growth < args.percent:
    print(clr.Fore.RED + clr.Style.BRIGHT + f"Not enough overall growth ({percent_growth*100:.2f}%)" + clr.Style.RESET_ALL)
    if args.growth:
        sys.exit(0)
    if not args.force_fit:
        addInfo(args.file)
        sys.exit(0)
else:
    print(clr.Fore.GREEN + clr.Style.BRIGHT + f"Overall growth: {percent_growth*100:.2f}%" + clr.Style.RESET_ALL)
    if args.growth:
        sys.exit(0)

# get the model
model = pwlf.PiecewiseLinFit(t, y)

model.fit(args.segments)
if not args.force:
    num_points_in_lines = np.zeros(model.slopes.shape[0])
    for i,slope in enumerate(model.slopes):
        num_points_in_lines[i] = len([idx for idx, k in enumerate(t) if k >= model.fit_breaks[i] and k <= model.fit_breaks[i+1]])
    while [i for i in num_points_in_lines if i <= 3] or [i for i in model.slopes if i > 1]:
        args.segments -= 1
        if args.segments == 1:
            args.segments += 1
            break
        model.fit(args.segments)
        num_points_in_lines = np.zeros(model.slopes.shape[0])
        for i,slope in enumerate(model.slopes):
            num_points_in_lines[i] = len([idx for idx, k in enumerate(t) if k >= model.fit_breaks[i] and k <= model.fit_breaks[i+1]])
slopes = model.slopes

num_points_in_lines = np.zeros(model.slopes.shape[0])
time_length = np.zeros(slopes.shape[0])
area_length = np.zeros(slopes.shape[0])
annotation = ""
ignored = []
if args.fit_range[0] != -1: # not the default
    keep_start = model.fit_breaks[args.fit_range[0]]
else:
    keep_start = -1
if args.fit_range[1] != -2: # not the default
    if args.fit_range[1] >= len(model.fit_breaks):
        print(f"Indicated fitting to endpoint {args.fit_range}, but only {len(model.fit_breaks)} total breaks. Using final endpoint.")
        keep_end = model.fit_breaks[-1]
    else:
        keep_end = model.fit_breaks[args.fit_range[1]]

# Check for outliers
outlier_upper = np.percentile(model.slopes, 75) + 1.5 * iqr(model.slopes) # outliers are defined as values outside of 1.5 times the interquartile range
outlier_lower = np.percentile(model.slopes, 25) - 1.5 * iqr(model.slopes)
cumulative_growth = 0
for i, slope in enumerate(slopes):
    valid_data = [idx for idx, k in enumerate(t) if k >= model.fit_breaks[i] and k <= model.fit_breaks[i + 1]]
    num_points_in_lines[i] = len(valid_data)
    time_length[i] = (t[valid_data[-1]] - t[valid_data[0]]) / t[-1] * 100 # percent of total time over the specified range
    area_length[i] = (y[valid_data[0]] - y[valid_data[-1]])   / y[0] * 100 # percent growth observed in the specified range
    print(f"Regime {i+1}: {time_length[i]:.2f}% of total time, {area_length[i]:.2f}% grain shrinkage ({int(num_points_in_lines[i])} points), slope = {slopes[i]:.2f}")
    annotation += rf'''$\left(\frac{{dA}}{{dt}}\right)_{{{i+1}}}^{{{int(valid_data[0])}:{int(valid_data[-1])}}} \approx$ {slope:.2f}
''' # we force a newline here
    if len(slopes) > 2:
        if ((slope > outlier_upper or slope < outlier_lower) or  # current slope is an outlier
           ((abs(slope) < args.percent * max(y) / max(t) and i > 0))  # current slope is not the first slope, and is less than the defined threshold
           and (cumulative_growth > 95 or observed_growth > 95)): # overall growth has not reached 95% of it's total growth.
            ignored.append(model.fit_breaks[i:i+2]) # store the start and end point of the ignored section (uses values of t)
            cumulative_growth += area_length[i]
            observed_growth = (y[0] - y[valid_data[-1]]) / (y[0] - y[-1])*100
        else:
            if args.fit_range[0] == -1 and keep_start == -1:
                keep_start = model.fit_breaks[i]
            if args.fit_range[1] == -2:
                keep_end = model.fit_breaks[i+1]
            cumulative_growth += area_length[i]
            observed_growth = (y[0] - y[valid_data[-1]]) / (y[0] - y[-1])*100
    else:
        # TODO: There is a bug here that can cause keep_end to equal keep_start, leading to an error in line 196 (processed is an empty list)
        if abs(slope) < args.percent * max(y) / max(t): # for the two-regime case, we ignore any section that has a growth rate that leads to less than 25% growth over the specified time
            ignored.append(model.fit_breaks[i:i+2])
        else:
            if args.fit_range[0] == -1 and keep_start == -1:
                keep_start = model.fit_breaks[i]
            if args.fit_range[1] == -2:
                keep_end = model.fit_breaks[i+1]
            keep_end = model.fit_breaks[i+1]

# Need a check where I see if all the data points are ignored.
tHat = np.linspace(min(t), max(t), num=10000)
yHat = model.predict(tHat)

index_start = next((id for id,pt in enumerate(t) if pt >= keep_start))
try:
    index_end = next((id for id,pt in enumerate(t) if pt > keep_end))
except StopIteration:
    index_end = len(t)

fig, ax = plt.subplots()

if args.interactive:
    fig.canvas.mpl_connect('key_press_event', press) # allows key presses to be tracked
    ax.set_title("Press y to save, n otherwise")

if ignored: # check for empty list
    tHat_avg_ignored = np.concatenate([np.linspace(i,j,int((j-i)/10)) for i,j in ignored]).ravel() # creates an array of the ignored values
    yHat_avg_ignored = model.predict(tHat_avg_ignored)
    # ax.plot(tHat_avg_ignored, yHat_avg_ignored, '-', color = "gray", alpha = 0.5, label="Excluded Data Fit") # the ignored data

processed = list(zip(t[index_start:index_end], y[index_start:index_end])) # the data we are keeping
skipped = list(set(zip(t,y)) - set(processed)) # get the data points not in the processed list
avg_params, avg_rsq, avg_slope_err = calculateFit(*zip(*processed))
tHat_avg_kept = np.linspace(keep_start,keep_end) # creates an array of the kept values
yHat_avg_kept = tHat_avg_kept * avg_params[0] + avg_params[1]
ax.plot(*zip(*processed), 'r.', label="Fitted Data") # the data
if skipped:
    ax.plot(*zip(*skipped), '.', color="gray", alpha=0.5, label="Non-fitted Data")
ax.plot(tHat, yHat, '--', color = "black", label=f"Piecewise Linear Fit ({args.segments} segments)") # all piecewise linear fits
ax.plot(tHat_avg_kept, yHat_avg_kept, 'k-', label="Full fitted slope") # the average slope without the ignored data
for breakpoint in model.fit_breaks[1:-1]:
    ax.axvline(x=breakpoint, color = 'gray', alpha = 0.1, ls=':', zorder = 1)
ax.set_ylim([0,1.05*max(y)])
ax.set_xlabel("Time (ps)")
ax.set_ylabel(r"Area ($\AA^2$)")

ax.annotate(annotation, (1,1), (5,0), xycoords='axes fraction', textcoords='offset points', va='top', zorder = 10)
legend = ax.legend(loc='best')
info_text = f"Slope: {avg_params[0]:.2f}$\pm${avg_slope_err:.2f}\n$r^2$: {avg_rsq:.4f}\nFitted points: {index_start + 1}:{index_end} ({index_end - index_start})"
fig.canvas.draw()
fig.tight_layout()
p = legend.get_window_extent() # gives two points, p0 is the lower left, p1 is the upper right
bbox = fig.get_window_extent()
if p.ymax < 0.5 * bbox.ymax: # legend is at the bottom of the plot
    if p.xmax < 0.80 * bbox.xmax: # left side of plot
        offset = mtext.OffsetFrom(legend, (0.0, 1.0)) # top left corner of legend
        ax.annotate(info_text, xy=(0,0), xycoords='figure fraction', xytext=(0,10), ha='left', textcoords=offset, zorder=9)
    else: # right side of plot
        offset = mtext.OffsetFrom(legend, (1.0, 1.0)) # top right corner of legend
        ax.annotate(info_text, xy=(0,0), xycoords='figure fraction', xytext=(0,10), ha='right', textcoords=offset, zorder=9)
else: # legend is at the top
    offset = mtext.OffsetFrom(legend, (1.0, 0.0)) # lower left corner of legend
    ax.annotate(info_text, xy=(0,0), xycoords='figure fraction', xytext=(0,-40), ha='right', textcoords=offset, zorder=9)

if not args.interactive:
    plt.savefig("mobility_plot.png", bbox_inches = 'tight')
    with open(args.output, 'a') as f:
        f.write(f"{args.file} dA/dt = {avg_params[0]:.2f} error = {avg_slope_err:.2f} r_sq = {avg_rsq:.4f} fit to points {index_start + 1}:{index_end + 1}\n")
else:
    plt.show()
