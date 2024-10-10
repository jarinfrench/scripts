#! /usr/bin/env python

import os
import numpy as np
from math import log
import scipy.optimize as optimize
import argparse
import matplotlib.pyplot as plt

def vf(data, m0, q, T0):
    x = np.array(data)
    return m0 * np.exp(-q / (8.617e-5 * (x - T0)))

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] file(s) [options]',
    description = "Calculates the Vogel-Fulcher parameters that fit the data in 'file'")
parser.add_argument('file', help = "The file to process")
parser.add_argument('-l', '--lower-bound', type = float, default = 0, help = "The lower temperature bound to fit to")
parser.add_argument('-u', '--upper-bound', type = float, default = 1000000, help = "The upper temperature bound to fit to")
parser.add_argument('-T', '--fix-T', type = float, help = "Fix the temperature during fitting to this value")
parser.add_argument('--compare', type = float, nargs = '+', help = "Plots the specified (fixed) T0 values on the same plot")
parser.add_argument('--ignore-error', action = 'store_true', help = "Ignore the error when calculating the fit.")

args = parser.parse_args()

with open(args.file, 'r') as f:
    x = []
    y = []
    y_err = []
    for line in f:
        data = line.split()
        try:
        # if float(data[0]) >= args.lower_bound and float(data[0]) <= args.upper_bound:
            y_err.append(float(data[2]))
            y.append(float(data[1]))
            x.append(float(data[0]))
        except:
            continue

    allowed = [idx for idx, i in enumerate(x) if i >= args.lower_bound and i <= args.upper_bound]
    x_allowed = [x[i] for i in allowed]
    y_allowed = [y[i] for i in allowed]
    y_err_allowed = [y_err[i] for i in allowed]
    guess = (1.0e-7, 0.1, 1000)
    guess_fix_T = (1.0e-7, 0.1)
    if args.compare:
        all_params = []
        for T0 in args.compare:
            if args.ignore_error:
                params, _ = optimize.curve_fit(lambda x, m0, q: vf(x, m0, q, T0), x_allowed, y_allowed, guess_fix_T, maxfev = 1000000)
            else:
                params, _ = optimize.curve_fit(lambda x, m0, q: vf(x, m0, q, T0), x_allowed, y_allowed, guess_fix_T, sigma = y_err_allowed, maxfev = 1000000)
            tmp = list(params)
            tmp.append(T0)
            params = tuple(tmp)
            all_params.append(params)
    else:
        if args.ignore_error:
            if not args.fix_T:
                params, _ = optimize.curve_fit(vf, x_allowed, y_allowed, guess, maxfev = 1000000)
            else:
                params, _ = optimize.curve_fit(lambda x, m0, q: vf(x, m0, q, args.fix_T), x_allowed, y_allowed, guess_fix_T, maxfev = 1000000)
        else:
            if not args.fix_T:
                params, _ = optimize.curve_fit(vf, x_allowed, y_allowed, guess, sigma = y_err_allowed, maxfev = 1000000)
            else:
                params, _ = optimize.curve_fit(lambda x, m0, q: vf(x, m0, q, args.fix_T), x_allowed, y_allowed, guess_fix_T, sigma = y_err_allowed, maxfev = 1000000)

        if args.fix_T:
            tmp = list(params)
            tmp.append(args.fix_T)
            params = tuple(tmp)

        print("M0 = {}, Q = {}, T0 = {}".format(*params))

    fig = plt.figure()
    fig.add_subplot(111)
    ax = plt.gca()
    x_conv = [1/(8.617e-5 * i) for i in x] # 1/kB*T
    x_fit = [i for i in np.arange(min(x), max(x)) if i >= args.lower_bound and i <= args.upper_bound]
    x_conv_fit = [1/(8.617e-5 * i) for i in x_fit]
    y_conv = [log(i) for i in y] # ln(M)
    if args.compare:
        ax.set_title("Vogel-Fulcher fit to GB Mobility $T_0$ Comparison")
        ax.set_xlabel("1/$k_B*T$ ($eV^{-1}$)")
        ax.set_ylabel("log(Mobility ($m^2/s$))")
        with open(os.path.splitext(args.file)[0] + "_fitted_params.txt", 'w') as f:
            for params in all_params:
                y_fit = vf(x_fit, *params)
                y_fit_conv = [log(i) for i in y_fit]
                residuals = y_allowed - vf(x_allowed, *params)
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((y_allowed - np.mean(y_allowed))**2)
                r_sq = 1 - (ss_res / ss_tot)
                ax.plot(x_conv, y_conv, 'ro', x_conv_fit, y_fit_conv, 'k-')
                f.write(f"M0 = {params[0]:.4e}, Q = {params[1]:.4f}, T0 (fixed) = {params[2]:.2f}, Rsq = {r_sq:.4f}\n")
        plt.savefig(os.path.splitext(args.file)[0] + "_T0_comparison.png")

    else:
        y_fit = vf(x_fit, *params)
        y_fit_conv = [log(i) for i in y_fit]
        residuals = y_allowed - vf(x_allowed, *params)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y_allowed - np.mean(y_allowed))**2)
        r_sq = 1 - (ss_res / ss_tot)
        txt = f"$M_0$ = {params[0]:.4e}\nQ = {params[1]:.4f}\n$T_0$ = {params[2]:.2f}\n$R^2$ = {r_sq:.4f}"
        ax.plot(x_conv, y_conv, 'ro', x_conv_fit, y_fit_conv, 'k-')
        ax.set_title("Vogel-Fulcher fit to GB Mobility")
        ax.set_xlabel("1/$k_B*T$ ($eV^{-1}$)")
        ax.set_ylabel("log(Mobility ($m^2/s$))")
        ax.text(0.05, 0.1, txt, transform = ax.transAxes)
        plt.savefig(os.path.splitext(args.file)[0] + ".png")
