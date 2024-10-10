#! /usr/bin/env python3

import argparse

import numba as nb
import numpy as np
import pandas as pd


# Based off of https://stackoverflow.com/a/73643899
@nb.jit(nopython=True)
def conditional_assignment(arr, res, threshold):
    length = len(arr)
    for i in range(length):
        if arr[i] >= threshold:
            res[i] = True
        else:
            res[i] = False
    return res


def percent2float(value: str):
    return float(value.strip("%")) / 100


def createHeaderRows(header_rows, new_columns):
    for i in range(len(header_rows)):
        header_rows[i].append(new_columns[0].split(",")[i])
        header_rows[i].append(new_columns[1].split(",")[i])

    return header_rows


parser = argparse.ArgumentParser(
    usage="%(prog)s [-h] file",
    description="Generate the mobility database given the specified file",
)
parser.add_argument("file", help="File containing the slope data")
parser.add_argument(
    "-t",
    "--threshold",
    type=float,
    default=0.25,
    help="Threshold for including slope values in calculations (default: 0.25)",
)
parser.add_argument(
    "--temperatures",
    nargs="+",
    type=int,
    default="-1",
    help="The desired temperatures to use",
)
parser.add_argument(
    "--radii", nargs="+", type=int, default="-1", help="The desired radii to use"
)

args = parser.parse_args()

converters = {
    "axis": int,
    "potential": str,
    "concentration": float,
    "impurity": str,
    "misorientation": float,
    "Temperature": int,
    "r0": int,
    "run_number": int,
    "slope": float,
    "growth": percent2float,
}

df = pd.read_csv(
    args.file,
    sep=" ",
    names=tuple(converters.keys()),
    converters=converters,
)
if args.temperatures == -1:
    args.temperatures = df.Temperature.unique()
if args.radii == -1:
    args.radii = df.r0.unique()

res = np.empty(len(df), dtype="float")
df["enough growth"] = conditional_assignment(
    df["growth"].values, res, args.threshold
).astype("bool")

df["mobility"] = abs(df["slope"] * 1e-8 / (2 * np.pi))

multi_index = [
    "axis",
    "potential",
    "concentration",
    "impurity",
    "misorientation",
    "Temperature",
    "r0",
    "run_number",
]
df = df.set_index(multi_index).sort_values(multi_index)

# get the average mobility per different r0
avg_by_indices_separate_r0 = [
    "axis",
    "potential",
    "concentration",
    "impurity",
    "misorientation",
    "Temperature",
    "r0",
    "enough growth",
]
avg_mobilities_separate_r0 = (
    df.groupby(avg_by_indices_separate_r0).mobility.mean().reset_index()
)
std_errs_separate_r0 = (
    df.groupby(avg_by_indices_separate_r0).mobility.sem().reset_index()
)

# Do the same, but averaging over all radii
avg_by_indices_averaged_r0 = [
    "axis",
    "potential",
    "concentration",
    "impurity",
    "misorientation",
    "Temperature",
    "enough growth",
]
avg_mobilities_averaged_r0 = (
    df.query("r0 in @args.radii")
    .groupby(avg_by_indices_averaged_r0)
    .mobility.mean()
    .reset_index()
)
std_errs_averaged_r0 = (
    df.query("r0 in @args.radii")
    .groupby(avg_by_indices_averaged_r0)
    .mobility.sem()
    .reset_index()
)

# Make a separate dataframe for the relevant data
# This dataframe has each r0 separated out
result_df_separate_r0 = pd.DataFrame(avg_mobilities_separate_r0)
result_df_separate_r0["StdErr"] = std_errs_separate_r0["mobility"]

# This dataframe is the average of all initial radii
result_df_averaged_r0 = pd.DataFrame(avg_mobilities_averaged_r0)
result_df_averaged_r0["StdErr"] = std_errs_averaged_r0["mobility"]

m_vs_T_indices = ["axis", "misorientation", "Temperature"]
r1 = result_df_separate_r0.set_index(m_vs_T_indices).sort_values(m_vs_T_indices)
r2 = result_df_averaged_r0.set_index(m_vs_T_indices).sort_values(m_vs_T_indices)


r1 = r1.query("Temperature in @args.temperatures and r0 in @args.radii")
r2 = r2.query("Temperature in @args.temperatures")
# Write the MvsT data for separate r0...
for ax in result_df_separate_r0.axis.unique():
    for mis in result_df_separate_r0.misorientation.unique():
        header_rows = [["Temperature"], ["K"], [""], [""], [""]]
        tmp_df = pd.DataFrame()
        for pot in r1.potential.unique():
            for c in r1.concentration.unique():
                for imp in r1.impurity.unique():
                    for r0 in np.sort(r1["r0"].unique())[::-1]:
                        c_label = (
                            "Pure UO2" if (c == 0 and imp == "Pure") else f"{c}%{imp}"
                        )
                        col_headers = [
                            i + f",m^2/s,{pot},{c_label},r0={r0}"
                            for i in ["M*", "StdErr"]
                        ]
                        condition = (
                            (r1["potential"] == pot)
                            & (r1["concentration"] == c)
                            & (r1["impurity"] == imp)
                            & (r1["r0"] == r0)
                            & (r1["enough growth"])
                        )
                        try:
                            if condition.sum() == 0:  # Check if all rows are False
                                continue
                            selected_df = r1[condition].loc[ax, mis]
                            header_rows = createHeaderRows(header_rows, col_headers)
                            tmp_df = pd.concat(
                                (
                                    tmp_df,
                                    pd.DataFrame(
                                        {
                                            col_headers[0]: selected_df.mobility,
                                            col_headers[1]: selected_df.StdErr,
                                        }
                                    ),
                                ),
                                axis=1,
                            )
                        except KeyError:
                            continue
        if tmp_df.empty:
            continue
        mis_rep = "sigma7" if mis == 38.20 else str(int(mis)) + "degree"
        with open(f"{ax}_{mis_rep}_M_vs_T_separate_r0.csv", "w") as f:
            for i in header_rows:
                f.write(",".join(i))
                f.write("\n")
        tmp_df.sort_index().to_csv(
            f"{ax}_{mis_rep}_M_vs_T_separate_r0.csv",
            header=False,
            na_rep="--",
            mode="a",
        )

# ... and averaged r0
for ax in result_df_averaged_r0.axis.unique():
    for mis in result_df_averaged_r0.misorientation.unique():
        header_rows = [["Temperature"], ["K"], [""], [""]]
        tmp_df = pd.DataFrame()
        for pot in r2.potential.unique():
            for c in r2.concentration.unique():
                for imp in r2.impurity.unique():
                    c_label = "Pure UO2" if (c == 0 and imp == "Pure") else f"{c}%{imp}"
                    col_headers = [
                        i + f",m^2/s,{pot},{c_label}" for i in ["M*", "StdErr"]
                    ]
                    condition = (
                        (r2["potential"] == pot)
                        & (r2["concentration"] == c)
                        & (r2["impurity"] == imp)
                        & (r2["enough growth"])
                    )
                    try:
                        if condition.sum() == 0:
                            continue
                        selected_df = r2[condition].loc[ax, mis]
                        header_rows = createHeaderRows(header_rows, col_headers)
                        tmp_df = pd.concat(
                            (
                                tmp_df,
                                pd.DataFrame(
                                    {
                                        col_headers[0]: selected_df.mobility,
                                        col_headers[1]: selected_df.StdErr,
                                    }
                                ),
                            ),
                            axis=1,
                        )

                    except KeyError:
                        continue
        if tmp_df.empty:
            continue
        mis_rep = "sigma7" if mis == 38.20 else str(int(mis)) + "degree"
        with open(f"{ax}_{mis_rep}_M_vs_T_avg_r0.csv", "w") as f:
            for i in header_rows:
                f.write(",".join(i))
                f.write("\n")
        tmp_df.sort_index().to_csv(
            f"{ax}_{mis_rep}_M_vs_T_avg_r0.csv", header=False, na_rep="--", mode="a"
        )

# This is for the M vs c data
m_vs_c_indices = ["axis", "misorientation", "concentration"]
r1 = result_df_separate_r0.set_index(m_vs_c_indices).sort_values(m_vs_c_indices)
r2 = result_df_averaged_r0.set_index(m_vs_c_indices).sort_values(m_vs_c_indices)
r1 = r1.query("Temperature in @args.temperatures")
r2 = r2.query("Temperature in @args.temperatures")
# Write the dataframes out
for ax in result_df_separate_r0.axis.unique():
    for mis in result_df_separate_r0.misorientation.unique():
        header_rows = [["Concentration"], ["at%"], [""], [""]]
        tmp_df = pd.DataFrame()
        for pot in r1.potential.unique():
            for T in r1.Temperature.unique():
                for r0 in np.sort(r1["r0"].unique())[::-1]:
                    col_headers = [
                        i + ",m^2/s," + f"{pot},{T}K,r0={r0}" for i in ["M*", "StdErr"]
                    ]
                    condition = (
                        (r1["potential"] == pot)
                        & (r1["Temperature"] == T)
                        & (r1["r0"] == r0)
                        & (r1["enough growth"])
                    )
                    try:
                        if condition.sum() == 0:
                            continue
                        selected_df = r1[condition].loc[ax, mis]
                        header_rows = createHeaderRows(header_rows, col_headers)
                        tmp_df = pd.concat(
                            (
                                tmp_df,
                                pd.DataFrame(
                                    {
                                        col_headers[0]: selected_df.mobility,
                                        col_headers[1]: selected_df.StdErr,
                                    }
                                ),
                            ),
                            axis=1,
                        )
                    except KeyError:
                        continue
    if tmp_df.empty:
        continue
    mis_rep = "sigma7" if mis == 38.20 else str(int(mis)) + "degree"
    with open(f"{ax}_{mis_rep}_M_vs_c_separate_r0.csv", "w") as f:
        for i in header_rows:
            f.write(",".join(i))
            f.write("\n")
    tmp_df.sort_index().to_csv(
        f"{ax}_{mis_rep}_M_vs_c_separate_r0.csv",
        header=False,
        na_rep="--",
        mode="a",
    )

for ax in result_df_averaged_r0.axis.unique():
    for mis in result_df_averaged_r0.misorientation.unique():
        header_rows = [["Concentration"], ["at%"], [""], [""]]
        tmp_df = pd.DataFrame()
        for pot in r2.potential.unique():
            for T in r2.Temperature.unique():
                col_headers = [i + ",m^2/s," + f"{pot},{T}K" for i in ["M*", "StdErr"]]
                condition = (
                    (r2["potential"] == pot)
                    & (r2["Temperature"] == T)
                    & (r2["enough growth"])
                )
                try:
                    if condition.sum() == 0:
                        continue
                    selected_df = r2[condition].loc[ax, mis]
                    header_rows = createHeaderRows(header_rows, col_headers)
                    tmp_df = pd.concat(
                        (
                            tmp_df,
                            pd.DataFrame(
                                {
                                    col_headers[0]: selected_df.mobility,
                                    col_headers[1]: selected_df.StdErr,
                                }
                            ),
                        ),
                        axis=1,
                    )
                except KeyError:
                    continue
        if tmp_df.empty:
            continue
        mis_rep = "sigma7" if mis == 38.20 else str(int(mis)) + "degree"
        with open(f"{ax}_{mis_rep}_M_vs_c_avg_r0.csv", "w") as f:
            for i in header_rows:
                f.write(",".join(i))
                f.write("\n")
        tmp_df.sort_index().to_csv(
            f"{ax}_{mis_rep}_M_vs_c_avg_r0.csv", header=False, na_rep="--", mode="a"
        )
