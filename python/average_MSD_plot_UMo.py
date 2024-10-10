#! /usr/bin/env python
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.style.use("/home/jarinf/projects/scripts/python/journal.mplstyle")
mpl.rcParams["xtick.minor.visible"] = False
mpl.rcParams["ytick.minor.visible"] = False

df_names = ("step", "msd_x", "msd_y", "msd_z", "msd_xyz")
df_dtypes = {
    "step": "int",
    "msd_x": "float",
    "msd_y": "float",
    "msd_z": "float",
    "msd_xyz": "float",
}

Ts = ["T" + str(i) for i in [1050, 1100, 1150, 1200, 1250, 1300, 1350]]

for T in Ts:
    pure_fn = f"/media/jarinf/Research1/U/grain_growth/gamma/100/ternary_eam/20degree/{T}/large_r/"
    impure_fn = f"20degree/{T}/large_r"
    data1 = [
        pd.read_csv(
            pure_fn + "dir_1/MSD_U.dat",
            sep=" ",
            names=df_names,
            dtype=df_dtypes,
            skiprows=2,
        ).assign(concentration=0)
    ] + [
        pd.read_csv(
            f"{c}at%/{impure_fn}/dir_final_1/MSD_all.dat",
            sep=" ",
            names=df_names,
            dtype=df_dtypes,
            skiprows=2,
        ).assign(concentration=c)
        for c in [1, 3]
    ]
    data2 = [
        pd.read_csv(
            pure_fn + "dir_2/MSD_U.dat",
            sep=" ",
            names=df_names,
            dtype=df_dtypes,
            skiprows=2,
        ).assign(concentration=0)
    ] + [
        pd.read_csv(
            f"{c}at%/{impure_fn}/dir_final_2/MSD_all.dat",
            sep=" ",
            names=df_names,
            dtype=df_dtypes,
            skiprows=2,
        ).assign(concentration=c)
        for c in [1, 3]
    ]
    data3 = [
        pd.read_csv(
            pure_fn + "dir_3/MSD_U.dat",
            sep=" ",
            names=df_names,
            dtype=df_dtypes,
            skiprows=2,
        ).assign(concentration=0)
    ] + [
        pd.read_csv(
            f"{c}at%/{impure_fn}/dir_final_3/MSD_all.dat",
            sep=" ",
            names=df_names,
            dtype=df_dtypes,
            skiprows=2,
        ).assign(concentration=c)
        for c in [1, 3]
    ]
    data1_U = [
        pd.read_csv(
            pure_fn + "dir_1/MSD_U.dat",
            sep=" ",
            names=df_names,
            dtype=df_dtypes,
            skiprows=2,
        ).assign(concentration=0)
    ] + [
        pd.read_csv(
            f"{c}at%/{impure_fn}/dir_final_1/MSD_U.dat",
            sep=" ",
            names=df_names,
            dtype=df_dtypes,
            skiprows=2,
        ).assign(concentration=c)
        for c in [1, 3]
    ]
    data2_U = [
        pd.read_csv(
            pure_fn + "dir_2/MSD_U.dat",
            sep=" ",
            names=df_names,
            dtype=df_dtypes,
            skiprows=2,
        ).assign(concentration=0)
    ] + [
        pd.read_csv(
            f"{c}at%/{impure_fn}/dir_final_2/MSD_U.dat",
            sep=" ",
            names=df_names,
            dtype=df_dtypes,
            skiprows=2,
        ).assign(concentration=c)
        for c in [1, 3]
    ]
    data3_U = [
        pd.read_csv(
            pure_fn + "dir_3/MSD_U.dat",
            sep=" ",
            names=df_names,
            dtype=df_dtypes,
            skiprows=2,
        ).assign(concentration=0)
    ] + [
        pd.read_csv(
            f"{c}at%/{impure_fn}/dir_final_3/MSD_U.dat",
            sep=" ",
            names=df_names,
            dtype=df_dtypes,
            skiprows=2,
        ).assign(concentration=c)
        for c in [1, 3]
    ]
    data1_Mo = [
        pd.read_csv(
            pure_fn + "dir_1/MSD_U.dat",
            sep=" ",
            names=df_names,
            dtype=df_dtypes,
            skiprows=2,
        ).assign(concentration=0)
    ] + [
        pd.read_csv(
            f"{c}at%/{impure_fn}/dir_final_1/MSD_Mo.dat",
            sep=" ",
            names=df_names,
            dtype=df_dtypes,
            skiprows=2,
        ).assign(concentration=c)
        for c in [1, 3]
    ]
    data2_Mo = [
        pd.read_csv(
            pure_fn + "dir_2/MSD_U.dat",
            sep=" ",
            names=df_names,
            dtype=df_dtypes,
            skiprows=2,
        ).assign(concentration=0)
    ] + [
        pd.read_csv(
            f"{c}at%/{impure_fn}/dir_final_2/MSD_Mo.dat",
            sep=" ",
            names=df_names,
            dtype=df_dtypes,
            skiprows=2,
        ).assign(concentration=c)
        for c in [1, 3]
    ]
    data3_Mo = [
        pd.read_csv(
            pure_fn + "dir_3/MSD_U.dat",
            sep=" ",
            names=df_names,
            dtype=df_dtypes,
            skiprows=2,
        ).assign(concentration=0)
    ] + [
        pd.read_csv(
            f"{c}at%/{impure_fn}/dir_final_3/MSD_Mo.dat",
            sep=" ",
            names=df_names,
            dtype=df_dtypes,
            skiprows=2,
        ).assign(concentration=c)
        for c in [1, 3]
    ]
    data_all = [
        pd.concat([data1[i], data2[i], data3[i]]).groupby(level=0).mean()
        for i in range(len(data1))
    ]
    data_U = [
        pd.concat([data1_U[i], data2_U[i], data3_U[i]]).groupby(level=0).mean()
        for i in range(len(data1_U))
    ]
    data_Mo = [
        pd.concat([data1_Mo[i], data2_Mo[i], data3_Mo[i]]).groupby(level=0).mean()
        for i in range(len(data1_Mo))
    ]
    data_all[0]["step"] = data_all[0]["step"] - 10000
    data_U[0]["step"] = data_U[0]["step"] - 10000
    data_Mo[0]["step"] = data_Mo[0]["step"] - 10000
    for df in data_all:
        df["msd_xy"] = df["msd_x"] + df["msd_y"]
        df["msd_xz"] = df["msd_x"] + df["msd_z"]
        df["msd_yz"] = df["msd_y"] + df["msd_z"]
    for df in data_U:
        df["msd_xy"] = df["msd_x"] + df["msd_y"]
        df["msd_xz"] = df["msd_x"] + df["msd_z"]
        df["msd_yz"] = df["msd_y"] + df["msd_z"]
    for df in data_Mo:
        df["msd_xy"] = df["msd_x"] + df["msd_y"]
        df["msd_xz"] = df["msd_x"] + df["msd_z"]
        df["msd_yz"] = df["msd_y"] + df["msd_z"]
    colors = [
        "#1B9E77",
        "#D95F02",
        "#7570B3",
        "#E7298A",
        "#66A61E",
        "#E6AB02",
        "#A6761D",
        "#666666",
    ]
    l1 = ["", "U_", "Mo_"]
    l2 = ["", "U", "Mo"]
    labels = {0: "Pure U", 1: "U1Mo", 3: "U3Mo"}
    # scale_factor = {0: 1333.33, 1: 1000, 3: 1433.33} # not sure what these represent
    resample = 300
    for d in ["x", "y", "z", "xy", "xz", "yz", "xyz"]:
        for j, ddata in enumerate([data_all, data_U, data_Mo]):
            cs = [int(np.mean(df["concentration"])) for df in ddata]
            fig, ax = plt.subplots()
            for i, df in enumerate(ddata):
                ax.plot(
                    df.iloc[::resample, :]["step"] * 0.002,  # / scale_factor[cs[i]],
                    df.iloc[::resample, :]["msd_" + d],
                    "o",
                    color=colors[i],
                    label=labels[cs[i]] + " " + d.upper() + " MSD",
                )
            # ax.axvline(x=1333.33, color = colors[0], linestyle = '--')
            # ax.axvline(x=1000.00, color = colors[1], linestyle = '--')
            # ax.axvline(x=1433.33, color = colors[2], linestyle = '--')
            ax.set_xlabel("Time (ps)")
            ax.set_ylabel(f"{l2[j]} MSD ($\AA^2$)")
            plt.legend(loc="best", frameon=False)
            plt.savefig(f"average_{l1[j]}MSD_{d}_comparison_{T}.png")
            plt.close(fig)
