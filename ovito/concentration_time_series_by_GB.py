#! /usr/bin/env python

from ovito.io import import_file, export_file
import ovito.modifiers as om
from ovito.data import DataCollection
import argparse
import numpy as np

last_center = [None, None, None]


def GrainCenterModifier(frame: int, data: DataCollection):
    global last_center
    timestep = data.attributes["Timestep"]
    center = [0, 0, 0]
    if last_center == [None, None, None]:
        last_center = [data.cell[0, 0] / 2, data.cell[1, 1] / 2, data.cell[2, 2] / 2]
    try:
        counts = 0

        with open(f"{timestep}_interface.dat") as f:
            var_names = f.readline().split(",")
            x_idx = [
                i for i, val in enumerate(var_names) if "X" in val and "Xu" not in val
            ][0]
            y_idx = [
                i for i, val in enumerate(var_names) if "Y" in val and "Yu" not in val
            ][0]
            z_idx = [
                i for i, val in enumerate(var_names) if "Z" in val and "Zu" not in val
            ][0]
            grain_num_idx = [
                i for i, val in enumerate(var_names) if "Grain Number" in val
            ][0]
            for line in f:
                if line.startswith("#"):
                    continue
                line_data = line.split()[0 : grain_num_idx + 1]

                grain_num = int(line_data[grain_num_idx])
                if grain_num == 2:
                    x = float(line_data[x_idx])
                    y = float(line_data[y_idx])
                    z = float(line_data[z_idx])
                    counts += 1
                    center = [i + j for i, j in zip(center, [x, y, z])]
            center = [val / counts for val in center]
            data.attributes["GrainCenter.X"] = center[0]
            data.attributes["GrainCenter.Y"] = center[1]
            data.attributes["GrainCenter.Z"] = center[2]
    except FileNotFoundError:
        data.attributes["GrainCenter.X"] = last_center[0]
        data.attributes["GrainCenter.Y"] = last_center[1]
        data.attributes["GrainCenter.Z"] = last_center[2]


def GrainRadiusModifier(frame: int, data: DataCollection):
    timestep = data.attributes["Timestep"]
    time = timestep * 0.002
    try:
        with open("area_data.txt", "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                line_data = line.split()
                if int(line_data[0]) == time:
                    data.attributes["GrainRadius"] = np.sqrt(
                        float(line_data[2]) / np.pi
                    )
    except FileNotFoundError:
        data.attributes["GrainRadius"] = None


parser = argparse.ArgumentParser(
    usage="%(prog)s [-h] r_in r_out",
    description="Calculate solute concentration in a ring over the whole simulation.",
)

parser.add_argument("width", type=float, help="Width of cylinder in angstroms")

args = parser.parse_args()
w = args.width

pipeline = import_file("*.dump")
ntypes = len(set(pipeline.compute().particles.particle_type))
if ntypes < 2:
    print("No solute atom types found (ntypes < 2)")
    exit()

# Calculate the center of the grain
pipeline.modifiers.append(GrainCenterModifier)
# Calculate the approximate radius of the grain
pipeline.modifiers.append(GrainRadiusModifier)
# Get the number of atoms in the ring
pipeline.modifiers.append(
    om.ExpressionSelectionModifier(
        expression=f"(Position.X-GrainCenter.X)^2+(Position.Y-GrainCenter.Y)^2"
        f" > (GrainRadius < {w/2} ? 0 : (GrainRadius-{w/2}))^2 && "
        f"(Position.X-GrainCenter.X)^2+(Position.Y-GrainCenter.Y)^2"
        f" < (GrainRadius < {w/2} ? {w} : (GrainRadius+{w/2}))^2"
    )
)
# Get the number of solute atoms in the ring
pipeline.modifiers.append(
    om.ExpressionSelectionModifier(
        expression=f"(Position.X-GrainCenter.X)^2+(Position.Y-GrainCenter.Y)^2"
        f" > (GrainRadius < {w/2} ? 0 : (GrainRadius-{w/2}))^2 && "
        f"(Position.X-GrainCenter.X)^2+(Position.Y-GrainCenter.Y)^2"
        f" < (GrainRadius < {w/2} ? {w} : (GrainRadius+{w/2}))^2 && ParticleType==2"
    )
)
# Get the time series evolution of these values
pipeline.modifiers.append(
    om.TimeSeriesModifier(
        operate_on=(
            "SelectExpression.num_selected",
            "SelectExpression.num_selected.2",
            "GrainCenter.X",
            "GrainCenter.Y",
            "GrainRadius",
        )
    )
)

data = pipeline.compute()
table = data.make_mutable(data.tables["time-series"])
series = table.y
table.y = table.create_property(
    "solute_percent", data=series[:, 1] / series[:, 0] * 100
)
export_file(table, "concentration_data.txt", "txt/table")
