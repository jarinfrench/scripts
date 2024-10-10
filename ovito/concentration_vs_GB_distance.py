from ovito.io import import_file, export_file
from ovito.data import DataCollection, DataTable
import numpy as np
import argparse

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
            data.attributes["GrainCenter"] = np.array(center)
    except FileNotFoundError:
        data.attributes["GrainCenter.X"] = last_center[0]
        data.attributes["GrainCenter.Y"] = last_center[1]
        data.attributes["GrainCenter.Z"] = last_center[2]
        data.attributes["GrainCenter"] = np.array(last_center)


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


def DistanceFromGBModifier(frame: int, data: DataCollection):
    grain_center = data.attributes["GrainCenter"]
    cur_grain_radius = data.attributes["GrainRadius"]
    distance_from_GB = [None] * data.particles.count

    # Loop over every particle and calculate the distance from the grain center
    # Then subtract the current grain radius
    for i in range(data.particles.count):
        distance_from_GB = (
            data.cell.delta_vector(grain_center, data.particles.positions[i])
            - cur_grain_radius
        )

    data.particles_.create_property("DistanceFromGB", data=distance_from_GB)


def ConcentrationRingModifier(frame: int, data: DataCollection):
    global width
    dr = data.attributes["GrainCenter"] - np.array(
        [
            data.cell[0, 0],
            data.cell[1, 1],
            data.cell[2, 2],
        ]
    )
    distances = []
    if dr[0] > 0:
        distances.append(data.attributes["GrainCenter.X"])
    else:
        distances.append(data.cell[0, 0] - data.attributes["GrainCenter.X"])
    if dr[1] > 0:
        distances.append(data.attributes["GrainCenter.Y"])
    else:
        distances.append(data.cell[1, 1] - data.attributes["GrainCenter.Y"])

    n_rings = int(round(max(distances) / width))
    concentrations = []
    distances_to_GB = []

    for i in range(n_rings):
        r_min = i * width
        r_max = (i + 1) * width
        r_mid = (i + 0.5) * width
        r_min_sq = r_min * r_min
        r_max_sq = r_max * r_max
        selection = (
            (data.particles["Position"][:, 0] - data.attributes["GrainCenter.X"]) ** 2
            + (data.particles["Position"][:, 1] - data.attributes["GrainCenter.Y"]) ** 2
            > r_min_sq
        ) & (
            (data.particles["Position"][:, 0] - data.attributes["GrainCenter.X"]) ** 2
            + (data.particles["Position"][:, 1] - data.attributes["GrainCenter.Y"]) ** 2
            < r_max_sq
        )
        n_in_ring = np.sum(selection)
        selection = selection & (data.particles["Particle Type"] == 2)
        n_solute_in_ring = np.sum(selection)
        concentrations.append(n_solute_in_ring / n_in_ring * 100)
        distances_to_GB.append(data.attributes["GrainRadius"] - r_mid)

    table = DataTable(
        title="Concentration vs Distance from GB",
        identifier="c-vs-r",
        plot_mode=DataTable.PlotMode.Scatter,
    )
    table.x = table.create_property("Distance", data=distances_to_GB)
    table.y = table.create_property("Concentration", data=concentrations)
    data.objects.append(table)


parser = argparse.ArgumentParser(
    usage="%(prog)s [-h] r_in r_out",
    description="Calculate solute concentration in a ring over the whole simulation.",
)

parser.add_argument("width", type=float, help="Width of cylinder in angstroms")

args = parser.parse_args()
width = args.width

pipeline = import_file("*.dump")

# Calculate the center of the grain
pipeline.modifiers.append(GrainCenterModifier)

# Calculate the approximate radius of the grain
pipeline.modifiers.append(GrainRadiusModifier)

# Calculate the concentration in each ring from the center of the grain to the
# edge of the cell
pipeline.modifiers.append(ConcentrationRingModifier)

for frame in range(pipeline.source.num_frames):
    data = pipeline.compute(frame)
    table = data.tables["c-vs-r"]
    export_file(
        table, f"concentration_vs_GB_distance.{frame}.txt", "txt/table", frame=frame
    )
