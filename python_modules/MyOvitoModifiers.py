from ovito.data import DataCollection
import numpy as np


# Calculates the center of the grain based on *_interface.dat files in the
# current directory. TODO: Make this independent of those files!
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


# Calculates the radius of an embedded cylindrical grain based on an area_data.txt file
# TODO: Make this independent of that file!
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


# Calculates the 1D, 2D, and 3D MSD for a given simulation
def CalculateMSDModifier(frame: int, data: DataCollection):
    displacement_magnitudes_all = data.particles["Displacement Magnitude"]
    displacement_magnitudes = data.particles["Displacement"]

    msd_all = np.sum(displacement_magnitudes_all**2) / len(
        displacement_magnitudes_all
    )
    num_elements = len(displacement_magnitudes)
    tmp = np.array([i[0] for i in displacement_magnitudes])
    msd_x = np.sum(np.square(tmp)) / num_elements
    tmp = np.array([i[1] for i in displacement_magnitudes])
    msd_y = np.sum(np.square(tmp)) / num_elements
    tmp = np.array([i[2] for i in displacement_magnitudes])
    msd_z = np.sum(np.square(tmp)) / num_elements

    data.attributes["MSD.XYZ"] = msd_all
    data.attributes["MSD.XY"] = msd_x + msd_y
    data.attributes["MSD.XZ"] = msd_x + msd_z
    data.attributes["MSD.YZ"] = msd_y + msd_z
    data.attributes["MSD.X"] = msd_x
    data.attributes["MSD.Y"] = msd_y
    data.attributes["MSD.Z"] = msd_z
