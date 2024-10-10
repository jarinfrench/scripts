#! /usr/bin/env python

from ovito.io import import_file, export_file
import ovito.modifiers as om
from ovito.pipeline import FileSource
import numpy as np


def calculate_MSDs(frame, data):
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

    data.attributes["MSD"] = msd_all
    data.attributes["MSD X"] = msd_x
    data.attributes["MSD Y"] = msd_y
    data.attributes["MSD Z"] = msd_z


pipeline = import_file("*.dump")  # loads all dump files in the current directory
pipeline.modifiers.append(om.SelectTypeModifier(types={2}))  # select O atoms
pipeline.modifiers.append(
    om.DeleteSelectedModifier(operate_on={"particles"})
)  # delete them
disp_mod = om.CalculateDisplacementsModifier()
disp_mod.reference = FileSource()
disp_mod.reference.load("10000.dump")  # uses this file as the reference configuration
pipeline.modifiers.append(disp_mod)  # calculate the U atom displacements
pipeline.modifiers.append(calculate_MSDs)
export_file(
    pipeline,
    "MSD_U.dat",
    format="txt/attr",
    columns=["Timestep", "MSD X", "MSD Y", "MSD Z", "MSD"],
    multiple_frames=True,
)

# pipeline.modifiers.clear()
# pipeline.modifiers.append(om.SelectTypeModifier(types = {1})) # select U atoms
# pipeline.modifiers.append(om.DeleteSelectedModifier(operate_on = {'particles'})) # delete them
# pipeline.modifiers.append(om.CalculateDisplacementsModifier()) # calculate the U atom displacements
# pipeline.modifiers.append(calculate_MSDs)
# export_file(pipeline, "MSD_O.dat", format = "txt/attr", columns = ["Timestep", "MSD X", "MSD Y", "MSD Z", "MSD"], multiple_frames = True)
#
# pipeline.modifiers.clear()
# pipeline.modifiers.append(om.CalculateDisplacementsModifier()) # calculate the U atom displacements
# pipeline.modifiers.append(calculate_MSDs)
# export_file(pipeline, "MSD_all.dat", format = "txt/attr", columns = ["Timestep", "MSD X", "MSD Y", "MSD Z", "MSD"], multiple_frames = True)
