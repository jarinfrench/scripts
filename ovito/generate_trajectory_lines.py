#! /usr/bin/env python3.9

from ovito.io import import_file
import ovito.modifiers as om
from ovito.vis import Viewport, ColorLegendOverlay
from ovito.qt_compat import QtCore

# Set up the viewport
vp = Viewport()
vp.type = Viewport.Type.Top
vp.zoom_all()

pipeline = import_file("*.dump")  # loads all dump files in the current directory
data = pipeline.compute()  # use the last frame
num_types = max(data.particles["Particle Type"].array)
print(f"Found {pipeline.source.num_frames} data files")

# set up the modifiers
# These four lines set up a SelectTypeModifer for selecting individual atom types
select_types_modifiers = []
for i in range(1, num_types + 1):
    select_types_modifiers.append(om.SelectTypeModifier(types={i}))

# This modifier will be used to generate trajectory lines for the selected atoms
trajectory_modifier = om.GenerateTrajectoryLinesModifier(only_selected=True)

for i in range(1, num_types):
    pipeline.modifiers.append(select_types_modifiers[i])  # selects type
    pipeline.modifiers.append(
        trajectory_modifier
    )  # generate trajectory lines (must generate() after adding to the pipeline)
    data = pipeline.compute()
    trajectory_modifier.generate()
    trajectory_modifier.vis.color_mapping_property = "Time"
    trajectory_modifier.vis.color_mapping_gradient = om.ColorCodingModifier.Rainbow()
    trajectory_modifier.vis.color_mapping_interval = (0, pipeline.source.num_frames)
    trajectory_modifier.vis.width = 1.0

    pipeline.add_to_scene()
    vis_particles = pipeline.source.data.particles.vis
    vis_particles.enabled = False
    vp.overlays.append(
        ColorLegendOverlay(
            color_mapping_source=trajectory_modifier.vis,
            title="Frame",
            alignment=QtCore.Qt.AlignHCenter ^ QtCore.Qt.AlignBottom,
            orientation=QtCore.Qt.Vertical,
            offset_x=0.42,
            offset_y=0.35,
        )
    )
    vp.zoom_all()
    vp.render_image(filename=f"trajectories_type{i + 1}.png")
    pipeline.modifiers.clear()
