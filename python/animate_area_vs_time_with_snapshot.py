import argparse
from os import mkdir
from os.path import exists, splitext

import cv2
import pandas as pd
import plotly.graph_objects as go
from ovito.io import import_file
from ovito.modifiers import DeleteSelectedModifier, SelectTypeModifier
from ovito.vis import Viewport

parser = argparse.ArgumentParser(
    usage="%(prog)s [-h] file [options]",
    description="Create animation of specified area vs time data.",
)
parser.add_argument("file", help="Datafile to plot.")
parser.add_argument("-o", "--output", default="*.gif", help="Output file name")
parser.add_argument("-i", "--ignore", nargs="+", type=int, help="Atom types to ignore")
parser.add_argument(
    "--dt", type=int, help="The number of ps between snapshots", default=20
)

args = parser.parse_args()

if args.output == "*.gif":
    args.output = splitext(args.file)[0] + ".gif"

if not exists(".tmp"):
    mkdir(".tmp")

# For the area_vs_time data
df = pd.read_csv(
    args.file,
    delimiter=" ",
    skiprows=1,
    names=("time", "area1", "area2"),
    dtype={"time": "int", "area1": "float", "area2": "float"},
)


for i in range(len(df)):
    data = df[df["time"] < (i - 1) * args.dt]
    newpoint = df[df["time"] == (i - 1) * args.dt]
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=data["time"], y=data["area2"], mode="markers"))
    fig.add_trace(
        go.Scatter(
            x=newpoint["time"],
            y=newpoint["area2"],
            mode="markers",
            marker={"color": "red", "size": 15},
        )
    )
    fig.update_yaxes(range=[-100, 1.1 * max(df["area2"])])
    fig.update_xaxes(range=[-10, 1.1 * max(df["time"])])
    fig.update_layout(width=1000, height=1000, showlegend=False)
    fig.write_image(f".tmp/area_data_{i:04}.png")


# For the snapshot animation
pipeline = import_file("*.dump")
data = pipeline.compute()
if len(set(data.particles.particle_type)) == 2:
    colors = [
        (85 / 255, 87 / 255, 83 / 255),
        (239 / 255, 41 / 255, 41 / 255),
    ]  # gray and red
    radii = [0.6, 1.5]
elif len(set(data.particles.particle_type)) == 3:
    colors = [
        (85 / 255, 87 / 255, 83 / 255),
        (0 / 255, 164 / 255, 123 / 255),
        (239 / 255, 41 / 255, 41 / 255),
    ]  # gray, green, and red
    radii = [0.6, 1.0, 1.5]
else:
    "Incorrect number of particle types"
    exit()


def particle_setup(frame, data):
    types = data.particles_.particle_types_
    for i in set(types):
        types.type_by_id_(i).radius = radii[i - 1]
        types.type_by_id_(i).color = colors[i - 1]


if args.ignore:
    pipeline.modifiers.append(
        SelectTypeModifier(property="Particle Type", types=set(args.ignore))
    )
    pipeline.modifiers.append(DeleteSelectedModifier())
pipeline.modifiers.append(particle_setup)
# particle_colors = [None] * data.particles.count
# particle_radius = [None] * data.particles.count
# for i in range(data.particles.count):
#     particle_colors[i] = colors[data.particles.particle_types[i] - 1]
#     particle_radius[i] = radii[data.particles.particle_types[i] - 1]
# data.particles_.create_property("Color", data=particle_colors)
# data.particles_.create_property("Radius", data=particle_radius)

pipeline.add_to_scene()
vp = Viewport()
vp.type = Viewport.Type.Top
vp.zoom_all()
vp.render_anim(
    size=(1000, 1000),
    filename=".tmp/system_video_.png",
)

for i in range(len(df)):
    img1 = cv2.imread(f".tmp/area_data_{i:04}.png")
    img2 = cv2.imread(f".tmp/system_video_{i:04}.png")
    im_h = cv2.hconcat([img1, img2])
    cv2.imwrite(f".tmp/combined.{i:04}.png", im_h)
