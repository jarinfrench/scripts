#! /usr/bin/env python

from ovito.io import import_file, export_file
import ovito.modifiers as om
import argparse

parser = argparse.ArgumentParser(
    usage="%(prog)s [-h] r_in r_out",
    description="Calculate solute concentration in a ring over the whole simulation.",
)

parser.add_argument("r_in", type=float, help="Inner radius of cylinder")
parser.add_argument("r_out", type=float, help="Outer radius of cylinder")
parser.add_argument(
    "--center",
    nargs=3,
    type=float,
    help="Center of the grain. Default is center of simulation cell",
)
args = parser.parse_args()
if not args.center:
    args.center = ["CellSize.X/2", "CellSize.Y/2", "CellSize.Z/2"]

pipeline = import_file("*.dump")
ntypes = len(set(pipeline.compute().particles.particle_type))
if ntypes < 2:
    print("No solute atom types found (ntypes < 2)")
    exit()

# Get the number of atoms in the ring
pipeline.modifiers.append(
    om.ExpressionSelectionModifier(
        expression=f"(Position.X-{args.center[0]})^2+(Position.Y-{args.center[1]})^2"
        f" > {args.r_in}^2 && "
        f"(Position.X-{args.center[0]})^2+(Position.Y-{args.center[1]})^2"
        f" < {args.r_out}^2"
    )
)
# Get the number of solute atoms in the ring
pipeline.modifiers.append(
    om.ExpressionSelectionModifier(
        expression=f"(Position.X-CellSize.X/2)^2+(Position.Y-CellSize.Y/2)^2"
        f" > {args.r_in}^2 && (Position.X -CellSize.X/2)^2+(Position.Y-CellSize.Y/2)^2"
        f" < {args.r_out}^2 && ParticleType==2"
    )
)
# Get the time series evolution of these values
pipeline.modifiers.append(
    om.TimeSeriesModifier(
        operate_on=("SelectExpression.num_selected", "SelectExpression.num_selected.2")
    )
)

data = pipeline.compute()
table = data.make_mutable(data.tables["time-series"])
series = table.y
table.y = table.create_property(
    "solute_percent", data=series[:, 1] / series[:, 0] * 100
)
export_file(table, "concentration_data.txt", "txt/table")
