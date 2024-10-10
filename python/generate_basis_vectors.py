#! /usr/bin/env python3

import argparse
import numpy as np
import pprint as pp

parser = argparse.ArgumentParser(
    usage="%(prog)s [-h] orientation misorientation a0",
    description="Generate the orientation matrix and basis vectors for the bicrystal "
    "systems generated for cylindrical grain growth",
)

parser.add_argument(
    "orientation",
    nargs=9,
    type=float,
    help="Alignment of the crystal system (z axis is the rotation axis)",
)
parser.add_argument(
    "misorientation", type=float, help="Misorientation angle (in degrees)"
)
parser.add_argument("a0", type=float, help="Lattice parameter (in Angstroms)")
parser.add_argument(
    "structure", choices=["fcc", "bcc"], help="Crystal structure of the system"
)
parser.add_argument(
    "--rotation-matrix",
    action="store_true",
    help="Show the rotation matrices of the grains",
)
parser.add_argument("--pretty", action="store_true", help="Pretty print the values")

args = parser.parse_args()
args.misorientation = np.pi / 180.0 * args.misorientation
cosin = np.cos(args.misorientation)
sinin = np.sin(args.misorientation)
eps = 1e-8  # rounding limit

if args.structure == "fcc":
    basis_vectors = (
        np.array(
            [
                np.array([0.5, 0.5, 0.0]),
                np.array([0.5, 0.0, 0.5]),
                np.array([0.0, 0.5, 0.5]),
            ]
        )
        * args.a0
    )
else:  # bcc
    basis_vectors = (
        np.array(
            [
                np.array([-0.5, 0.5, 0.5]),
                np.array([0.5, -0.5, 0.5]),
                np.array([0.5, 0.5, -0.5]),
            ]
        )
        * args.a0
    )

rotation_matrix = np.array(
    [
        np.array([cosin, -sinin, 0]),
        np.array([sinin, cosin, 0]),
        np.array([0.0, 0.0, 1.0]),
    ]
)

matrix_grain_orientation = [
    [i for idx, i in enumerate(args.orientation) if idx % 3 == a] for a in range(3)
]
matrix_grain_orientation = [i / np.linalg.norm(i) for i in matrix_grain_orientation]

embedded_grain_orientation = np.dot(rotation_matrix, matrix_grain_orientation)

matrix_grain_basis = np.dot(matrix_grain_orientation, basis_vectors).T
embedded_grain_basis = np.dot(embedded_grain_orientation, basis_vectors).T

if args.rotation_matrix:
    print("Matrix grain orientation matrix:")
    if args.pretty:
        pp.pprint(matrix_grain_orientation)
    else:
        print(
            " ".join(
                [
                    str(i) if abs(i) > eps else "0"
                    for i in matrix_grain_orientation.flatten()
                ]
            )
        )

    print("Embedded grain orientation matrix:")
    if args.pretty:
        pp.pprint(embedded_grain_orientation)
    else:
        print(
            " ".join(
                [
                    str(i) if abs(i) > eps else "0"
                    for i in embedded_grain_orientation.flatten()
                ]
            )
        )

print("Matrix grain basis vectors: ")
if args.pretty:
    pp.pprint(matrix_grain_basis)
else:
    print(
        " ".join(
            ([str(i) if abs(i) > eps else "0" for i in matrix_grain_basis.flatten()])
        )
    )

print("Embedded grain basis vectors: ")
if args.pretty:
    pp.pprint(embedded_grain_basis)
else:
    print(
        " ".join(
            ([str(i) if abs(i) > eps else "0" for i in embedded_grain_basis.flatten()])
        )
    )
