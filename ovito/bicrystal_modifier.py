#! /usr/bin/env python3.8
"""This script is meant to be used as an Ovito modifier to assign
atoms to one of two grains (or a grain boundary region) given the
(initial) orientation of the two grains, and a cutoff parameter"""

from ovito.data import DataCollection, CutoffNeighborFinder
from ovito.io import import_file, export_file
import numpy as np
import argparse
from tqdm import tqdm
import concurrent.futures as cf


# adapted from LAMMPS via fix_orient_eco.cpp
# commit 27da7168525980caeeb52496273817e24ceb53a6
def get_reciprocal_vectors(dir_vec: list) -> list:
    # volume of unit cell A
    vol = (
        0.5
        * (
            dir_vec[0][0]
            * (dir_vec[1][1] * dir_vec[2][2] - dir_vec[2][1] * dir_vec[1][2])
            + dir_vec[1][0]
            * (dir_vec[2][1] * dir_vec[0][2] - dir_vec[0][1] * dir_vec[2][2])
            + dir_vec[2][0]
            * (dir_vec[0][1] * dir_vec[1][2] - dir_vec[1][1] * dir_vec[0][2])
        )
        / np.pi
    )
    i_vol = 1.0 / vol
    reciprocal_vectors = [
        [[None for i in range(3)] for j in range(3)] for k in range(2)
    ]

    # grain A: reciprocal_vectors 0
    reciprocal_vectors[0][0][0] = (
        dir_vec[1][1] * dir_vec[2][2] - dir_vec[2][1] * dir_vec[1][2]
    ) * i_vol
    reciprocal_vectors[0][0][1] = (
        dir_vec[1][2] * dir_vec[2][0] - dir_vec[2][2] * dir_vec[1][0]
    ) * i_vol
    reciprocal_vectors[0][0][2] = (
        dir_vec[1][0] * dir_vec[2][1] - dir_vec[2][0] * dir_vec[1][1]
    ) * i_vol

    # grain A: reciprocal_vectors 1
    reciprocal_vectors[0][1][0] = (
        dir_vec[2][1] * dir_vec[0][2] - dir_vec[0][1] * dir_vec[2][2]
    ) * i_vol
    reciprocal_vectors[0][1][1] = (
        dir_vec[2][2] * dir_vec[0][0] - dir_vec[0][2] * dir_vec[2][0]
    ) * i_vol
    reciprocal_vectors[0][1][2] = (
        dir_vec[2][0] * dir_vec[0][1] - dir_vec[0][0] * dir_vec[2][1]
    ) * i_vol

    # grain A: reciprocal_vectors 2
    reciprocal_vectors[0][2][0] = (
        dir_vec[0][1] * dir_vec[1][2] - dir_vec[1][1] * dir_vec[0][2]
    ) * i_vol
    reciprocal_vectors[0][2][1] = (
        dir_vec[0][2] * dir_vec[1][0] - dir_vec[1][2] * dir_vec[0][0]
    ) * i_vol
    reciprocal_vectors[0][2][2] = (
        dir_vec[0][0] * dir_vec[1][1] - dir_vec[1][0] * dir_vec[0][1]
    ) * i_vol

    # volume of unit cell B
    vol = (
        0.5
        * (
            dir_vec[3][0]
            * (dir_vec[4][1] * dir_vec[5][2] - dir_vec[5][1] * dir_vec[4][2])
            + dir_vec[4][0]
            * (dir_vec[5][1] * dir_vec[3][2] - dir_vec[3][1] * dir_vec[5][2])
            + dir_vec[5][0]
            * (dir_vec[3][1] * dir_vec[4][2] - dir_vec[4][1] * dir_vec[3][2])
        )
        / np.pi
    )
    i_vol = 1.0 / vol

    # grain B: reciprocal_vectors 0
    reciprocal_vectors[1][0][0] = (
        dir_vec[4][1] * dir_vec[5][2] - dir_vec[5][1] * dir_vec[4][2]
    ) * i_vol
    reciprocal_vectors[1][0][1] = (
        dir_vec[4][2] * dir_vec[5][0] - dir_vec[5][2] * dir_vec[4][0]
    ) * i_vol
    reciprocal_vectors[1][0][2] = (
        dir_vec[4][0] * dir_vec[5][1] - dir_vec[5][0] * dir_vec[4][1]
    ) * i_vol

    # grain B: reciprocal_vectors 1
    reciprocal_vectors[1][1][0] = (
        dir_vec[5][1] * dir_vec[3][2] - dir_vec[3][1] * dir_vec[5][2]
    ) * i_vol
    reciprocal_vectors[1][1][1] = (
        dir_vec[5][2] * dir_vec[3][0] - dir_vec[3][2] * dir_vec[5][0]
    ) * i_vol
    reciprocal_vectors[1][1][2] = (
        dir_vec[5][0] * dir_vec[3][1] - dir_vec[3][0] * dir_vec[5][1]
    ) * i_vol

    # grain B: reciprocal_vectors 2
    reciprocal_vectors[1][2][0] = (
        dir_vec[3][1] * dir_vec[4][2] - dir_vec[4][1] * dir_vec[3][2]
    ) * i_vol
    reciprocal_vectors[1][2][1] = (
        dir_vec[3][2] * dir_vec[4][0] - dir_vec[4][2] * dir_vec[3][0]
    ) * i_vol
    reciprocal_vectors[1][2][2] = (
        dir_vec[3][0] * dir_vec[4][1] - dir_vec[4][0] * dir_vec[3][1]
    ) * i_vol
    return reciprocal_vectors


def get_norm(
    dir_vec: list,
    reciprocal_vectors: list,
    squared_cutoff: float,
    inv_squared_cutoff: float,
) -> int:
    max_co = 4  # will produce wrong results for rcut > 3 * lattice constant
    neigh = 0
    delta = [None, None, None]  # relative position
    wsum = 0.0  # sum of all weight functions
    reesum = [0.0, 0.0, 0.0]  # sum of real part
    imesum = [0.0, 0.0, 0.0]  # sum of imaginary part
    for i in range(-max_co, max_co + 1):
        for j in range(-max_co, max_co + 1):
            for k in range(-max_co, max_co + 1):
                for ll in range(3):
                    delta[ll] = (
                        dir_vec[0][ll] * i + dir_vec[1][ll] * j + dir_vec[2][ll] * k
                    )
                squared_distance = np.dot(delta, delta)

                # Check if atom is within cutoff region
                if not squared_distance == 0 and squared_distance < squared_cutoff:
                    neigh += 1
                    squared_distance *= inv_squared_cutoff

                    weight = squared_distance * (squared_distance - 2.0) + 1.0
                    wsum += weight

                    # three reciprocal directions
                    for m in range(3):
                        scalar_product = np.dot(reciprocal_vectors[1][m], delta)
                        reesum[m] += weight * np.cos(scalar_product)
                        imesum[m] -= weight * np.sin(scalar_product)
    norm_fac = 3.0 * wsum * wsum
    for k in range(3):
        norm_fac -= reesum[k] * reesum[k] + imesum[k] * imesum[k]
    return (neigh, norm_fac)


def calculateOrderParameter(frame: int, data: DataCollection):
    def getOrderParam(i: int):
        chi = 0.0
        real_phi = [
            [[0.0 for a in range(3)] for b in range(2)]
            for c in range(data.particles.count)
        ]
        imag_phi = [
            [[0.0 for a in range(3)] for b in range(2)]
            for c in range(data.particles.count)
        ]

        for neigh in neigh_finder.find(i):
            squared_distance = neigh.distance_squared
            if squared_distance < squared_cutoff:
                squared_distance *= inv_squared_cutoff
                weight = squared_distance * (squared_distance - 2.0) + 1.0

                for lam in range(2):
                    for k in range(3):
                        scalar_product = np.dot(reciprocal_vectors[lam][k], neigh.delta)
                        real_phi[i][lam][k] += weight * np.cos(scalar_product)
                        imag_phi[i][lam][k] += weight * np.sin(scalar_product)

        for k in range(3):
            chi += (
                real_phi[i][0][k] * real_phi[i][0][k]
                + imag_phi[i][0][k] * imag_phi[i][0][k]
                - real_phi[i][1][k] * real_phi[i][1][k]
                - imag_phi[i][1][k] * imag_phi[i][1][k]
            )
        chi *= inv_norm_fac
        order_params[i][0] = chi

        if chi > args.eta:
            order_params[i][1] = 1
        elif chi < -args.eta:
            order_params[i][1] = 2
        else:
            omega = np.pi / 2.0 * inv_eta * chi
            order_params[i][1] = np.sin(omega)

    yield "Calculating order parameters"

    # Create output particle property
    # first value is the order parameter, second value is the grain number
    # (1 is chi > eta, -1 is chi < eta, other values are GB)
    order_params = data.particles_.create_property(
        "Order Parameter", dtype=float, components=2
    )

    # Adapted from LAMMPS via fix_orient_eco.cpp::post_force
    # commit: 27da7168525980caeeb52496273817e24ceb53a6
    neigh_finder = CutoffNeighborFinder(args.rcut, data)
    with order_params:
        ll = len(range(data.particles.count))
        with tqdm(total=ll) as pbar:
            with cf.ThreadPoolExecutor(max_workers=10) as executor:
                futures = {
                    executor.submit(getOrderParam, item): item
                    for item in range(data.particles.count)
                }
                for future in cf.as_completed(futures):
                    # arg = futures[future]
                    pbar.update(1)


if __name__ == "__main__":
    # Needs the following parameters:
    # Basis vectors of grain 1 (x, y, and z direction vectors = 9 parameters)
    # Basis vectors of grain 2 (x, y and z direction vectors = 9 parameters)
    # Note: orientation matrices uniquely determined by 2 vectors, save effort here?
    # cutoff distance (1 parameter)

    parser = argparse.ArgumentParser(
        usage="$(prog)s [-h] b1 b2 rcut eta",
        description="Assign atoms to one of two grains (or a GB region) based on the "
        "initial orientations of the grains",
    )
    parser.add_argument("b1", nargs=9, type=float, help="basis vectors of grain 1")
    parser.add_argument("b2", nargs=9, type=float, help="basis vectors of grain 2")
    parser.add_argument("rcut", type=float, help="Cutoff parameter")
    parser.add_argument(
        "eta", type=float, help="Truncating parameter for ignoring thermal fluctuations"
    )

    args = parser.parse_args()
    args.b1 = [[i for idx, i in enumerate(args.b1) if idx % 3 == a] for a in range(3)]
    args.b2 = [[i for idx, i in enumerate(args.b2) if idx % 3 == a] for a in range(3)]
    dir_vec = args.b1 + args.b2
    reciprocal_vectors = get_reciprocal_vectors(dir_vec)
    squared_cutoff = args.rcut * args.rcut
    inv_squared_cutoff = 1 / squared_cutoff
    inv_eta = 1.0 / args.eta
    neigh, norm_factor = get_norm(
        dir_vec, reciprocal_vectors, squared_cutoff, inv_squared_cutoff
    )
    inv_norm_fac = 1.0 / norm_factor

    pipeline = import_file("*.dump")
    pipeline.modifiers.append(calculateOrderParameter)
    export_file(
        pipeline,
        "order_params_*.txt",
        "xyz",
        columns=[
            "Particle Identifier",
            "Particle Type",
            "Position.X",
            "Position.Y",
            "Position.Z",
            "Order Parameter",
        ],
        multiple_frames=True,
    )
