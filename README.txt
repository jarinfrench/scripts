This directory contains a set of scripts that have been developed by Jarin French
for furthering research goals.

angles_generator.py: A simple script that generates a list of angles for use
in other scripts (i.e. genOrientationMatrices.sh).

atom.*: C++ class that handles atoms as output by LAMMPS dump files, etc.  Each
atom is assigned an ID number, a type, a charge (can be neutral), and a
position.  Code processing the atoms have access to an additional parameter called
"mark" which is an integer useful for defining various properties.  Note that
position is taken care of explicitly, but there is a class written that could
separate out the position into its own class.

calculate_GBE.cpp: From the file header: "This script reads a file as generated
by extract_all_energies.sh, and calculates the per atom energy at the grain
boundary.  The filename can be passed in via command line.  Note that this
program expects a filename in the format of <description>_total_energies.csv
IMPORTANT PARAMETERS
  grain radius: the radius of the grain that was specified when rotating the grain.
    note that this assumes a cylindrical grain.  Given in Angstroms.
  grain thickness: the thickness of the grain (determined by Lz).  This is given
    in Angstroms."
In other words, this script looks at the energy output from another script, and
calculates a grain boundary energy.  Note that the script expects to find the
string "_total" in the filename given to it.

check_all_distances.sh: Simple script that checks the distances between atoms
as given by LAMMPS .dat files.  Used in the development of rotate_and_remove_atoms.cpp

check_distances.cpp: C++ script used by check_all_distances.sh.

convertOmatToEuler.py: Script developed at INL to convert an orientation matrix
to a set of Euler angles.  Not fully developed - still some bugs to work out.

csv2tecplot.cpp: Script to easily convert a csv file to a file readable by
tecplot.  As tecplot is basically csv without the commas, it's an easy adjustment.
Note that the csv files that are processed are assumed to have four lines of a
header.

Euler_angle_generator.py: Another script developed at INL based off of work done
by John-Michael Bradley.  This script generates a set of Euler angles based on
input of a maximum angle, an axis, and the type of grain boundary (twist or tilt).
Still incomplete, but should handle specific cases well.

extract_all_energies.sh: Uses the script extract_energy.cpp to parse LAMMPS output
files for the minimum energy value of a minimized atomic configuration. Only
parses the directory in which the script is called.

extract_energy.cpp: A script that parses a LAMMPS output file for the minimum
energy from a minimizing simulation.

find_grains.cpp: A script that parses either a LAMMPS dump file, or a LAMMPS
input file, and assigns atoms to a grain based on a symmetry parameter.

gen_*_grains.sh: These scripts (Cu and UO2) generate the rotated structures for
cylindrical grain boundaries.  As the script that generates the rotated structures
creates 3 data files, the different files are moved to specific directories.

gen_pbs_script(_Cu).sh: Generates the input files for UO2 (copper) to be used in
a LAMMPS simulation.  Requires a base file for the LAMMPS input file with a very
specific format, but this can be changed.

genOrientationMatrices.sh: A script developed at INL to generate the orientation
matrices for angles as specified by a csv file containing the angle and energy
associated with it.  Note that the energy is not used in any calculations, so
it is a meaningless number in this context.

grain_finding.sh: A script that parses the current directory and finds the files
that end in .dump for processing to find the grains (uses the script find_grains.cpp).

makefile: Compiles the various C++ scripts.

orientation_matrix.py: This script was also developed at INL to create the P and
Q matrices (see Bulatov et al. Acta Materialia 65 (2014) 161-175) used in the
GB5DOF Matlab script.  Still some error to work out, but does work for the most
part.

parse_all_dumps.sh: A simple script to parse all of the LAMMPS dump files in the
current directory into a format that LAMMPS can read as input.

parse_lammps_dump.cpp: A script utilized by parse_all_dumps.sh to convert a dump
file from LAMMPS into a LAMMPS input file.

parse_lammps_output.cpp: A script that converts output from LAMMPS to a csv file.

parseEuler.py: A script developed at INL to parse the orientation matrices database
(as generated by orientation_matrix.py) to grab the Euler angles.

planeNormal.py: Simple script to generate images representing boundary planes with
their unit normal vectors.

plot_data.py: A script that reads a csv file (or series of csv files) and plots
the data on one plot.

plot_LAMMPS_data.py: From the file: "This script is should be utilized AFTER
using the parse_lammps_output.cpp script.  This script reads in the data from
that file, and then prompts the user as to which plots to show."

position.*: C++ class that contains positions (x, y, and z).  Not used in anything.

rotate_and_remove_atom.cpp: Script that rotates atoms within a radius r of
the center of the structure, and then removes atoms that are too close to each
other.

rotation_matrix.py: Generates rotation matrices to for use in Bulatov's GB5DOF
MATLAB script.  Required if the P and Q matrices from orientation_matrix.py are
used.

RSW.py: Generates a graphic that describes what an RSW function is.

uo2_potential.py: Plots the Basak potential with it's individual parts for UO2.

xyz2dat.cpp: A script that converts data from GBstudio (or a .xyz file) to a
LAMMPS input file (.dat)
