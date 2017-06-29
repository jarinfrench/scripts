all : rotate_and_remove calculate_GBE calculate_mobility check_distances \
csv2tecplot extract_energy find_grains parse_lammps_dump parse_lammps_output \
xyz2dat

debug : rotate_and_remove_dbg calculate_GBE_dbg calculate_mobility_dbg \
check_distances_dbg csv2tecplot_dbg extract_energy_dbg find_grains_dbg \
parse_lammps_dump_dbg parse_lammps_output_dbg xyz2dat_dbg

rotate_and_remove : rotate_and_remove_atom.cpp atom.o
	g++ -O3 -o rotate_and_remove rotate_and_remove_atom.cpp atom.o

atom.o : atom.h atom.cpp
	g++ -c atom.cpp

calculate_GBE : calculate_GBE.cpp
	g++ -O3 -o calculate_GBE calculate_GBE.cpp

calculate_mobility : calculate_mobility.cpp
	g++ -O3 -o calculate_mobility calculate_mobility.cpp

check_distances : check_distances.cpp atom.cpp atom.o
	g++ -O3 -o check_distances check_distances.cpp atom.o

csv2tecplot : csv2tecplot.cpp
	g++ -O3 -o csv2tecplot csv2tecplot.cpp

extract_energy : extract_energy.cpp
	g++ -O3 -o extract_energy extract_energy.cpp

find_grains : find_grains.cpp atom.cpp atom.o
	g++ -O3 -o find_grains find_grains.cpp atom.o

parse_lammps_dump : parse_lammps_dump.cpp
	g++ -O3 -o parse_lammps_dump parse_lammps_dump.cpp

parse_lammps_output : parse_lammps_output.cpp
	g++ -O3 -o parse_lammps_output parse_lammps_output.cpp

xyz2dat : xyz2dat.cpp
	g++ -O3 -o xyz2dat xyz2dat.cpp

rotate_and_remove_dbg : rotate_and_remove_atom.cpp atom.o
	g++ -ggdb -o rotate_and_remove_dbg rotate_and_remove_atom.cpp atom.o

calculate_GBE_dbg : calculate_GBE.cpp
	g++ -ggdb -o calculate_GBE_dbg calculate_GBE.cpp

calculate_mobility_dbg : calculate_mobility.cpp
	g++ -ggdb -o calculate_mobility_dbg calculate_mobility.cpp

check_distances_dbg : check_distances.cpp atom.cpp atom.o
	g++ -ggdb -o check_distances_dbg check_distances.cpp atom.o

csv2tecplot_dbg : csv2tecplot.cpp
	g++ -ggdb -o csv2tecplot_dbg csv2tecplot.cpp

extract_energy_dbg : extract_energy.cpp
	g++ -ggdb -o extract_energy_dbg extract_energy.cpp

find_grains_dbg : find_grains.cpp atom.cpp atom.o
	g++ -ggdb -o find_grains_dbg find_grains.cpp atom.o

parse_lammps_dump_dbg : parse_lammps_dump.cpp
	g++ -ggdb -o parse_lammps_dump_dbg parse_lammps_dump.cpp

parse_lammps_output_dbg : parse_lammps_output.cpp
	g++ -ggdb -o parse_lammps_output_dbg parse_lammps_output.cpp

xyz2dat_dbg : xyz2dat.cpp
	g++ -ggdb -o xyz2dat_dbg xyz2dat.cpp

clean :
	rm rotate_and_remove calculate_GBE calculate_mobility check_distances \
	csv2tecplot extract_energy find_grains parse_lammps_dump parse_lammps_output \
	xyz2dat *_dbg
