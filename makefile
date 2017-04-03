all : rotate_and_remove calculate_GBE check_distances extract_energy parse_lammps xyz2dat

rotate_and_remove : rotate_and_remove_atom.cpp atom.o
	g++ -O3 -o rotate_and_remove rotate_and_remove_atom.cpp atom.o

atom.o : atom.h atom.cpp
	g++ -c atom.cpp

calculate_GBE : calculate_GBE.cpp
	g++ -O3 -o calculate_GBE calculate_GBE.cpp

check_distances : check_distances.cpp
	g++ -O3 -o checked_distances check_distances.cpp

extract_energy : extract_energy.cpp
	g++ -O3 -o extract_energy extract_energy.cpp

parse_lammps : parse_lammps.cpp
	g++ -O3 -o parse_lammps parse_lammps.cpp

xyz2dat : xyz2dat.cpp
	g++ -O3 -o xyz2dat xyz2dat.cpp
