all : rotate_and_remove calculate_GBE calculate_grain_area calculate_MSD \
csv2tecplot extract_energy find_grains dump2LAMMPS_input dump2tec parse_lammps_output \
xyz2dat generate_impurities
	mv rotate_and_remove calculate_GBE calculate_grain_area calculate_MSD \
	csv2tecplot extract_energy find_grains dump2LAMMPS_input dump2tec parse_lammps_output \
	xyz2dat generate_impurities bin/

debug : rotate_and_remove_dbg calculate_GBE_dbg calculate_grain_area_dbg \
calculate_MSD_dbg csv2tecplot_dbg extract_energy_dbg find_grains_dbg \
dump2LAMMPS_input_dbg dump2tec_dbg parse_lammps_output_dbg xyz2dat_dbg generate_impurities_dbg
	mv *_dbg bin/

rotate_and_remove : rotate_and_remove_atom.cpp atom.o
	g++ -O3 -o rotate_and_remove rotate_and_remove_atom.cpp atom.o

atom.o : atom.h atom.cpp
	g++ -c atom.cpp

calculate_GBE : calculate_GBE.cpp
	g++ -O3 -o calculate_GBE calculate_GBE.cpp

calculate_grain_area : calculate_grain_area.cpp
	g++ -O3 -o calculate_grain_area calculate_grain_area.cpp

calculate_MSD : calculate_MSD.cpp atom.o
	g++ -O3 -o calculate_MSD calculate_MSD.cpp atom.o

csv2tecplot : csv2tecplot.cpp
	g++ -O3 -o csv2tecplot csv2tecplot.cpp

dump2LAMMPS_input : dump2LAMMPS_input.cpp
	g++ -O3 -o dump2LAMMPS_input dump2LAMMPS_input.cpp

dump2tec : dump2tec.cpp
	g++ -O3 -o dump2tec dump2tec.cpp

extract_energy : extract_energy.cpp
	g++ -O3 -o extract_energy extract_energy.cpp

find_grains : find_grains.cpp atom.o
	g++ -O3 -o find_grains find_grains.cpp atom.o

generate_impurities : generate_impurities.cpp atom.cpp atom.o
	g++ -O3 -o generate_impurities generate_impurities.cpp atom.o

parse_lammps_output : parse_lammps_output.cpp
	g++ -O3 -o parse_lammps_output parse_lammps_output.cpp

xyz2dat : xyz2dat.cpp
	g++ -O3 -o xyz2dat xyz2dat.cpp

rotate_and_remove_dbg : rotate_and_remove_atom.cpp atom.o
	g++ -ggdb -g -o rotate_and_remove_dbg rotate_and_remove_atom.cpp atom.o

calculate_GBE_dbg : calculate_GBE.cpp
	g++ -ggdb -g -o calculate_GBE_dbg calculate_GBE.cpp

calculate_grain_area_dbg : calculate_grain_area.cpp
	g++ -ggdb -g -o calculate_grain_area_dbg calculate_grain_area.cpp

calculate_MSD_dbg : calculate_MSD.cpp atom.o
	g++ -ggdb -g -o calculate_MSD_dbg calculate_MSD.cpp atom.o

csv2tecplot_dbg : csv2tecplot.cpp
	g++ -ggdb -g -o csv2tecplot_dbg csv2tecplot.cpp

dump2LAMMPS_input_dbg : dump2LAMMPS_input.cpp
	g++ -ggdb -g -o dump2LAMMPS_input_dbg dump2LAMMPS_input.cpp

dump2tec_dbg : dump2tec.cpp
	g++ -ggdb -g -o dump2tec_dbg dump2tec.cpp

extract_energy_dbg : extract_energy.cpp
	g++ -ggdb -g -o extract_energy_dbg extract_energy.cpp

find_grains_dbg : find_grains.cpp atom.cpp atom.o
	g++ -ggdb -g -o find_grains_dbg find_grains.cpp atom.o

generate_impurities_dbg : generate_impurities.cpp atom.cpp atom.o
	g++ -ggdb -g -o generate_impurities_dbg generate_impurities.cpp atom.o

parse_lammps_output_dbg : parse_lammps_output.cpp
	g++ -ggdb -g -o parse_lammps_output_dbg parse_lammps_output.cpp

xyz2dat_dbg : xyz2dat.cpp
	g++ -ggdb -g -o xyz2dat_dbg xyz2dat.cpp

clean :
	cd bin/ && ls -I '*.*' | xargs rm && cd .. && rm *.o
