FLAGS = -O3 -std=c++11 -o
DBGFLAGS = -ggdb -g -std=c++11 -o
BIN=../bin

all : rotate_and_remove calculate_GBE calculate_grain_area calculate_MSD calculate_displacement \
csv2tecplot extract_energy find_grains find_new_positions dump2LAMMPS_input dump2tec parse_lammps_output \
xyz2dat generate_impurities debug
	mv rotate_and_remove calculate_GBE calculate_grain_area calculate_MSD calculate_displacement \
	csv2tecplot extract_energy find_grains dump2LAMMPS_input dump2tec parse_lammps_output \
	xyz2dat generate_impurities find_new_positions ../bin/

debug : rotate_and_remove_dbg calculate_GBE_dbg calculate_grain_area_dbg calculate_displacement_dbg \
calculate_MSD_dbg csv2tecplot_dbg extract_energy_dbg find_grains_dbg find_new_positions_dbg\
dump2LAMMPS_input_dbg dump2tec_dbg parse_lammps_output_dbg xyz2dat_dbg generate_impurities_dbg
	mv *_dbg ../bin/

rotate_and_remove : rotate_and_remove_atom.cpp atom.o error_code_defines.h
	g++ ${FLAGS} rotate_and_remove rotate_and_remove_atom.cpp atom.o

atom.o : atom.h atom.cpp
	g++ -O3 -c atom.cpp

calculate_displacement : calculate_displacement.cpp atom.o
	g++ ${FLAGS} calculate_displacement calculate_displacement.cpp atom.o

calculate_GBE : calculate_GBE.cpp
	g++ ${FLAGS} calculate_GBE calculate_GBE.cpp

calculate_grain_area : calculate_grain_area.cpp
	g++ ${FLAGS} calculate_grain_area calculate_grain_area.cpp

calculate_MSD : calculate_MSD.cpp atom.o
	g++ ${FLAGS} calculate_MSD calculate_MSD.cpp atom.o

csv2tecplot : csv2tecplot.cpp
	g++ ${FLAGS} csv2tecplot csv2tecplot.cpp

dump2LAMMPS_input : dump2LAMMPS_input.cpp
	g++ ${FLAGS} dump2LAMMPS_input dump2LAMMPS_input.cpp

dump2tec : dump2tec.cpp
	g++ ${FLAGS} dump2tec dump2tec.cpp

extract_energy : extract_energy.cpp
	g++ ${FLAGS} extract_energy extract_energy.cpp

find_grains : find_grains.cpp atom.o error_code_defines.h
	g++ ${FLAGS} find_grains find_grains.cpp atom.o

find_new_positions : find_new_positions.cpp atom.o
	g++ ${FLAGS} find_new_positions find_new_positions.cpp atom.o

generate_impurities : generate_impurities.cpp atom.cpp atom.o
	g++ ${FLAGS} generate_impurities generate_impurities.cpp atom.o

parse_lammps_output : parse_lammps_output.cpp
	g++ ${FLAGS} parse_lammps_output parse_lammps_output.cpp

xyz2dat : xyz2dat.cpp atom.o error_code_defines.h
	g++ ${FLAGS} xyz2dat xyz2dat.cpp atom.o

rotate_and_remove_dbg : rotate_and_remove_atom.cpp atom.o error_code_defines.h
	g++ ${DBGFLAGS} rotate_and_remove_dbg rotate_and_remove_atom.cpp atom.o

calculate_displacement_dbg : calculate_displacement.cpp atom.o
	g++ ${DBGFLAGS} calculate_displacement_dbg calculate_displacement.cpp atom.o

calculate_GBE_dbg : calculate_GBE.cpp
	g++ ${DBGFLAGS} calculate_GBE_dbg calculate_GBE.cpp

calculate_grain_area_dbg : calculate_grain_area.cpp
	g++ ${DBGFLAGS} calculate_grain_area_dbg calculate_grain_area.cpp

calculate_MSD_dbg : calculate_MSD.cpp atom.o
	g++ ${DBGFLAGS} calculate_MSD_dbg calculate_MSD.cpp atom.o

csv2tecplot_dbg : csv2tecplot.cpp
	g++ ${DBGFLAGS} csv2tecplot_dbg csv2tecplot.cpp

dump2LAMMPS_input_dbg : dump2LAMMPS_input.cpp
	g++ ${DBGFLAGS} dump2LAMMPS_input_dbg dump2LAMMPS_input.cpp

dump2tec_dbg : dump2tec.cpp
	g++ ${DBGFLAGS} dump2tec_dbg dump2tec.cpp

extract_energy_dbg : extract_energy.cpp
	g++ ${DBGFLAGS} extract_energy_dbg extract_energy.cpp

find_grains_dbg : find_grains.cpp atom.cpp atom.o error_code_defines.h
	g++ ${DBGFLAGS} find_grains_dbg find_grains.cpp atom.o

find_new_positions_dbg : find_new_positions.cpp atom.o
	g++ ${DBGFLAGS} find_new_positions_dbg find_new_positions.cpp atom.o

generate_impurities_dbg : generate_impurities.cpp atom.cpp atom.o
	g++ ${DBGFLAGS} generate_impurities_dbg generate_impurities.cpp atom.o

parse_lammps_output_dbg : parse_lammps_output.cpp
	g++ ${DBGFLAGS} parse_lammps_output_dbg parse_lammps_output.cpp

xyz2dat_dbg : xyz2dat.cpp atom.o error_code_defines.h
	g++ ${DBGFLAGS} xyz2dat_dbg xyz2dat.cpp atom.o

clean :
	cd ../bin/ && ls -I '*.*' | xargs rm && cd ../src && rm *.o

# TARGETDIR = ../bin
# PROGS     = rotate_and_remove calculate_GBE calculate_grain_area calculate_MSD calculate_displacement \
# csv2tecplot extract_energy find_grains find_new_positions dump2LAMMPS_input dump2tec parse_lammps_output \
# xyz2dat generate_impurities
# DBGPROGS  = rotate_and_remove_dbg calculate_GBE_dbg calculate_grain_area_dbg calculate_displacement_dbg \
# calculate_MSD_dbg csv2tecplot_dbg extract_energy_dbg find_grains_dbg find_new_positions_dbg \
# dump2LAMMPS_input_dbg dump2tec_dbg parse_lammps_output_dbg xyz2dat_dbg generate_impurities_dbg
#
# OBJPROG     = $(addprefix $(TARGETDIR)/, $(PROGS))
# OBJDBGPROGS = $(addprefix $(TARGETDIR)/, $(DBGPROGS))
#
# all : $(OBJPROG) $(OBJDBGPROGS)
#
# debug : $(OBJDBGPROGS)
#
# $(TARGETDIR)/rotate_and_remove : rotate_and_remove_atom.cpp atom.o
# 	g++ ${FLAGS} rotate_and_remove rotate_and_remove_atom.cpp atom.o
#
# atom.o : atom.h atom.cpp
# 	g++ -c atom.cpp
#
# $(TARGETDIR)/calculate_displacement : calculate_displacement.cpp atom.o
# 	g++ ${FLAGS} calculate_displacement calculate_displacement.cpp atom.o
#
# $(TARGETDIR)/calculate_GBE : calculate_GBE.cpp
# 	g++ ${FLAGS} calculate_GBE calculate_GBE.cpp
#
# $(TARGETDIR)/calculate_grain_area : calculate_grain_area.cpp
# 	g++ ${FLAGS} calculate_grain_area calculate_grain_area.cpp
#
# $(TARGETDIR)/calculate_MSD : calculate_MSD.cpp atom.o
# 	g++ ${FLAGS} calculate_MSD calculate_MSD.cpp atom.o
#
# $(TARGETDIR)/csv2tecplot : csv2tecplot.cpp
# 	g++ ${FLAGS} csv2tecplot csv2tecplot.cpp
#
# $(TARGETDIR)/dump2LAMMPS_input : dump2LAMMPS_input.cpp
# 	g++ ${FLAGS} dump2LAMMPS_input dump2LAMMPS_input.cpp
#
# $(TARGETDIR)/dump2tec : dump2tec.cpp
# 	g++ ${FLAGS} dump2tec dump2tec.cpp
#
# $(TARGETDIR)/extract_energy : extract_energy.cpp
# 	g++ ${FLAGS} extract_energy extract_energy.cpp
#
# $(TARGETDIR)/find_grains : find_grains.cpp atom.o
# 	g++ ${FLAGS} find_grains find_grains.cpp atom.o
#
# $(TARGETDIR)/find_new_positions : find_new_positions.cpp atom.o
# 	g++ ${FLAGS} find_new_positions find_new_positions.cpp atom.o
#
# $(TARGETDIR)/generate_impurities : generate_impurities.cpp atom.cpp atom.o
# 	g++ ${FLAGS} generate_impurities generate_impurities.cpp atom.o
#
# $(TARGETDIR)/parse_lammps_output : parse_lammps_output.cpp
# 	g++ ${FLAGS} parse_lammps_output parse_lammps_output.cpp
#
# $(TARGETDIR)/xyz2dat : xyz2dat.cpp
# 	g++ ${FLAGS} xyz2dat xyz2dat.cpp
#
# $(TARGETDIR)/rotate_and_remove_dbg : rotate_and_remove_atom.cpp atom.o
# 	g++ ${DBGFLAGS} rotate_and_remove_dbg rotate_and_remove_atom.cpp atom.o
#
# $(TARGETDIR)/calculate_displacement_dbg : calculate_displacement.cpp atom.o
# 	g++ ${DBGFLAGS} calculate_displacement_dbg calculate_displacement.cpp atom.o
#
# $(TARGETDIR)/calculate_GBE_dbg : calculate_GBE.cpp
# 	g++ ${DBGFLAGS} calculate_GBE_dbg calculate_GBE.cpp
#
# $(TARGETDIR)/calculate_grain_area_dbg : calculate_grain_area.cpp
# 	g++ ${DBGFLAGS} calculate_grain_area_dbg calculate_grain_area.cpp
#
# $(TARGETDIR)/calculate_MSD_dbg : calculate_MSD.cpp atom.o
# 	g++ ${DBGFLAGS} calculate_MSD_dbg calculate_MSD.cpp atom.o
#
# $(TARGETDIR)/csv2tecplot_dbg : csv2tecplot.cpp
# 	g++ ${DBGFLAGS} csv2tecplot_dbg csv2tecplot.cpp
#
# $(TARGETDIR)/dump2LAMMPS_input_dbg : dump2LAMMPS_input.cpp
# 	g++ ${DBGFLAGS} dump2LAMMPS_input_dbg dump2LAMMPS_input.cpp
#
# $(TARGETDIR)/dump2tec_dbg : dump2tec.cpp
# 	g++ ${DBGFLAGS} dump2tec_dbg dump2tec.cpp
#
# $(TARGETDIR)/extract_energy_dbg : extract_energy.cpp
# 	g++ ${DBGFLAGS} extract_energy_dbg extract_energy.cpp
#
# $(TARGETDIR)/find_grains_dbg : find_grains.cpp atom.cpp atom.o
# 	g++ ${DBGFLAGS} find_grains_dbg find_grains.cpp atom.o
#
# $(TARGETDIR)/find_new_positions_dbg : find_new_positions.cpp atom.o
# 	g++ ${DBGFLAGS} find_new_positions_dbg find_new_positions.cpp atom.o
#
# $(TARGETDIR)/generate_impurities_dbg : generate_impurities.cpp atom.cpp atom.o
# 	g++ ${DBGFLAGS} generate_impurities_dbg generate_impurities.cpp atom.o
#
# $(TARGETDIR)/parse_lammps_output_dbg : parse_lammps_output.cpp
# 	g++ ${DBGFLAGS} parse_lammps_output_dbg parse_lammps_output.cpp
#
# $(TARGETDIR)/xyz2dat_dbg : xyz2dat.cpp
# 	g++ ${DBGFLAGS} xyz2dat_dbg xyz2dat.cpp
#
# clean :
# 	cd ../bin/ && ls -I '*.*' | xargs rm && cd ../src && rm *.o
