CC = g++
FLAGS = -O3 -std=c++11 -o
DBGFLAGS = -ggdb -g -std=c++11 -o
BIN =../bin

ATOMPROGS = rotate_and_remove calculate_displacement calculate_MSD find_grains \
find_new_positions generate_impurities xyz2dat
ATOMDBGPROGS = rotate_and_remove_dbg calculate_displacement_dbg calculate_MSD_dbg \
find_grains_dbg find_new_positions_dbg generate_impurities_dbg xyz2dat_dbg
NOATOMPROGS = calculate_GBE calculate_grain_area csv2tecplot dump2LAMMPS_input \
dump2tec extract_energy parse_lammps_output
NOATOMDBGPROGS = calculate_GBE_dbg calculate_grain_area_dbg csv2tecplot_dbg \
dump2LAMMPS_input_dbg dump2tec_dbg extract_energy_dbg parse_lammps_output_dbg

ERRORCODES = error_code_defines.h

all: atom_builds non_atom_builds debug

debug: atom_builds_dbg non_atom_builds_dbg

atom_builds: TARGETS=atom.o
atom_builds: $(addprefix $(BIN)/, $(ATOMPROGS)) atom.o

non_atom_builds: TARGETS=
non_atom_builds: $(addprefix $(BIN)/, $(NOATOMPROGS))

atom_builds_dbg: TARGETS=atom.o
atom_builds_dbg: $(addprefix $(BIN)/, $(ATOMDBGPROGS)) atom.o

non_atom_builds_dbg: TARGETS=
non_atom_builds_dbg: $(addprefix $(BIN)/, $(NOATOMDBGPROGS))

atom.o : atom.h atom.cpp
	$(CC) -O3 -std=c++11 -c atom.cpp

$(BIN)/% : $(addsuffix .cpp, $(@F)) $(TARGETS) $(ERRORCODES)
	$(CC) $(FLAGS) $(addprefix $(BIN)/, $(@F)) $(addsuffix .cpp, $(@F)) $(TARGETS)

clean :
	cd ../bin/ && ls -I '*.*' | xargs rm && cd ../src && rm *.o
