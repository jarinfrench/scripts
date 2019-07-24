.SECONDEXPANSION:

UNAME = $(shell uname -s)

ifeq ($(UNAME), Darwin)
	CC = clang++
	CLIB = -stdlib=libc++
	ADDITIONAL_LDFLAGS =
	NODBGFLAGS = -O3 -std=gnu++11 -o
	DBGFLAGS = -ggdb -g -std=gnu++11 -o
endif
ifeq ($(UNAME), Linux)
	CC = g++
	CLIB =
	ADDITIONAL_LDFLAGS =
	NODBGFLAGS = -O3 -std=c++11 -o
	DBGFLAGS = -ggdb -g -std=c++11 -o
endif

BIN =../bin

ATOMPROGS = rotate_and_remove calculate_displacement calculate_MSD  dump2dat \
find_grains track_atoms generate_impurities xyz2dat shift_upper_grain
ATOMDBGPROGS = rotate_and_remove_dbg calculate_displacement_dbg calculate_MSD_dbg \
dump2dat_dbg find_grains_dbg track_atoms_dbg generate_impurities_dbg xyz2dat_dbg \
shift_upper_grain_dbg
NOATOMPROGS = calculate_GBE calculate_grain_area csv2tecplot \
dump2tec extract_energy parse_lammps_output grains2slice csl
NOATOMDBGPROGS = calculate_GBE_dbg calculate_grain_area_dbg csv2tecplot_dbg \
dump2tec_dbg extract_energy_dbg parse_lammps_output_dbg grains2slice_dbg csl_dbg
POSITIONPROGS = phase_field
POSITIONDBGPROGS = phase_field_dbg

ERRORCODES = error_code_defines.h

all: FLAGS = $(CLIB) $(NODBGFLAGS)
all: atom.o position.o atom_builds non_atom_builds position_builds

debug: FLAGS = $(CLIB) $(DBGFLAGS)
debug: atom_builds_dbg non_atom_builds_dbg position_dbg_builds

# The extra line for each rule allows a different value to be passed to the pattern-matched rule.
atom_builds: TARGETS=atom.o
atom_builds: $(addprefix $(BIN)/, $(ATOMPROGS)) atom.o

non_atom_builds: TARGETS=
non_atom_builds: $(addprefix $(BIN)/, $(NOATOMPROGS))

position_builds: TARGETS=position.o
position_builds: $(addprefix $(BIN)/, $(POSITIONPROGS)) position.o

atom_builds_dbg: TARGETS=atom.o
atom_builds_dbg: $(addprefix $(BIN)/, $(ATOMDBGPROGS)) atom.o

non_atom_builds_dbg: TARGETS=
non_atom_builds_dbg: $(addprefix $(BIN)/, $(NOATOMDBGPROGS))

position_dbg_builds: TARGETS=position.o
position_dbg_builds: $(addprefix $(BIN)/, $(POSITIONDBGPROGS)) position.o

atom.o : atom.h atom.cpp
	$(CC) -O3 -std=c++11 $(CLIB) -c atom.cpp $(ADDITIONAL_LDFLAGS)

position.o : position.h position.cpp
	$(CC) -O3 -std=c++11 $(CLIB) -c position.cpp $(ADDITIONAL_LDFLAGS)

# Explanation of this rule:
# the % sign is simply what caused this rule to be used.
# $(@F) gives the filename (which in this case is equal to the %, though I couldn't use that)
# subst _dbg,, : subst is the command to swap strings.  the first argument (_dbg)
# is the string to be replaced, the second argument (between the two commas) is
# the replacement string, and the third argument ($(@F) in this case) is the filename
$(BIN)/% : $$(addsuffix .cpp, $$(subst _dbg,,$$(@F))) $$(TARGETS) $$(ERRORCODES)
	$(CC) $(FLAGS) $(addprefix $(BIN)/, $(@F)) $(addsuffix .cpp, $(subst _dbg,,$(@F))) $(TARGETS) $(ADDITIONAL_LDFLAGS)

# ls -I '*.*' lists all files that do NOT match the pattern '*.*' - all of the
# executables created from this makefile.  BEWARE CHANGING FORMATS FOR OTHER FILE TYPES
clean :
	cd ../bin/ && find . -type f ! -name "*.*" | xargs rm && cd ../src && rm *.o

# This allows variables to be printed: i.e. make print-FLAGS will show the value of FLAGS.
print-%  : ; @echo $* = $($*)

.PHONY : all debug clean
