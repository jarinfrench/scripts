.SECONDEXPANSION:

UNAME = $(shell uname -s)

ifeq ($(UNAME), Darwin)
	CC = clang++
	CLIB = -stdlib=libc++
	STDLIB = -std=gnu++11
	ADDITIONAL_LDFLAGS =
	NODBGFLAGS = -O3 $(STDLIB) -o
	DBGFLAGS = -ggdb -g $(STDLIB) -o
	GREP = ggrep
endif
ifeq ($(UNAME), Linux)
	CC = g++
	CLIB =
	STDLIB = -std=c++11
	ADDITIONAL_LDFLAGS =
	NODBGFLAGS = -O3 $(STDLIB) -o
	DBGFLAGS = -ggdb -g $(STDLIB) -D_GLIBCXX_DEBUG -o
	GREP = grep
	PYTHON33CFLAGS = -I/usr/include/python3.5m -I/usr/include/python3.5m  -Wno-unused-result -Wsign-compare -g -fstack-protector-strong -Wformat -Werror=format-security  -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes
	PYTHON33LDFLAGS = -L/usr/lib/python3.5/config-3.5m-x86_64-linux-gnu -L/usr/lib -lpython3.5m -lpthread -ldl  -lutil -lm  -Xlinker -export-dynamic -Wl,-O1 -Wl,-Bsymbolic-functions
endif

BIN =../bin

# This may need some fine tuning between Mac and Linux
# -L prints the name of each input file that does NOT have a match
# -r searches the search directory recursively
# -Z outputs a zero byte (NULL character) instead of the character that normally follows the file name
# -l prints the name of each input file that does have a match
OMPFILES := $(shell $(GREP) -lr "\#include <omp.h>" *.cpp)
ATOMFILES := $(shell $(GREP) -lr --exclude="atom.cpp" --exclude="verifyNewFile.cpp" "\#include \"atom.h\"" *.cpp)
UTILFILES := $(shell $(GREP) -LrZ "\#include \"atom.h\"" --exclude="atom.cpp" --exclude="position.cpp" --exclude="verifyNewFile.cpp" --exclude="phase_field.cpp" --exclude="Box.cpp" --exclude="LAMMPSDataFileReader.cpp" *.cpp | xargs -0 $(GREP) -Lr "\#include \"position.h\"")
POSITIONFILES := $(shell $(GREP) -LrZ "\#include \"atom.h\"" *.cpp --exclude="phase_field.cpp" --exclude="position.cpp" --exclude="Box.cpp" | xargs -0 grep -lr "\#include \"position.h\"")

ATOMPROGS := $(shell $(GREP) -lr --exclude="atom.cpp" --exclude="verifyNewFile.cpp" "\#include \"atom.h\"" *.cpp | sed "s/.cpp//g")
ATOMDBGPROGS := $(shell $(GREP) -lr --exclude="atom.cpp" --exclude="verifyNewFile.cpp" "\#include \"atom.h\"" *.cpp | sed 's/.cpp/_dbg /g')
UTILS := $(shell $(GREP) -LrZ "\#include \"atom.h\"" --exclude="atom.cpp" --exclude="position.cpp" --exclude="verifyNewFile.cpp" --exclude="phase_field.cpp" --exclude="Box.cpp" --exclude="LAMMPSDataFileReader.cpp" *.cpp | xargs -0 $(GREP) -Lr "\#include \"position.h\"" | sed 's/.cpp//g')
UTILSDBG := $(shell $(GREP) -LrZ "\#include \"atom.h\"" --exclude="atom.cpp" --exclude="position.cpp" --exclude="verifyNewFile.cpp" --exclude="phase_field.cpp" --exclude="Box.cpp" --exclude="LAMMPSDataFileReader.cpp" *.cpp | xargs -0 $(GREP) -Lr "\#include \"position.h\"" | sed 's/.cpp/_dbg/g')
POSITIONPROGS = $(shell $(GREP) -LrZ "\#include \"atom.h\"" *.cpp --exclude="phase_field.cpp" --exclude="position.cpp" --exclude="verifyNewFile.cpp" | xargs -0 grep -lr "\#include \"position.h\"" | sed 's/.cpp//g')
POSITIONDBGPROGS = $(shell $(GREP) -LrZ "\#include \"atom.h\"" *.cpp --exclude="phase_field.cpp" --exclude="position.cpp" --exclude="verifyNewFile.cpp" | xargs -0 grep -lr "\#include \"position.h\"" | sed 's/.cpp/_dbg/g')
ALLEXECUTABLES = $(ATOMPROGS) $(ATOMDBGPROGS) $(UTILS) $(UTILSDBG) $(POSITIONPROGS) $(POSITIONDBGPROGS)
ALLDBGEXECUTABLES = $(ATOMDBGPROGS) $(UTILSDBG) $(POSITIONDBGPROGS)

ERRORCODES = error_code_defines.h

all: FLAGS = $(CLIB) $(NODBGFLAGS)
all: atom.o position.o atom_builds utils position_builds phase_field Box.o verifyNewFile.o LAMMPSDataFileReader.o

debug: FLAGS = $(CLIB) $(DBGFLAGS)
debug: atom.dbg.o position.dbg.o verifyNewFile.dbg.o atom_builds_dbg utils_dbg position_builds_dbg phase_field_dbg

phase_field: FLAGS = $(CLIB) -fopenmp $(NODBGFLAGS) # sets the FLAGS variable
phase_field: TARGETS = position.o # sets the TARGETS variable
phase_field: $(addprefix $(BIN)/, phase_field) position.o

phase_field_dbg: FLAGS = $(CLIB) -fopenmp $(DBGFLAGS)
phase_field_dbg: TARGETS = position.dbg.o
phase_field_dbg: $(addprefix $(BIN)/, phase_field_dbg) position.dbg.o

# The extra line for each rule allows a different value to be passed to the pattern-matched rule.
atom_builds: TARGETS=atom.o position.o verifyNewFile.o
atom_builds: $(addprefix $(BIN)/, $(ATOMPROGS)) atom.o

utils: TARGETS=verifyNewFile.o
utils: $(addprefix $(BIN)/, $(UTILS))

position_builds: TARGETS=position.o verifyNewFile.o
position_builds: $(addprefix $(BIN)/, $(POSITIONPROGS)) position.o

atom_builds_dbg: TARGETS=atom.dbg.o position.dbg.o verifyNewFile.dbg.o
atom_builds_dbg: $(addprefix $(BIN)/, $(ATOMDBGPROGS)) atom.dbg.o position.dbg.o

utils_dbg: TARGETS=verifyNewFile.dbg.o
utils_dbg: $(addprefix $(BIN)/, $(UTILSDBG))

position_builds_dbg: TARGETS=position.dbg.o verifyNewFile.dbg.o
position_builds_dbg: $(addprefix $(BIN)/, $(POSITIONDBGPROGS)) position.dbg.o

atom.o : atom.h atom.cpp
	$(CC) -O3 $(CXXFLAGS) $(STDLIB) $(CLIB) -c atom.cpp $(ADDITIONAL_LDFLAGS)

atom.dbg.o : atom.h atom.cpp
	$(CC) -ggdb -g -D_GLIBCXX_DEBUG $(CXXFLAGS) $(STDLIB) $(CLIB) -c atom.cpp $(ADDITIONAL_LDFLAGS) -o atom.dbg.o

position.o : position.h position.cpp
	$(CC) -O3 $(CXXFLAGS) $(STDLIB) $(CLIB) -c position.cpp $(ADDITIONAL_LDFLAGS)

position.dbg.o : position.h position.cpp
	$(CC) -ggdb -g -D_GLIBCXX_DEBUG $(CXXFLAGS) $(STDLIB) $(CLIB) -c position.cpp $(ADDITIONAL_LDFLAGS) -o position.dbg.o

verifyNewFile.o : verifyNewFile.h verifyNewFile.cpp
	$(CC) -O3 $(CXXFLAGS) $(STDLIB) $(CLIB) -c verifyNewFile.cpp $(ADDITIONAL_LDFLAGS)

verifyNewFile.dbg.o : verifyNewFile.h verifyNewFile.cpp
	$(CC) -ggdb -g -D_GLIBCXX_DEBUG $(CXXFLAGS) $(STDLIB) $(CLIB) -c verifyNewFile.cpp $(ADDITIONAL_LDFLAGS) -o verifyNewFile.dbg.o

Box.o : Box.h Box.cpp
	$(CC) -O3 $(CXXFLAGS) $(STDLIB) $(CLIB) -c Box.cpp $(ADDITIONAL_LDFLAGS)

Box.dbg.o : Box.h Box.cpp
	$(CC) -ggdb -g -D_GLIBCXX_DEBUG $(CXXFLAGS) $(STDLIB) $(CLIB) -c Box.cpp $(ADDITIONAL_LDFLAGS) -o Box.dbg.o

LAMMPSDataFileReader.o : LAMMPSDataFileReader.h LAMMPSDataFileReader.cpp
	$(CC) -O3 $(CXXFLAGS) $(STDLIB) $(CLIB) -c LAMMPSDataFileReader.cpp $(ADDITIONAL_LDFLAGS)

LAMMPSDataFileReader.dbg.o : LAMMPSDataFileReader.h LAMMPSDataFileReader.cpp
	$(CC) -ggdb -g -D_GLIBCXX_DEBUG $(CXXFLAGS) $(STDLIB) $(CLIB) -c LAMMPSDataFileReader.cpp $(ADDITIONAL_LDFLAGS) -o LAMMPSDataFileReader.dbg.o

$(BIN)/convert2coreshell : $$(addsuffix .cpp, $$(subst _dbg,,$$(@F))) atom.o position.o verifyNewFile.o Box.o $$(ERRORCODES)
	$(CC) $(strip $(FLAGS)) $(addprefix $(BIN)/, $(@F)) $(CXXFLAGS) $(addsuffix .cpp, $(subst _dbg,,$(@F))) atom.o position.o verifyNewFile.o Box.o $(ADDITIONAL_LDFLAGS)

$(BIN)/convert2coreshell_dbg : $$(addsuffix .cpp, $$(subst _dbg,,$$(@F))) atom.dbg.o position.dbg.o verifyNewFile.dbg.o Box.dbg.o $$(ERRORCODES)
	$(CC) $(strip $(FLAGS)) $(addprefix $(BIN)/, $(@F)) $(CXXFLAGS) $(addsuffix .cpp, $(subst _dbg,,$(@F))) atom.dbg.o position.dbg.o verifyNewFile.dbg.o Box.dbg.o $(ADDITIONAL_LDFLAGS)

# Explanation of this rule:
# the % sign is simply what caused this rule to be used.
# $(@F) gives the filename (which in this case is equal to the %, though I couldn't use that)
# subst _dbg,, : subst is the command to swap strings.  the first argument (_dbg)
# is the string to be replaced, the second argument (between the two commas) is
# the replacement string, and the third argument ($(@F) in this case) is the string
# to perform the replacement on.
$(BIN)/% : $$(addsuffix .cpp, $$(subst _dbg,,$$(@F))) $$(TARGETS) $$(ERRORCODES)
	$(CC) $(strip $(FLAGS)) $(addprefix $(BIN)/, $(@F)) $(CXXFLAGS) $(addsuffix .cpp, $(subst _dbg,,$(@F))) $(TARGETS) $(ADDITIONAL_LDFLAGS)

# ls -I '*.*' lists all files that do NOT match the pattern '*.*' - all of the
# executables created from this makefile.  BEWARE CHANGING FORMATS FOR OTHER FILE TYPES
clean :
	rm $(addprefix $(BIN)/, $(ALLEXECUTABLES)); rm *.o

# This allows variables to be printed: e.g. make print-FLAGS will show the value of FLAGS.
print-%  : ; @echo $* = $($*)

.PHONY : all debug phase_field phase_field_dbg atom_builds utils position_builds \
atom_builds_dbg utils_dbg position_dbg_builds clean list
