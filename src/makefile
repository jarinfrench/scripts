.SECONDEXPANSION:

UNAME = $(shell uname -s)

ifeq ($(UNAME), Darwin)
	CC = clang++
	CLIB = -stdlib=libc++
	STDLIB = -std=gnu++11
	ADDITIONAL_LDFLAGS =
	NODBGFLAGS = -O3 $(STDLIB) -fopenmp -o
	DBGFLAGS = -ggdb -g $(STDLIB) -fopenmp -o
	GREP = ggrep
endif
ifeq ($(UNAME), Linux)
	CC = g++
	CLIB =
	STDLIB = -std=c++11
	ADDITIONAL_LDFLAGS =
	NODBGFLAGS = -O3 $(STDLIB) -fopenmp -o
	DBGFLAGS = -ggdb -g $(STDLIB) -fopenmp -o
	GREP = grep
endif

BIN =../bin

# This may need some fine tuning between Mac and Linux
ATOMPROGS := $(shell $(GREP) -lr --exclude="atom.cpp" "\#include \"atom.h\"" *.cpp | sed "s/.cpp//g")
ATOMDBGPROGS := $(shell $(GREP) -lr --exclude="atom.cpp" "\#include \"atom.h\"" *.cpp | sed 's/.cpp/_dbg /g')
UTILS := $(shell $(GREP) -LrZ "\#include \"atom.h\"" --exclude="atom.cpp" --exclude="position.cpp" *.cpp | xargs -0 $(GREP) -Lr "\#include \"position.h\"" | sed 's/.cpp//g')
UTILSDBG := $(shell $(GREP) -LrZ "\#include \"atom.h\"" --exclude="atom.cpp" --exclude="position.cpp" *.cpp | xargs -0 $(GREP) -Lr "\#include \"position.h\"" | sed 's/.cpp/_dbg/g')
POSITIONPROGS = phase_field
POSITIONDBGPROGS = phase_field_dbg
ALLEXECUTABLES = $(ATOMPROGS) $(ATOMDBGPROGS) $(UTILS) $(UTILSDBG) $(POSITIONPROGS) $(POSITIONDBGPROGS)

ERRORCODES = error_code_defines.h

all: FLAGS = $(CLIB) $(NODBGFLAGS)
all: atom.o position.o atom_builds utils position_builds

debug: FLAGS = $(CLIB) $(DBGFLAGS)
debug: atom_builds_dbg utils_dbg position_dbg_builds

# The extra line for each rule allows a different value to be passed to the pattern-matched rule.
atom_builds: TARGETS=atom.o position.o
atom_builds: $(addprefix $(BIN)/, $(ATOMPROGS)) atom.o

utils: TARGETS=
utils: $(addprefix $(BIN)/, $(UTILS))

position_builds: TARGETS=position.o
position_builds: $(addprefix $(BIN)/, $(POSITIONPROGS)) position.o

atom_builds_dbg: TARGETS=atom.o position.o
atom_builds_dbg: $(addprefix $(BIN)/, $(ATOMDBGPROGS)) atom.o

utils_dbg: TARGETS=
utils_dbg: $(addprefix $(BIN)/, $(UTILSDBG))

position_dbg_builds: TARGETS=position.o
position_dbg_builds: $(addprefix $(BIN)/, $(POSITIONDBGPROGS)) position.o

atom.o : atom.h atom.cpp
	$(CC) -O3 $(STDLIB) $(CLIB) -c atom.cpp $(ADDITIONAL_LDFLAGS)

position.o : position.h position.cpp
	$(CC) -O3 $(STDLIB) $(CLIB) -c position.cpp $(ADDITIONAL_LDFLAGS)

# Explanation of this rule:
# the % sign is simply what caused this rule to be used.
# $(@F) gives the filename (which in this case is equal to the %, though I couldn't use that)
# subst _dbg,, : subst is the command to swap strings.  the first argument (_dbg)
# is the string to be replaced, the second argument (between the two commas) is
# the replacement string, and the third argument ($(@F) in this case) is the string
# to perform the replacement on.
$(BIN)/% : $$(addsuffix .cpp, $$(subst _dbg,,$$(@F))) $$(TARGETS) $$(ERRORCODES)
	$(CC) $(strip $(FLAGS)) $(addprefix $(BIN)/, $(@F)) $(addsuffix .cpp, $(subst _dbg,,$(@F))) $(TARGETS) $(ADDITIONAL_LDFLAGS)

# ls -I '*.*' lists all files that do NOT match the pattern '*.*' - all of the
# executables created from this makefile.  BEWARE CHANGING FORMATS FOR OTHER FILE TYPES
clean :
	rm $(addprefix $(BIN)/, $(ALLEXECUTABLES)); rm *.o

# This allows variables to be printed: e.g. make print-FLAGS will show the value of FLAGS.
print-%  : ; @echo $* = $($*)

.PHONY : all debug clean list
