.SECONDEXPANSION:

UNAME=$(shell uname -s)

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
	DBGFLAGS = -ggdb -g $(STDLIB) -o
	GREP = grep
endif

OBJ_DIR := ./obj
# This command finds all the files in the current directory that have duplicate basenames (without extension)
# In essence, we are looking for files with BOTH a .cpp and .h file
SRC_FILES := $(wildcard ./*.cpp)
OBJ_FILES := $(patsubst ./%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

all: $(OBJ_FILES)

$(OBJ_DIR)/%.o: ./%.cpp
	$(CC) -O3 $(CXXFLAGS) $(STDLIB) $(CLIB) -c -o $@ $<

clean:
	rm *.o

print-%  : ; @echo $* = $($*)

.PHONY: all clean
