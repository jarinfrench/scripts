rotate_and_remove : rotate_and_remove_atom.cpp atom.o
	g++ -O3 -o rotate_and_remove rotate_and_remove_atom.cpp atom.o

atom.o : atom.h atom.cpp
	g++ -c atom.cpp
