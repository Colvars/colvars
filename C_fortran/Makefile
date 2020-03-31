# -*- makefile -*-

# This Makefile is useful for compilation tests only.  To use the library,
# please use instead the build recipes of the packages that include it.

COLVARS_LIB = ../src/libcolvars.a
COLVARS_SRC_DIR = ../src
IF_OBJ = colvarproxy_C.o

CXXFLAGS := -std=c++11 -pedantic -g -O2 -fPIC

.PHONY: default $(COLVARS_LIB)

default: fortran_test C_test libcolvars_C.so

fortran_test: fortran_test.o $(COLVARS_LIB) $(IF_OBJ)
	g++ -o $@ $<  $(IF_OBJ) $(COLVARS_LIB) -lgfortran

fortran_test.o: fortran_test.f90
	gfortran -c $<

C_test: C_test.o $(COLVARS_LIB) $(IF_OBJ)
	g++ -o $@ $<  $(IF_OBJ) $(COLVARS_LIB)

C_test.o: C_test.c
	gcc -c $< -pedantic -g  -I$(COLVARS_SRC_DIR)

$(COLVARS_LIB):
	make -f $(COLVARS_SRC_DIR)/Makefile -C $(COLVARS_SRC_DIR) libcolvars.a -j20

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<  -I$(COLVARS_SRC_DIR)

libcolvars_C.so: $(COLVARS_LIB) $(IF_OBJ)
	g++ -shared -o $@ $(IF_OBJ) $(COLVARS_LIB)
