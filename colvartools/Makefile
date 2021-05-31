SHELL=/bin/sh
#########################
# adjust as needed
# settings for linux with GCC
COLVARS_SRC=../src
CXX=g++
CXXFLAGS=-Wall -O2 -I$(COLVARS_SRC)
EXT=
# settings for mingw cross-compiler to windows (32-bit)
#CXX=i688-w64-mingw64-g++
#CXXFLAGS=-Wall -O2
#EXT=.exe
#########################

all: poisson_integrator$(EXT) abf_integrate$(EXT)

clean:
	-rm *~ *.o abf_integrate$(EXT) poisson_integrator$(EXT) *.exe

abf_integrate$(EXT): abf_integrate.o abf_data.o
	$(CXX) -o $@ $(CXXFLAGS) $^

poisson_integrator$(EXT): poisson_integrator.o libcolvars
	$(CXX) -o $@ $(CXXFLAGS) poisson_integrator.o $(COLVARS_SRC)/libcolvars.a

poisson_integrator_conv: poisson_integrator_conv.o libcolvars
	$(CXX) -o $@ $(CXXFLAGS) poisson_integrator_conv.o $(COLVARS_SRC)/libcolvars.a

%.o: %.cpp
	$(CXX) -o $@ -c $(CXXFLAGS) $<

libcolvars:
	$(MAKE) -C $(COLVARS_SRC) libcolvars.a

# dependencies
abf_integrate.o: abf_integrate.cpp abf_data.h
abf_data.o: abf_data.cpp abf_data.h

