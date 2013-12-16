SHELL=/bin/sh
#########################
# adjust as needed
# settings for linux with GCC
CXX=g++
CXXFLAGS=-Wall -O2
EXT=
# settings for mingw cross-compiler to windows (32-bit)
#CXX=i688-w64-mingw64-g++
#CXXFLAGS=-Wall -O2
#EXT=.exe
#########################

all: abf_integrate$(EXT)

clean:
	-rm *~ *.o abf_integrate$(EXT) *.exe

abf_integrate$(EXT): abf_integrate.o abf_data.o
	$(CXX) -o $@ $(CXXFLAGS) $^

%.o: %.cpp
	$(CXX) -o $@ -c $(CXXFLAGS) $<

# dependencies
abf_integrate.o: abf_integrate.cpp abf_data.h
abf_data.o: abf_data.cpp abf_data.h

