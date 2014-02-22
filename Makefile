

CXX = g++
LAPACK = /Users/mmccull/lib
#OPTS = -O3 -ftree-vectorize -lm -lgomp
OPTS = -O3 -ftree-vectorize -fopenmp -lm
#-ftree-parallelize-loops

dx: psf2dx.cpp psflib.h psflib.cpp pdblib.h pdblib.c
	$(CXX) -c psf2dx.cpp pdblib.c psflib.cpp $(OPTS) 
	$(CXX) psf2dx.o pdblib.o psflib.o $(OPTS) -o psf2dx

clean:
	rm -f *.o
