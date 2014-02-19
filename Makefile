

CXX = g++
LAPACK = /Users/mmccull/lib
OPTS = -O3 -ftree-vectorize -fopenmp 
#-ftree-parallelize-loops

dx: psf2dx.cpp psflib.h psflib.cpp stringlib.h stringlib.cpp pdblib.h pdblib.c
	$(CXX) -c psf2dx.cpp stringlib.cpp pdblib.c psflib.cpp $(OPTS) 
	$(CXX) psf2dx.o stringlib.o pdblib.o psflib.o $(OPTS) -L$(LAPACK) -lclapack -lcblas -lf2c -lm -o psf2dx

