SHELL = /bin/bash
CXX=g++
CXXFLAGS= -O3 -fopenmp
LDLIBS= -lgsl -lgslcblas -lm -fopenmp
LDFLAGS= -O3 
OBJS= TriCubicInterpolator.o FN_Knot.o
DEPS=FN_Knot.h FN_Constants.h TriCubicInterpolator.h

all: preprocessing executable clean
release: preprocessing executable clean
debug: CXXFLAGS = -Wall -g
debug: LDFLAGS = -Wall -g
debug: preprocessing executable clean

preprocessing:
	INSERT_INITIALISATION_TYPE=FROM_CURVE_FILE ; \
	INSERT_SURFACE_FILENAME=test_three1 ; \
	INSERT_UV_FILENAME=uv_plot0.vtk ; \
	INSERT_RUNTIME=21 ; \
	INSERT_UVPRINTTIME=5 ; \
	INSERT_VELOCITYPRINTTIME=10 ; \
	INSERT_FREQUENTPRINTTIME=10 ; \
	INSERT_SKIPTIME=0 ; \
	INSERT_GRIDSPACING=0.6 ; \
	INSERT_NX=200 ; \
	INSERT_NY=200 ; \
	INSERT_NZ=200 ; \
	INSERT_TIMESTEP=0.02 ; \
	names="INSERT_INITIALISATION_TYPE INSERT_SURFACE_FILENAME INSERT_UV_FILENAME INSERT_RUNTIME INSERT_UVPRINTTIME INSERT_VELOCITYPRINTTIME INSERT_FREQUENTPRINTTIME INSERT_SKIPTIME INSERT_GRIDSPACING INSERT_NX INSERT_NY INSERT_NZ INSERT_TIMESTEP" ; \
	cp FN_Constants.h FN_Constants_backup.h ; \
	for name in $${names} ; do  value=$${!name} ; sed "s/$${name}/$${value}/" FN_Constants.h > temp ; mv -f temp FN_Constants.h ; done 
executable:$(OBJS)
	$(CXX) -o FN_Knot $(OBJS) $(LDLIBS) $(LDFLAGS)

%.o: %.c $(DEPS)
	$(CXX) $(CXXFLAGS)-c -o $@ $< 

clean:
	rm -f *.o
	cp -f FN_Constants.h FN_Constants_written.h
	mv -f FN_Constants_backup.h FN_Constants.h

.PHONY: clean
	
# a handy list of commands I was using for different architectures
#icpc -O3 -qopenmp FN_Knot.cpp TriCubicInterpolator.cpp -lgsl -lgslcblas -lm	
