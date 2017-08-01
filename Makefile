SHELL = /bin/bash
CXX=g++
CXXFLAGS= -O3 -fopenmp
LDLIBS= -lm -fopenmp
LDFLAGS= -O3 
OBJS= SolidAngleGaussBonnet.o
DEPS=SolidAngleGaussBonnet.h Constants.h 

all: executable clean
release: executable clean
debug: CXXFLAGS = -Wall -g
debug: LDFLAGS = -Wall -g
debug: executable clean

executable:$(OBJS)
	$(CXX) -o SolidAngleGaussBonnet $(OBJS) $(LDLIBS) $(LDFLAGS)

%.o: %.c $(DEPS)
	$(CXX) $(CXXFLAGS)-c -o $@ $< 

clean:
	rm -f *.o

.PHONY: clean
	
# a handy list of commands I was using for different architectures
#icpc -O3 -qopenmp FN_Knot.cpp TriCubicInterpolator.cpp -lgsl -lgslcblas -lm	
