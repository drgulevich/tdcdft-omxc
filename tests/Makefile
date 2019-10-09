.PHONY : clean

LIBS_LAPACK = -llapacke -llapack -lblas -lgfortran
LIBS = -lm -fopenmp

PLATFORM := $(shell uname)
ifeq ($(PLATFORM),Linux)
$(info Linux OS detected)
SONAME_FLAG = -soname
endif
ifeq ($(PLATFORM),Darwin)
$(info Mac OS detected)
SONAME_FLAG = -install_name
endif


DEBUG = -DARMA_EXTRA_DEBUG
## Uncomment the above line to enable low-level debugging.
## Lots of debugging information will be printed when a compiled program is run.
## Please enable this option when reporting bugs.

FINAL = -DARMA_NO_DEBUG
## Uncomment the above line to disable Armadillo's checks.
## Not recommended unless your code has been first thoroughly tested!

CXXFLAGS  = $(FINAL) -march=native -Ofast -std=c++11 -Wall
#CXXFLAGS  = $(DEBUG) -g -Wall
#CXXFLAGS  = $(FINAL) -g -Wall
#CXXFLAGS  = $(FINAL) -pg # option for gprof profiler


### Tests

OBJ = ../xc.o ../dft.o ../tddft.o main.o test1.o test2.o test3.o test4.o
KERNELS = kernels

tests: $(OBJ) 
	g++ $(CXXFLAGS) $(OBJ) -o tests -larmadillo

main.o: main.cpp
	g++ $(CXXFLAGS) -c $<

test%.o: test%.cpp ../qwmodels.hpp
	g++ $(CXXFLAGS) -c $<

#test%: test%.o $(OBJ)
#	g++ $(CXXFLAGS) $(OBJ) $< -o $@ -larmadillo


### 

dft.o: ../dft.cpp ../dft.hpp ../systems.hpp
	g++ $(CXXFLAGS) -c ../dft.cpp -std=c++17

tddft.o: ../tddft.cpp ../tddft.hpp ../systems.hpp $(KERNELS)
	g++ $(CXXFLAGS) -c ../tddft.cpp -std=c++17

xc.o: ../xc.cpp ../xc.hpp
	g++ $(CXXFLAGS) -c ../xc.cpp -std=c++17

clean :
	-rm -vf $(OBJ) tests
