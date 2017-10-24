SHELL = /bin/sh

# Path to LAMMPS extraction directory
LAMMPS_ROOT = ../../..
LAMMPS_SRC = $(LAMMPS_ROOT)/src

FC = /opt/local/bin/mpif90  # replace with your Fortran compiler
CXX = /opt/local/bin/mpicxx # replace with your C++ compiler

# Flags for Fortran compiler, C++ compiler, and C preprocessor, respectively
FFLAGS = -O3 -fPIC
CXXFLAGS = -O3 -fPIC
CPPFLAGS =

all : liblammps_fortran.so simple.o driver.x 

driver.x :
	$(FC) $(FFLAGS) simple.o LAMMPS.o LAMMPS-wrapper.o \
	    -L${LAMMPS_SRC} -llammps

simple.o : simple.f90 liblammps_fortran.so 
	$(FC) $(FFLAGS) -I${LAMMPS_ROOT}/examples/COUPLE/fortran2 -c simple.f90

liblammps_fortran.so : LAMMPS.o LAMMPS-wrapper.o
	$(FC) $(FFLAGS) -shared -o liblammps_fortran.so LAMMPS.o LAMMPS-wrapper.o -L${LAMMPS_SRC} -llammps

LAMMPS.o lammps.mod : LAMMPS.F90
	$(FC) $(CPPFLAGS) $(FFLAGS) -I${LAMMPS_SRC} -c LAMMPS.F90

LAMMPS-wrapper.o : LAMMPS-wrapper.cpp LAMMPS-wrapper.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I${LAMMPS_SRC} -c LAMMPS-wrapper.cpp

clean :
	$(RM) *.o a.out *.x *.mod liblammps_fortran.so
