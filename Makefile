#Name of Fortran compiler
FPATH = 
F77 = $(FPATH)ifort
F90 = $(FPATH)ifort
#F77 = gfortran
#compiling object file flags
OPT1 = -O3 
OPT2 = 
#OPT1 = -O0 -g -CB -warn 
#OPT2 = 
FFLAGS = -c $(OPT1) $(OPT2)
#FFLAGS = -c -O0 -g -CB -warn $(OPT2)
#linking flags
LFLAGS = $(OPT1) $(OPT2) 
#LFLAGS = -O0 -g -CB -warn $(OPT2)
#testing flags
#LFLAGS = -O0 -g -CB
#Compiler flag for MPI
#Directory where executable are placed
BIN = ./bin/
#utils source directory
UTILS = ./utils/

#Listing of programs to create.

all: claretquad_tess 

limbpriorsincl = precision.o calcldprior.o locate.o trilinear.o lininterp.o
claretquad_tess: $(UTILS)claretquad_tess.f90 $(limbpriorsincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ $(UTILS)claretquad_tess.f90 $(limbpriorsincl)

#building object libraries
%.o : $(UTILS)%.f90
	$(F90) $(FFLAGS) -o $@ $<

%.o : $(UTILS)%.f
	$(F77) $(FFLAGS) -o $@ $<

# Removing object files
.PHONY : clean
clean :
	rm *.o
	rm *.mod
