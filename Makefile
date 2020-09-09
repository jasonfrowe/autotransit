#Name of Fortran compiler
FPATH = 
F77 = $(FPATH)gfortran
F90 = $(FPATH)gfortran
F2PYFLAG = 
#F2PYFLAG = --fcompiler=intelem
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

all: claretquad_tess clean

limbpriorsincl = precision.o calcldprior.o locate.o trilinear.o lininterp.o
claretquad_tess: $(UTILS)claretquad_tess.f90 $(limbpriorsincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ $(UTILS)claretquad_tess.f90 $(limbpriorsincl)
	cd $(UTILS)
	f2py3 -c $(UTILS)detrend5.pyf $(UTILS)polyfilter_ramp.f $(UTILS)rqsort.f $(UTILS)gaussj.f $(UTILS)stdev.f $(F2PYFLAG)
	f2py3 -c $(UTILS)tfit5.pyf $(UTILS)transitmodel.f $(UTILS)keplerian.f $(UTILS)ttcor.f $(UTILS)occultquad.f $(UTILS)mandelagol.f $(UTILS)rqsort.f $(UTILS)transitdur.f $(F2PYFLAG)

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
