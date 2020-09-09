# autotransit
Code to model photometry of transiting exoplanets  

There is a bit of setup.  We need to compile a few python libraries.  

For gfortran users:

```
make
```

This will create a binary in the bin directory called `claretquad_tess` and
two shared (.so) libraries `detrend5` and `tfit5`.  

Then AutoTransit.ipynb should work.

If you are using ifort (my preferred compiler). Then edit the Makefile to have lines 2-5 as,
```
F77 = $(FPATH)ifort
F90 = $(FPATH)ifort
#F2PYFLAG =
F2PYFLAG = --fcompiler=intelem
```

Enjoy!
