SRC_DOI=../
# You must also obtain the HydroGrid library from
# https://github.com/stochasticHydroTools/HydroGrid
# and point to the directory here:
SRC_HYDROGRID=../../HydroGrid

SRC_HYDROLIB=$(SRC_HYDROGRID)/src

#======================================
# Choose compilers/optimization flags/directories
#======================================

# C compiler:
CC=gcc -O3 # Optimized
#CC=gcc -g # Debugging

# Fortran compiler:

# GNU:
# We require newer versions of gfortran (>4.6) to allow for C interoperability
FC_O=gfortran -O3 -Wall -ffree-line-length-none  # Optimized
FC_g=gfortran -g -fbounds-check -Wall -ffree-line-length-none # Debugging

# Choose level of optimization:
FC=$(FC_O) # For OPTIMIZED code, use for production code
#FC=$(FC_g) # For DEBUGGING code, use this during development of code

# Some paths for included sources:
HYDROLIB_PATH=$(SRC_HYDROLIB)

#======================================
# Make rules/targets/etc.
#======================================

# Look for sources in the repository:
vpath %.c $(SRC_DOI)
vpath %.h $(SRC_DOI)
vpath %.f90 $(SRC_DOI)
vpath %.c $(SRC_HYDROLIB)
vpath %.h $(SRC_HYDROLIB)
vpath %.f90 $(SRC_HYDROLIB)

# Libraries tolink
# Kishore: Update the directory path /usr/lib64 to whereever FFTW library is
LIBRARIES=-L/usr/lib64/ -lfftw3 -lfftw3f -lm

# List of additional source codes with modules/libraries (not main programs)
CSOURCES=
FSOURCES=

# From HydroGrid library:
CSOURCES += RNGs.c visit_writer.c
FSOURCES += Precision.f90 RNGEngine.f90 RNG.f90 NURNGs.f90 Random.f90 FFTW.f90 VisitWriter.f90

OBJECTS=$(CSOURCES:.c=.o) $(FSOURCES:.f90=.o) 

all:	COM_TwoParticles.exe COM_OneParticle.exe COM_TwoParticlesSpring.exe

clean:
	rm -rf *.o *.mod COM_OneParticle.exe DimerTest.exe

# HydroGrid make rules:
#--------------------------
Precision.o: Precision.f90
	$(FC) -c $<

Random.o: Random.f90 Precision.o
	$(FC) -c $<

RNGEngine.o: RNGEngine.f90
	$(FC) -c $<

NURNGs.o: NURNGs.f90 RNGEngine.o
	$(FC) -c $<

RNG.o: RNG.f90 Random.o RNGEngine.o NURNGs.o
	$(FC) -c $<

RNGs.o: RNGs.c
	$(CC) -c -DBL_FORT_USE_UNDERSCORE $<
                        
visit_writer.o: visit_writer.c
	$(CC) -c $<

FFTW.o: FFTW.f90 Precision.o
	$(FC) -c $<

VisitWriter.o: VisitWriter.f90
	$(FC) -c $<

# Executables:
#--------------------------

DimerTest.exe: DimerTest.f90 $(OBJECTS)
	$(FC) -o $@ $< $(OBJECTS) $(LIBRARIES) 

COM_OneParticle.exe: COM_OneParticle.f90 $(OBJECTS)
	$(FC) -o $@ $< $(OBJECTS) $(LIBRARIES) 

COM_TwoParticles.exe: COM_TwoParticles.f90 $(OBJECTS)
	$(FC) -o $@ $< $(OBJECTS) $(LIBRARIES)

COM_TwoParticlesSpring.exe: COM_TwoParticlesSpring.f90 $(OBJECTS)
	$(FC) -o $@ $< $(OBJECTS) $(LIBRARIES) 
 
