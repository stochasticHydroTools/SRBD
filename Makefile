SRC_DOI=./
# You must also obtain the HydroGrid library from
# https://github.com/stochasticHydroTools/HydroGrid
# and point to the directory here:
SRC_HYDROGRID=../HydroGrid

SRC_HYDROLIB=$(SRC_HYDROGRID)/src

#======================================
# Choose compilers

# C compiler:
CC=gcc -O3 # Optimized
#CC=gcc -g # Debugging

# Fortran compiler:

# GNU:
# We require newer versions of gfortran (>4.6) to allow for C interoperability
FC_O=gfortran -O3 -Wall -ffree-line-length-none  # Optimized
FC_g=gfortran -Og -g -fbounds-check -Wall -ffree-line-length-none # Debugging

# Choose level of optimization:
#FC=$(FC_O) # For OPTIMIZED code, use for production code
FC=$(FC_g) # For DEBUGGING code, use this during development of code

# We no longer support the old OpenGL-enabled HydroGrid code:
include $(SRC_DOI)/MakefileCommon
