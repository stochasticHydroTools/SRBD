# SRBD: Stochastic Reactive Brownian Dynamics

This repository provides codes for the SRBD particle method for reaction-diffusion problems described in the paper:

**Efficient Reactive Brownian Dynamics**, A. Donev, C.-Y. Yang and C. Kim.  arXiv:1710.02232 (2017) [arXiv](https://arxiv.org/abs/1710.02232).

This documentation is not very complete and the codes are research codes written for a summer undergraduate research project. Nevertheless, the code is rather efficient and sufficiently flexible that we hope you find it of use -- contact [Aleksandar Donev](http://cims.nyu.edu/~donev) with any questions. Note that this code only works for **periodic** systems.

The code uses in a non-essential way the [_HydroGrid_ library](https://github.com/stochasticHydroTools/HydroGrid) to generate random numbers and analyze the results, so you must get that library first. It is not difficult to remove this dependency and customize the code to your own needs.

## 1. Codes

The code is written in Fortran 95 with some extensions from Fortran 2003, and can be compiled with gfortran or any modern Fortran compiler. We provide a `Makefile` so that simply executing `make` should work once you edit the top of the `Makefile` to select your compiler and the location of `HydroGrid`. Successful compilation will generate the executable `DOI.exe`.

The code was written on top of a DSMC code and it is not as modular as it could be. Most importantly, there are some constants that are hard-wired as Fortran `parameter` values, and need to specified at compile time and not at runtime. This is not hard to improve and if you do it please let us know. In particular, you need to edit the following lines at the top of the main source code `DoiBoxModule.f90' each time you want to study a new problem:

---

```
integer, parameter :: nSpecies = 3, nReactions = 2 ! This is for 0->A+B, C->A+B
integer, parameter :: nDimensions = 2 ! Run in 1, 2 or 3 dimensions
integer, dimension (0:nMaxDimensions), parameter :: neighborhoodSize = (/9, 1, 1, 0/) ! 2D
```

---

Note that we provide the 3 required values of `neighborhoodSize`, so just remove the comment mark `!` for the line corresponding to `nDimensions` and comment out the remaining lines.

The code uses the `HydroGrid` library -- the routines `UniformRNG` and `UniformRNGVec` to generate standard random numbers uniformly distributed in [0,1), the routine `NormalRNG` to provide normally-distributed numbers, and the routine `UniformInteger` to generate a random integer from a given interval. The routine `PoissonRNG` is only used for initialization. By replacing these routines with your own you can switch to another Random Number Generator (RNG).

Note that this particle code is serial and it is nontrivial to parallelize it due to its event-driven nature. There are comments in the code marked `PARALLEL` for those parts that are easy to parallelize using OpenMP; if you do that please let us know.

WARNING: The names used in the code are not the same as those used in the paper, for example, the types of reactions. Notably, SRBD is called IRDME in the code (for Isotropic Reaction-Diffusion Master Equation, after the I-DSMC algorithm that the method is based on), and what the paper calles S-BD-RME is called (somewhat misleadingly) RDME. The key difference with the RDME is that in our code particles diffuse off latice. Nevertheless, for completeness, we have also implemented particles hopping on a lattice, however, this code is not efficient for this purpose. Instead, for an efficient RDME code please consult the fluctuating hydrodynamics code [FHD_ReactDiff](https://github.com/BoxLib-Codes/FHD_ReactDiff) which is also parallel. 


## 2. Input file

The input parameters are set via an input file that takes the form of a Fortran _namelist_ file. An example with some comments is provided in `DoiDriverOptions.nml`. Most of the input fields should be obvious but not all. The namelist `DoiDriverOptions` is read by `main.f90` and sets some basic parameters like number of time steps and how often to save statistics. The main namelist is `DoiBoxOptions` which sets all of the parameters for the SRBD algorithm.

One thing that is likely to be confusing at first is the fact that the code distinguishes between two kind of grids and thus grid cells. The first is a _hydro_ or _sampling_ grid used in initialization and for collecting statistics. It is how particle data is converted into data on a grid for potential coupling to continuum fluid solvers and for output of plot files via `HydroGrid`. This grid has nothing to do per se with the grid used to process Doi collisions, which can be finer. These two grids must be commensurate in the sense that:

`nBlockingSample/nBlockingCollisions = DoiCellLength / sampleCellLength`

That is, `nBlockingSample` sampling/hydro cells equal `nBlockingCollisions` collision cells

Note that the size of the domain and the number of cells can be set independently along each Cartesian dimension.

## 3. Examples

There are several example input files in the directory `Examples`, which are related to test problems studied in the paper describing the SRBD algorithm. For example, to run the reversible association model `A+B<->C`, use

`echo ABC-reversible-ACF-2D | ../DOI.exe`

Remember that the code has to be recompiled each time you switch to a new reaction model or change the dimensionality!
