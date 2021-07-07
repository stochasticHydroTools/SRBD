program main
   use DoiBoxModule
   use MinHeapModule
   use HydroGridModule
   use BoxLibRNGs
   implicit none

   type (DoiBox), target :: box
   type (HydroGrid), target :: grid

   ! Inputs:
   integer :: nEquilibrationStep=0, nSteps = 100, nSampleStep = 10, nOutputStep=10, nStatsStep=10, seed=0
   namelist / DoiDriverOptions / basename, nSteps, nSampleStep, nOutputStep, nEquilibrationStep, nStatsStep, seed
   
   ! Local variables
   integer :: nameListUnit
   integer :: iStep, particle, species
   real (wp) :: timestep

   character(LEN=1024) :: basename="" ! Base of filename used for reading/writing files
   character(LEN=1024) :: filename="" ! Doi output is saved to this file
   
   ! Executable statements
   !----------------------------------

   nameListUnit = 15
   
   write(*,*) "Enter the basename of the namelist to read (empty for default)"
   read(*,"(A)") basename
   !basename = "" ! Donev: read from default file
   
   if (basename == "") then ! Default values
      write(*,*) "Reading parameters from namelist"
      open (nameListUnit, file=  "input/DoiDriverOptions.nml", status="old", action="read")
   else
      filename=trim(basename) // ".nml"
      write(*,*) "Reading parameters from ", TRIM(filename)
      open (nameListUnit, file=  TRIM(filename), status="old", action="read")
   end if
   read (nameListUnit, nml=DoiDriverOptions)

   if (basename == "") then ! Default values
      basename = "DoiHydro"
      filename="DoiOutput.dat"
   else
      ! The prefix for Doi output files:
      filename=trim(basename) // ".output"   
   end if

   call SeedRNG (seed)

   ! Continue reading the input file
   ! The Fortran standard requires that the namelist be the next non-blank record in the file
   call readDoiParameters (nameListUnit)

   call initializeDoiParameters() ! Called once per processor for the whole Doi domain

   call createDoiBox (box, boundaryType=wallBCs%wallType, lbSample=(/1,1,1/), ubSample=nSampleCells(1:))
   write(*,*) 'nSampleCells=', nSampleCells(1:)

   call initializeDoiBox(box)

   timestep=inputTimestep
   
   if(nSampleStep > 0) then ! Create 
      ! Note: We read the hydro namelist from the same file
      
      ! We pretend like there is a ficticious last species here since values are never calculated for it!
      call createHydroAnalysis (grid, nCells = nSampleCells, nSpecies = nSpecies+1, &
         isSingleFluid = .true., nVelocityDimensions = 0, &
         systemLength = real(nSampleCells*sampleCellLength, dp), &
         timestep = real(nSampleStep*timestep, dp), fileUnit=nameListUnit, fileNamePrefix=TRIM(basename))

   end if

   ! We have now read all the namelists   
   close(nameListUnit)

   write(*,*) 'Starting equilibration loop'
   
   ! Skip a number of steps in the beginning
   do iStep = -nEquilibrationStep, -1 ! Count these as negative time steps
      call updateDoiBox(box,timestep)  ! moves, sorts, and reacts particles.
   end do
   
   write(*,*) 'Starting active loop'

   do iStep = 0, nSteps ! Go one more step so as to write files at the end
               
      if (nOutputStep > 0) then
      if (Mod(iStep, nOutputStep) == 0) then ! Output some values to screen
      ! Donev TODO: Open a file with a proper file name and write to it instead of fort.21 here (see DSMC code for example)
      if(.true.) then ! Output minimal stuff to minimize size of files
         write(9,*) box%globalTime, box%nParticles(1)/domainVolume ! Only first species is enough for stuff like A+B<->C due to conservation laws
      else   
         write(9,*) box%globalTime, box%nParticles(1:)/domainVolume ! Average number density       
         write(21,*) box%globalTime, box%nParticles ! Total number of particles
         write(22,*) box%globalTime, box%reactionCount(1:totalReactions) ! Donev: also write total number of reactions up to this time in a file 
      end if       
      end if
      end if
      
      if (nSampleStep > 0) then
      if (Mod(iStep, nSampleStep) == 0) then
         call getState (box)
         ! Donev: The temperature can be used to save some other scalar field such as an order parameter / collective coordinate
         ! For the A+B->0 a good (conserved) order parameter is the difference A-B         
         ! One can also plot the fraction of the first component but note that this fails if some cell is empty!
         call updateHydroAnalysisPrimitive (grid, density=box%numberDensity(:, :, :, 0), temperature=&
                 !box%numberDensity(:, :, :, 1)/max(1.0_wp/DoiCellVolume, box%numberDensity(:, :, :, 0)), & ! Fraction of first component in non-empty cells
                 !box%numberDensity(:, :, :, 1)-box%numberDensity(:, :, :, 2), & ! For A+B->0 write A-B
                 box%numberDensity(:, :, :, 1)+box%numberDensity(:, :, :, 2)+2*box%numberDensity(:, :, :, 3), & ! For A+B->C write A+B+2*C
                 concentration=box%numberDensity(:, :, :, 1:nSpecies))
      end if
      end if
      
      if (nStatsStep > 0) then
      if (Mod(iStep, nStatsStep) == 0) then ! Write results to files
         write(*,*) iStep, "t=", box%globalTime, " n_particles=", box%nParticles(0), &
            " fractions=", real(box%nParticles(1:nSpecies))/real(box%nParticles(0))
         write(*,*) "Number of reactions is:", box%reactionCount(0), "=", box%reactionCount(1:totalReactions)
         ! This includes statistics gathered from the last time we were here up until now
         if(nSampleStep > 0) call writeToFiles(grid, iStep) ! This will write files from HydroGrid analysis
                  
      end if
      end if

      ! Move to the next point in time:     
      if(iStep<nSteps) call updateDoiBox(box,timestep)  ! moves, sorts, and reacts particles.
       
   end do
   write(*,*) 'Completed time loop'
   write(*,*) 'Total reaction count is ', box%reactionCount(0), "=", box%reactionCount(1:totalReactions)  
   if(IRDME_test) then
      write (*,*) 'Number of fake reactions in 3 categories per unit time is ', box%IRDME_test_number/box%globalTime
   end if 
   
   if(nSampleStep > 0) then
      if(nStatsStep <= 0) call writeToFiles(grid) ! Only write statistics once at the end of the run
      call destroyHydroAnalysis(grid)
   end if   

end program main
