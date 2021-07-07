! Note: Reservoirs are not really working at present, only periodic BCs are working
module DoiBoxModule
   use MinHeapModule
   use BoxLibRNGs ! Random number generation
   use DiffusionCLs
   implicit none
   save
   
   ! Change for different reaction schemes and recompile:
   !integer, parameter :: nSpecies = 3, nReactions = 7 ! BPM model
   integer, parameter :: nSpecies = 3, nReactions = 2 ! A<->B+C test in IRDME Kishore: For our purposes is it enough to change nSpecies = 4?
   !integer, parameter :: nSpecies = 2, nReactions = 2 ! A+B->B, 0->A equilibrium (from Radek Erban)
   !integer, parameter :: nSpecies = 2, nReactions = 2 ! A+B->0, A+A->0, rate test in IRDME
   !integer, parameter :: nSpecies = 3, nReactions = 4 ! 0->A, A->0, A+B->0, A+A->0, test if selection is uniform for each cell/particle/pair
   !integer, parameter :: nSpecies = 3, nReactions = 2 ! 0->A+B, C->A+B
   !integer, parameter :: nSpecies = 1, nReactions = 2 ! A+A->A, 0->A

!--------------------------------------
! Supported reaction types as recorded in reactionTable and reactionType
!--------------------------------------
!
! 1, Birth:  products are created at randomly-sampled positions:
!    (1)  0 --> A (2) 0 --> A + B (3) 0 --> A + A
!    For RDME, the products are placed uniformly inside the cell
!    For IRDME, the products are placed within a reactive distance of each other:
!       one of the particles is uniformly distributed inside the cell, and the other, if there exists one, is within a reactive sphere
! 2, Annihilation:  no product is created:
!    (1) A --> 0 (2) A + A --> 0 (3) A + B --> 0
! 3, Replacement:  one particle gets replaced by one or two particles:
!    (1) A --> B (2) A --> B + C
!    The particle A is replaced by another species, and an extra particle, if required, is randomly placed
!       For RDME, the second particle is placed uniformly inside the cell
!       For IRDME, the second particle is placed within a reactive radius of the original A
! 4, Catalitic Birth: one particle remains unchanged in the reaction and creates a new particle, inverse of catalitic annihilation
!    (1) A -> A + A (2) A -> A + B
!    We simply create a new particle in this reaction near an existing one of species A
!    For RDME, the new particle is created randomly inside the same cell.
!    For IRDME, the new particle is created within a reactive sphere of the A.
! 5, Merge:  one particle disappears and the species of the other particle changes:
!    (1) A + A --> B (2) A + B --> C
!    We randomly choose one of the reagents to disappear and let the other one be replaced by a new particle with
!    the same species as the product
! 6, Catalitic Annihilation (and also Coagulation): one particle disappears when there is another particle nearby:
!    (1) A + B --> A (2) A + A --> A
!    For type (1), we keep the A where it is and make the B disappear
!    For type (2), we randomly choose one of the two reagents to disappear and the other one remains unchanged.
! 7, Transform:  the species of two particles change (randomly choose which one of the two becomes which): 
!                    (1) A + A --> B + C (2) A + B --> C + D (3) A + A --> A + B (4) A + B --> A + A 
!    We randomly let one of the reagent be replaced by one of the product and the other reagent be replaced by the other product.
! 8, Catalyst:    a particle changes species when another one is nearby: A + B --> A + C
!    We keep the catylist A unchanged and turn the B into a C

   integer, parameter :: nMaxDimensions = 3 ! Up to 3 dimensions

   ! Constants that require recompilation
   !---------------------
   ! The way this works is that you comment one of these two lines out (just as above) -- *don't* change the actual lines
   ! If you want something different add a new line with a comment
   
   integer, parameter :: dp=kind(0.0d0), sp=kind(0.0)
   ! The precision is now chosen in MinHeapModule.f90
   !integer, parameter :: wp = dp ! Double precision
   !integer, parameter :: wp = sp ! Single precision
   
   ! Choose dimensionality:
   integer, parameter :: nDimensions = 3 ! 1D, 2D and 3D work but remember to change lines below also
   !integer, dimension (0:nMaxDimensions), parameter :: neighborhoodSize = (/3, 1, 0, 0/) ! 1D
   integer, dimension (0:nMaxDimensions), parameter :: neighborhoodSize = (/9, 1, 1, 0/) ! 2D
   !integer, dimension (0:nMaxDimensions), parameter :: neighborhoodSize = (/27, 1, 1, 1/) ! 3D
   
   ! IMPORTANT: Turn off assertions and tracing/testing for getting optimized code!
   !logical, parameter :: isAssertsOn =.true., printDebug=.true. ! Debug
   !logical, parameter :: isAssertsOn =.true., printDebug=.false. ! Test but not print too much stuff (slow!)
   logical, parameter :: isAssertsOn =.false., printDebug=.false. ! Optimized for production runs
   
   ! For testing the code we do fake reactions and just count to compute the rate in the absence of spatial correlations
   ! Tracing is used to confirm that each pair of overlapping particles is selected with the same probability
   ! IRDME_trace distinguishes two kinds of testing: full history trace and just counting:
   !logical, parameter :: IRDME_test=.true., IRDME_trace=.true. ! Trace the number of times each pair is selected
   !logical, parameter :: IRDME_test=.true., IRDME_trace=.false. ! Count number of reactions only
   logical, parameter :: IRDME_test=.false., IRDME_trace=.false. ! For production runs after testing is complete
   integer :: nMaxTestPairs = 1000000 ! How much space to allocate to keep a history of reacting pairs
   integer :: MaxReactionsPerParticle = 10000 ! How many reactions are allowed for each particle in fake run

   ! Enumerators:
   !---------------------
   integer, parameter :: PERIODIC = 1, FREE = 0, RESERVOIR = -1 ! Types of boundaries. Only PERIODIC works at present
   ! Different types of chemical reactions:
   integer, parameter :: R_NONE=0, R_BIRTH = 1, R_ANNIHILATION=2, R_REPLACEMENT = 3, &
                         R_CATA_BIRTH = 4,  R_MERGE = 5, R_CATA_ANNIHILATION = 6, R_TRANSFORM = 7, R_CATALYST = 8
   integer, parameter :: RDME = 0, IRDME = 1 ! Type of reaction scheme
   
   ! Constants
   real (wp), parameter :: PI = 3.14159265358979_wp

   ! Derived types:
   !----------------
   
   type ParticleType
      real (wp) :: position (nMaxDimensions) ! We keep here three positions even for 2D to simplify indexing
      integer :: species = 0
   end type ParticleType

   type WallType ! Boundary condition for one side of the Doi box
      integer :: wallType = PERIODIC ! What kind of wall this is
      real (wp) :: number_density (1:nSpecies) = -1.0_wp ! For reservoirs, the target number density
   end type WallType

   type DoiBox
      
      integer*8 :: nTimesteps = 0 ! More accurate to count integers than just elapsed time
      real(wp) :: globalTime=0.0_wp ! Keep track of elapsed time
      real(wp) :: elapsedTime=0.0_wp ! In IRDME, used to track elapsed time in each timestep

      ! BCs for this specific box (will be different from the global BCs below if running in parallel)
      integer, dimension (2,nDimensions) :: boundaryType ! Type of boundary for each side of the Doi box      
      real (wp), dimension (2, nDimensions) :: boundaryLocation            
      
      ! IRDME neighrboring cells
      integer, allocatable :: neighboringCells (:,:)

      ! IRDME heap
      type (MinHeap) :: heap

      ! Particle data
      integer :: nMaxParticles ! Computed later
      type (ParticleType), allocatable :: particle (:)
      integer, allocatable, dimension(:)  :: nextParticle, particlesSortedByCell
      !integer, allocatable, dimension(:)  :: previousParticle ! We no longer use a doubly-linked list, we use a singly-linked list

      ! Species data
      integer, dimension(0:nSpecies) :: speciesSum
      integer :: nParticles (0:nSpecies) = 0 ! Number of particles of each species

      ! Reaction cells:
      integer, dimension (nMaxDimensions) :: lb, ub ! Lower and upper bounds for Doi collision cell arrays
      integer, dimension (0:nMaxDimensions) :: extent ! Extent of grid in each dimension (0 for total number of cells)
      integer, allocatable, dimension(:, :, :, :) :: speciesPointers ! Used in RDME for sorting in cells
      integer, pointer, dimension(:,:) :: speciesPointersList ! Alias to above of dimension (0:nSpecies,nCells)
      integer, allocatable, dimension(:, :) :: numberPerCellList, firstParticleList ! Of dimension (0:nSpecies,nCells)
      ! Aliased pointers to numberPerCellList, firstParticleList of dimension (0:nSpecies,nCells(1),..,nCells(nDimensions))
      integer, pointer, dimension(:, :, :, :) :: numberPerCell, firstParticle
      integer :: freeParticle = 1 ! Used to find new particles for insertion

      ! Number densities are calculated for each species on a hydro sampling grid
      ! Use index 0 for the total over all species
      integer, dimension (nMaxDimensions) :: lbSample, ubSample
      integer :: nHydroCells(0:nMaxDimensions)
      ! This next one needs to be double precision since HydroGrid expects that
      real (dp), allocatable :: numberDensity (:, :, :, :) ! Particles per unit volume in each species

      ! Some counters for various statistics
      integer*8, dimension(0:2*nReactions) :: reactionCount=0 ! Counts the number of reactions that happened for each reaction
      
      ! For testing:
      ! QY: The following records the number of particles of each species at the beginning. It is used for debugging only.
      integer, dimension(0:nSpecies):: initialParticleCount      
      ! QY: temporary testing, for recording the participating particles in each reaction
      integer, allocatable :: IRDME_test_number_1(:) ! reactions in each cell
      integer, allocatable :: IRDME_test_number_2(:) ! reactions in each particle
      integer :: IRDME_test_number(4) ! QY: number of each fake reaction in test (1=birth, 2=unary death, 3=A+B->0, 4=A+A->0)
      integer, dimension(:,:), allocatable :: reactionParticle_AA, reactionParticle_AB ! History log of particles in each reaction, AB for A+B->0, 
      !AA for A+A->0
      real (wp), allocatable :: timeOfReactions_AA (:), timeOfReactions_AB (:) ! History of the time of each pairwise reaction recorded in reactionParticle
      ! AA for A+A->0, AB for A+B->0
      integer, dimension(:), allocatable :: count

   end type DoiBox
   
   ! Global Doi box parameters:
   !----------------

   ! Species:
   real (wp) :: speciesDiameter (nSpecies) = 1.0_wp
   real (wp) :: speciesDiffusivity (nSpecies) = 1.0_wp ! In units of length^2 / unit_time

   ! Global reference or default equilibrium values:
   real (wp) :: numberDensity (nSpecies) = 0.0_wp ! Initial/default number density for each species

   ! Reactions:

   ! Reaction network (see routine classifyReactionNetwork below)
   ! We need to multiply by 2 for use in IRDME
   integer, dimension(2*nReactions,4) :: reactionTable=0 ! Records the type of reagent and product particles of each reaction
      ! For example, let the third reaction be 'iSpecies + 0 --> jSpecies + kSpecies', then we set:
      ! reactionType(3)=R_REPLACEMENT
      ! reactionTable(3) to be (i,0,j,k).
      ! For IRME we duplicate A+B->? also as B+A->? with half the rate and schedule them independently
   integer, dimension(nReactions*2) :: reactionType = 0 ! Classification of all reactions

   ! We assume here a classical law of mass action and the rates are the rates appearing in the LMA
   ! The LMA is expressed here in *number* densities
   ! Note that the units of the number densities are always 1/m^3 even in 1D and 2D
   ! For IRDME we assume here that the system is well-mixed
   ! and later we convert these to Poisson rates (per unit time) for binary reactions
   real (wp) :: reactionRate(nReactions) = 0.0_wp
   real (wp) :: IRDME_reactionRate(2*nReactions) = 0.0_wp ! Duplicate copy of above for binary reactions to distinguish A+B-> from B+A->

   integer :: totalReactions = nReactions ! Total number of reactions including duplicates for IRDME
   ! Stochiometric coefficients for products and reactants separately (nu_plus and nu_minus):
   integer, dimension(nSpecies,2,nReactions) :: reactionNetwork = 0 
   integer :: reactionScheme = IRDME ! 0 for RDME, 1 for IRDME
   
   ! Doi parameters:
   real (wp) :: DoiCellLength (nMaxDimensions) = 1.0_wp
   real (wp) :: DoiCellVolume
   real (wp) :: fractionExtraParticles = 0.1_wp ! Add 10% more particles
   integer   :: nMaxParticlesPerCell=0 ! This will be recomputed later if less than or equal to zero
   real (wp) :: inputTimestep = 0.0_wp
   real (wp), dimension(nSpecies) :: diffusiveCLF = 0.0_wp

   integer :: problem_type=0 ! 0=uniform, 1=half A's, half B's (interface), 2=stripe

   logical :: usePoisson = .true. ! Use Poisson or Gaussian RNGs
   logical :: addDensityFluctuations = .true. ! Initialize and fill reservoir cells with fluctuations or without?
   logical :: strangSplitting = .true. ! It is more accurate to use Strang splitting here
   logical :: randomShift = .true. ! Randomly shift the collision grid each time step
   integer :: diffuseByHopping = 0 ! If 0, diffusion is done by continuous random walk
                                   ! If 1, diffusion is done by random hopping but positions are not on lattice (RDME or IRDME)
                                   ! If 2, positions remain strictly on a lattice at the center of the cells for RDME only!
   logical :: writeOutput = .true. ! Is this an I/O process?

   ! For sampling the number densities on a coarser grid:
   ! Note: nBlockingSample/nBlockingCollisions = DoiCellLength / sampleCellLength
   ! That is, nBlockingSample sampling/hydro cells equal nBlockingCollisions collision cells
   ! Sampling/hydrodynamic cell parameters:
   integer, dimension (nMaxDimensions) :: nSampleCells = 1
   real (wp), dimension (nMaxDimensions) :: sampleCellLength =1.0_wp
   real (wp) :: sampleCellVolume
   
   ! Global system dimensions:
   real (wp), dimension (1:nMaxDimensions) :: domainLength=1.0_wp
   integer, dimension(nMaxDimensions) :: nBlockingSample=1, nBlockingCollision=1
   real (wp) :: domainVolume
   integer, dimension (0:nMaxDimensions) :: nDoiCells = 1
   
   ! Global boundary conditions on the whole system  
   type (WallType) :: wallBCs (2, nDimensions) ! BCs for non-periodic boundaries
   integer :: reservoirThickness = 1 ! Number of reservoir hydro cells to use
   
contains

!--------------------------------------
! Initialization routines for the whole Doi domain (global)
! These need to be called once per processor per Doi domain (level)

subroutine readDoiParameters (fileUnit) ! Read parameters from a NAMELIST file
   integer, intent(in), optional :: fileUnit ! If present, the unit to read the namelist from

   integer :: nameListFile = 37

   namelist / DoiBoxOptions / reactionRate, reactionNetwork, reactionScheme, numberDensity, &
      speciesDiffusivity, speciesDiameter, sampleCellLength, nSampleCells, &
      nBlockingSample, nBlockingCollision, fractionExtraParticles, &
      usePoisson, wallBCs, addDensityFluctuations, reservoirThickness, &
      inputTimestep, strangSplitting, randomShift, &
      nMaxParticlesPerCell, diffuseByHopping, problem_type
   
   if (.not.present(fileUnit)) then ! Default
      if(writeOutput) write(*,*) "Reading Doi parameters from namelist file"
      open (nameListFile, file=  "DoiBoxOptions.nml", status="old", action="read")
   else 
      nameListFile=fileUnit
   end if
   read (nameListFile, nml=DoiBoxOptions)
   if (.not.present(fileUnit)) close (nameListFile)
      
   ! Kishore/Donev: Read cross-linker parameters
   call readCLParameters(fileUnit)

   if(any(nSampleCells(nDimensions+1:nMaxDimensions)/=1)) stop "Must have only one sampling cell in all trivial dimensions"

end subroutine readDoiParameters


! Prepare all the required values for the whole Doi domain (once per processor)
subroutine initializeDoiParameters (writeToOutput)
   logical, intent(in), optional :: writeToOutput
   ! Before calling this, one must set sampleCellLength, nSampleCells(1:nMaxDimensions), the blocking widths
   if(present(writeToOutput)) writeOutput=writeToOutput

   domainLength = nSampleCells * sampleCellLength
   domainVolume = product(domainLength)
   sampleCellVolume = product (sampleCellLength)

   DoiCellLength = real(nBlockingSample,wp) / real(nBlockingCollision,wp) * sampleCellLength
   DoiCellVolume = product(DoiCellLength)

   nDoiCells(1:nMaxDimensions) = nint (domainLength / DoiCellLength)
   if(any(nDoiCells(nDimensions+1:nMaxDimensions)/=1)) stop "Must have a single particle cell in trivial dimensions"
   nDoiCells(0) = product (nDoiCells(1:nMaxDimensions))

   ! CFL is only meaningful for the non-trivial dimensions
   diffusiveCLF = speciesDiffusivity * inputTimestep / (minval(DoiCellLength(1:nDimensions)**2))
   if(writeOutput) write(*,*) 'The diffusive CLF number for each species is ', diffusiveCLF, " max=", maxval(diffusiveCLF)
   if(writeOutput) write(*,*) 'Domain covered by ', nDoiCells(1:nDimensions), ' or a total of ', nDoiCells(0), &
      ' collision cells of length ', DoiCellLength
   if(writeOutput) write(*,*) 'The domain has volume', domainVolume, 'and each cell has volume', DoiCellVolume
   if(writeOutPUT) write(*,*) 'Sample cell length:', sampleCellLength
   if(writeOutput) write(*,*) 'Doi timestep = ', inputTimestep

   call initializeCLs()    ! Kishore: To initialize cross linker parameters.

end subroutine initializeDoiParameters

!--------------------------------------

!--------------------------------------
! Initialization routines for a Doi box (local)
! These are called after the Doi domain is initialized, once per computational box

! Allocate arrays associated with one box
subroutine createDoiBox (box, boundaryType, lbSample, ubSample)
   type (DoiBox), target, intent (inout) :: box
   integer, intent(in) :: boundaryType (2,nDimensions) ! These must be set externally
   integer, intent(in), dimension(nMaxDimensions) :: lbSample, ubSample ! Bounds of this box
   
   ! Local variables:
   logical :: periodicOnly = .true.
   integer :: i, j, k, index, nCells, iDimension, iSpecies
   integer :: cellLocation(1:nMaxDimensions), cellIndex

   ! Copy the inputs: 
   box%boundaryType = boundaryType
   if (any(box%boundaryType /= PERIODIC)) periodicOnly = .false.
   if (.not. periodicOnly) randomShift = .false.
   if (reactionScheme == IRDME) randomShift = .false.

   ! Initialize reaction types.
   call classifyReactionNetwork(box)
   write (*,*) 'Reactions are'
   select case (reactionScheme)
   case(RDME)
      do i = 1, nReactions
         write (*,*) 'Reaction No.', i, 'is', reactionTable(i,1), '+', reactionTable(i,2), '-->', &
            reactionTable(i,3), '+', reactionTable(i,4), 'reaction type is', reactionType(i), 'reaction rate is', reactionRate(i)
      end do
   case(IRDME)
      do i = 1, totalReactions
         write (*,*) 'Reaction No.', i, 'is', reactionTable(i,1), '+', reactionTable(i,2), '-->', &
            reactionTable(i,3), '+', reactionTable(i,4), 'reaction type is', reactionType(i), &
            'reaction rate is', IRDME_reactionRate(i)
      end do
   case default
      stop 'Wrong Reaction Scheme'
   end select

   ! count reaction numbers
   box%reactionCount = 0
   
   box%lb = nint ( (lbSample - 1) * sampleCellLength  / DoiCellLength + 1 )
   box%ub = nint ( ubSample * sampleCellLength / DoiCellLength )
   box%lb(nDimensions+1:nMaxDimensions) = 1
   box%ub(nDimensions+1:nMaxDimensions) = 1
   box%extent(1:nMaxDimensions) = box%ub - box%lb + 1
   box%extent(0) = product(box%extent(1:nDimensions))

   if(printDebug) then
      write(*,*) "---------------------------------"
      write(*,*) "lbSample    = ", lbSample
      write(*,*) "lbCollision = ",  box%lb
      write(*,*) "ubSample    = ", ubSample
      write(*,*) "ubCollision = ", box%ub
   end if

   ! Calculate boundary locations
   box%boundaryLocation(1,:) = (lbSample(1:nDimensions) - 1) * sampleCellLength(1:nDimensions)
   box%boundaryLocation(2,:) = ubSample(1:nDimensions) * sampleCellLength(1:nDimensions)

   box%lbSample = lbSample
   box%ubSample = ubSample
   box%nHydroCells(1:nMaxDimensions)=(box%ubSample-box%lbSample+1)
   box%nHydroCells(0)=product(box%nHydroCells(1:))

   if(writeOutput) write(*,*) "Creating box with lb=", box%lb, " ub=", box%ub, &
      " lbSample=", box%lbSample, " ubSample=", box%ubSample, &
      " BC=", box%boundaryType

   ! Allocate reaction cell data:
   allocate (box%numberPerCell(0:nSpecies, box%lb(1) :box%ub(1), box%lb(2) :box%ub(2), box%lb(3) :box%ub(3)))
   
   select case(reactionScheme)
   case(RDME)
      allocate (box%speciesPointers(nSpecies, box%lb(1) :box%ub(1), box%lb(2) :box%ub(2), box%lb(3) :box%ub(3)))
      box%speciesPointersList(1:nSpecies,1:box%extent(0)) => box%speciesPointers
   case(IRDME)
      nCells = box%extent(0)
      call initializeHeap(box%heap,nCells)
      allocate (box%numberPerCellList(0:nSpecies,1:box%extent(0)))
      allocate (box%firstParticleList(nSpecies, nCells)) ! For LLCs
      ! This will be used to pre-compute and store the list of neighbors
      allocate (box%neighboringCells(neighborhoodSize(0),box%extent(0)))
      ! Alias a pointer to the same storage using Fortran 2003 pointer rank-remapping facility:
      box%firstParticle(1:nSpecies, box%lb(1) :box%ub(1), box%lb(2) :box%ub(2), box%lb(3) :box%ub(3)) => box%firstParticleList
      box%numberPerCell(0:nSpecies, box%lb(1) :box%ub(1), box%lb(2) :box%ub(2), box%lb(3) :box%ub(3)) => box%numberPerCellList
   case default
      stop "reactionScheme can only be one of 0 (RDME) or 1 (IRDME)"
   end select
   
   ! We calculate these for each species, and use index 0 for the total:
   allocate (box%numberDensity(box%lbSample(1) :box%ubSample(1), box%lbSample(2) :box%ubSample(2), &
      box%lbSample(3) :box%ubSample(3), 0:nSpecies))
         
end subroutine createDoiBox



subroutine initializeDoiBox (box) ! Initialize with uniform equilibrium values
   type (DoiBox), target, intent (inout) :: box

   integer, dimension(nMaxDimensions) :: iCell
   integer :: nParticles, specie,i,j,k,nCells, particle, iParticle
   real(wp), dimension(nSpecies) :: numberDensityCopy
   real(wp) :: stripe_width

   ! Estimate the average number of particles in each species:
   do specie=1, nSpecies
      if(problem_type==3) then
         box%nParticles (specie) = nint (numberDensity(specie)*sampleCellVolume)
      else
         box%nParticles (specie) = nint (numberDensity(specie)*domainVolume)      
      end if   
      if(writeOutput) write (*,*) "Estimated number of particles of species ", specie, " is ", box%nParticles(specie)
   end do
   box%nParticles(0) = sum(box%nParticles(1:nSpecies)) ! Total number of particles

   nParticles = ceiling ((1+fractionExtraParticles)*box%nParticles (0))
   if(writeOutput) write(*,*) "Estimated total number of particles = ", box%nParticles(0), " allocated=", nParticles
   box%nMaxParticles = nParticles
   
   if(nMaxParticlesPerCell<=0) then
      nMaxParticlesPerCell = ceiling ( (7.0_wp * nParticles) / nDoiCells(0) )
      write(*,*) "Setting nMaxParticlesPerCell=", nMaxParticlesPerCell
   end if

   allocate ( box%particle(nParticles) )
   if (reactionScheme==RDME) then
      allocate ( box%particlesSortedByCell(nParticles) )
      box%freeParticle=1
   else ! We use an LLC data structure
      !allocate ( box%previousParticle(nParticles) )
      allocate ( box%nextParticle(nParticles) )
      box%freeParticle=1 ! Start of list of free particles
      do particle=1, box%nMaxParticles-1
         box%nextParticle(particle)=particle+1
      end do
      box%nextParticle(box%nMaxParticles)=0 ! End of list of free particles
   end if
   
   ! Now actually initialize the domain with particles
   do iParticle = 1, nParticles
      box%particle(iParticle)%position = 0.0_wp
   end do
   
   write(*,*) "Starting to fill domain with particles"
   
   ! PARALLEL: Potential to parallelize this loop as operation is purely local
   !    However, this is only done once for initialization so benefit is minor
   ! We fill cell by cell:
   box%nParticles=0
   do k = box%lbSample(3), box%ubSample(3)
   do j = box%lbSample(2), box%ubSample(2)
   do i = box%lbSample(1), box%ubSample(1)
      iCell = (/ i, j, k/)
      
      select case(problem_type)
      case(1) ! Generate a biased configuration where all A's are in half the domain and B in the other half
         ! Donev: Warning, this assumes sampling cells 
         numberDensityCopy=numberDensity
         if( (i-box%lbSample(1) <= box%nHydroCells(1)/3) .or. (i-box%lbSample(1) >= 2*box%nHydroCells(1)/3) ) then
            numberDensityCopy(2)=0.0_wp
         else
            numberDensityCopy(1)=0.0_wp
         end if   
         call fillSamplingCell(box, iCell, numberDensityCopy)
      case(2) ! Generate a biased configuration where all particles are in a narrow stripe
         stripe_width = 1.0d0/16.0d0 ! Adjust this to match needs
         if( (i-box%lbSample(1) <  (1-stripe_width)*box%nHydroCells(1)/2) .or. &
             (i-box%lbSample(1) >= (1+stripe_width)*box%nHydroCells(1)/2) ) then
            numberDensityCopy=0.0_wp
         else
            numberDensityCopy=numberDensity
         end if   
         call fillSamplingCell(box, iCell, numberDensityCopy)         
      case(3) ! Only put particles in the center cell -- for detailed testing and comparison to FHD
         if(all(iCell==1+(box%nHydroCells(1:nMaxDimensions)-1)/2)) then
            numberDensityCopy=numberDensity
         else
            numberDensityCopy=0.0_wp
         end if               
         call fillSamplingCell(box, iCell, numberDensityCopy)         
      case default
         call fillSamplingCell(box, iCell, numberDensity)
      end select
      
   end do
   end do
   end do
   if (reactionScheme==IRDME) then
      call initializeIRDME(box)
   else
      call sorter(box)
   end if   
   write (*,*) "Total number of particles at the beginning is ", box%nParticles(0), "=", box%nParticles(1:nSpecies)
   
   if (isAssertsOn) then
      box%initialParticleCount = box%nParticles
   end if

end subroutine initializeDoiBox



subroutine classifyReactionNetwork(box)
   type (DoiBox), target, intent(inout) :: box
   integer :: iReaction, iSpecies, reagentCount, productCount, &
      firstReagent, secondReagent, firstProduct, secondProduct, indexReaction, IRDME_index
   real (wp) :: thickness, effectiveArea, factor, reactiveRadius
   logical :: catalyst

   IRDME_index = 0
   IRDME_reactionRate = 0
   indexReaction = 0
   do iReaction = 1, nReactions
      reagentCount = 0
      productCount = 0
      do iSpecies = 1, nSpecies
         if (reactionNetwork(iSpecies,1,iReaction)>0) then
            reagentCount = reagentCount + reactionNetwork(iSpecies,1,iReaction)
         end if
         if (reactionNetwork(iSpecies,2,iReaction)>0) then
            productCount = productCount + reactionNetwork(iSpecies,2,iReaction)
         end if
      end do

      if (reagentCount <0 .or. reagentCount >2) stop 'Unsupported Reagent Number'
      if (productCount <0 .or. productCount >2) stop 'Unsupported Product Number'

      firstReagent = 0; secondReagent= 0
      firstProduct = 0; secondProduct= 0
      do iSpecies = 1, nSpecies
         if (reactionNetwork(iSpecies,1,iReaction) > 0) then
            if (reactionNetwork(iSpecies,1,iReaction)==2) then
               firstReagent = iSpecies
               secondReagent= iSpecies
               cycle
            else if (firstReagent == 0) then
               firstReagent = iSpecies
            else
               secondReagent = iSpecies
            end if
         end if
      end do
      do iSpecies = 1, nSpecies
         if (reactionNetwork(iSpecies,2,iReaction) > 0) then
            if (reactionNetwork(iSpecies,2,iReaction)==2) then
               firstProduct = iSpecies
               secondProduct= iSpecies
               cycle
            else if (firstProduct == 0) then
               firstProduct = iSpecies
            else
               secondProduct = iSpecies
            end if
         end if
      end do
      indexReaction = indexReaction + 1
      reactionTable(indexReaction,1)=firstReagent
      reactionTable(indexReaction,2)=secondReagent
      reactionTable(indexReaction,3)=firstProduct
      reactionTable(indexReaction,4)=secondProduct
      call classifyReaction()
      if (reactionScheme == IRDME) then
      
         ! We convert the input rates here to have units of inverse time (Poisson rates)
         ! This assumes a well-mixed system and is only an approximation
         ! But it ensures that we can use the same input files for RDME and IRDME
         if (secondReagent/=0) then
            reactiveRadius = (speciesDiameter(firstReagent)+speciesDiameter(secondReagent))/2
            if (nDimensions==3) then
               thickness = 1.0_wp
               effectiveArea = 4.0_wp/3.0_wp*PI*(reactiveRadius**3)
            else if (nDimensions==2) then
               thickness = domainLength(3)
               effectiveArea = PI*(reactiveRadius**2)
            else if (nDimensions==1) then
               thickness = domainLength(2)*domainLength(3)
               effectiveArea = 2.0_wp*reactiveRadius
            end if
            factor = 1.0_wp/(thickness*effectiveArea)
            write(*,*) "IRDME: Multiplying input reactionRate for reaction ", iReaction, " by factor=", factor
         else
            factor = 1.0_wp            
         end if
         
         ! For binary reactions in IRDME with different species, we need to split into A+B-> and B+A->
         ! For A+B we have k = lambda * n_A*n_B * reactive_volume
         ! But for A+A we have k = lambda/2 * n_A^2 * reactive_volume since pairs are counted twice
         !    -- this is why we multiply the rate below by 2
         ! What we compute below is the microscopic rate "lambda" in the Doi model
         ! Since each pair of particles can be selected twice for both for A+A-> and A+B->,
         ! regardless of whether the two reacting particles are in the same cell or in different cells,
         ! we multiply the rate by 1/2 later in calculatePropensity_IRDME
         if ((reagentCount==2) .and. (firstReagent/=secondReagent)) then ! A+B->?
            IRDME_reactionRate(indexReaction) = reactionRate(iReaction)*factor
            ! Now also duplicate as the reverse reaction
            indexReaction = indexReaction + 1
            reactionTable(indexReaction,1)=secondReagent
            reactionTable(indexReaction,2)=firstReagent
            reactionTable(indexReaction,3)=firstProduct
            reactionTable(indexReaction,4)=secondProduct
            call classifyReaction()
            IRDME_reactionRate(indexReaction) = reactionRate(iReaction)*factor
         else if(reagentCount==2) then ! A+A->?
            IRDME_reactionRate(indexReaction) = 2*reactionRate(iReaction)*factor
         else ! Here factor=1 but let us leave it in for future use
            IRDME_reactionRate(indexReaction) = reactionRate(iReaction)*factor   
         end if
      end if
   end do

   totalReactions = indexReaction

contains

   subroutine classifyReaction()

      select case(reagentCount)

      case(0)
      reactionType(indexReaction)=R_BIRTH

      case(1)
      if (productCount==0) then
         reactionType(indexReaction)=R_ANNIHILATION
      else if (productCount==1) then
         reactionType(indexReaction)=R_REPLACEMENT
      else if (productCount==2) then
         if ((firstReagent==firstProduct).or.(firstReagent==secondProduct)) then
            reactionType(indexReaction) = R_CATA_BIRTH
         else
            reactionType(indexReaction)=R_REPLACEMENT
         end if
      else
         stop 'Wrong reaction type.'
      end if

      case(2)
      catalyst = .false.
      
      if ((firstReagent == firstProduct) .or. (firstReagent == secondProduct) &
         .or. (secondReagent == firstProduct) .or. (secondReagent == secondProduct)) catalyst = .true.
      if (productCount==0) then
         reactionType(indexReaction)=R_ANNIHILATION
      else if (productCount==1) then
         if (catalyst) then
            reactionType(indexReaction) = R_CATA_ANNIHILATION
         else
            reactionType(indexReaction) = R_MERGE
         end if
      else if (productCount==2) then
         ! Even though the reaction A+A->A+B is catalytic, we treat it as one with type TRANSFORM since
         ! they share the same way of position changing.
         if ((.not. catalyst) .or. ((firstReagent==secondReagent) .or. (firstProduct==secondProduct))) then
            reactionType(indexReaction) = R_TRANSFORM
         else
            reactionType(indexReaction) = R_CATALYST
         end if
      else
         stop 'Unsupported reaction type'
      end if

      case default
         stop 'Wrong reaction type: wrong reagent number'
      end select
   end subroutine classifyReaction

end subroutine classifyReactionNetwork

            
! This is the subroutine that calls mover, sorter, and reactor.
! It distinguishes between the case where we need to do strang splitting and the other case.
subroutine updateDoiBox(box,timestep)
   type (DoiBox), target :: box
   real (wp), intent(in) :: timestep
   real (wp), dimension(nDimensions) :: randomShiftDist

   integer :: consumed, produced, temp, difference, iSpecies, iReaction, iParticle, nParticles(0:nSpecies)
   integer, dimension(nMaxDimensions) :: cellIndices

   call fillReservoirCells(box) ! Donev: This is a no-op at present since only periodic BCs work

   ! In general we want to do strang splitting even for IRDME since it is more accurate than
   ! move by dt, react by dt
   ! which is called Lie splitting
   ! http://www.staff.science.uu.nl/~frank011/Classes/numwisk/ch13.pdf

   if(randomShift.and.(reactionScheme == RDME)) then
      ! Randomly displace all particles by up to half a grid cell each time step to minimize spatial artifacts
      call UniformRNGVec(randomShiftDist,nDimensions)
      randomShiftDist = randomShiftDist - 0.5_wp
   else
      randomShiftDist = 0.0_wp
   end if

   if (strangSplitting) then
      call mover(box,0.5_wp*timestep,randomShiftDist)
      ! call cleanParticles(box) ! Do this now even if not necessary      
      if (nReactions>0) call reactor(box,timestep)
      call mover(box,0.5_wp*timestep,-randomShiftDist)
      call cleanParticles(box) ! Remove empty spaces
   else ! No strang splitting below

      call mover(box,timestep,randomShiftDist)
      !call cleanParticles(box) ! Do this now even if not necessary      
      if (nReactions>0) call reactor(box,timestep)
      call mover(box,0.0_wp,-randomShiftDist)
      call cleanParticles(box) ! Remove empty spaces

   end if
   
   !box%globalTime = box%globalTime + timestep ! Looses accuracy in single precision
   box%nTimesteps = box%nTimesteps + 1
   box%globalTime = box%nTimesteps * timestep
   
   ! Checks if the number of particle change matches with reaction count.
   if(isAssertsOn) then
      nParticles=box%nParticles
      box%nParticles=0
      do iParticle = lbound(box%particle,1), ubound(box%particle,1)
         iSpecies = box%particle(iParticle)%species
         if (iSpecies <= 0) cycle
         box%nParticles(0) = box%nParticles(0) + 1
         box%nParticles(iSpecies) = box%nParticles(iSpecies) + 1
      end do
      if(any(box%nParticles/=nParticles)) stop "Wrong count of particles at end of time step"
   end if
   
   if (isAssertsOn.and.(.not.IRDME_test)) then
      call testPositionInFakeDimension(box)
      do iSpecies = 1, nSpecies
         difference = box%nParticles(iSpecies)-box%initialParticleCount(iSpecies)
         temp = 0
         do iReaction = 1, totalReactions
            consumed = 0
            produced = 0
            if (reactionTable(iReaction,1)==iSpecies) consumed = consumed + 1
            if (reactionTable(iReaction,2)==iSpecies) consumed = consumed + 1
            if (reactionTable(iReaction,3)==iSpecies) produced = produced + 1
            if (reactionTable(iReaction,4)==iSpecies) produced = produced + 1
            temp = temp - box%reactionCount(iReaction)*consumed + box%reactionCount(iReaction)*produced
         end do
         if ((all(wallBCs%wallType==PERIODIC)).and.(temp /= difference)) then
            write (*,*) 'Difference is ', difference, 'Temp is ', temp
            write (*,*) box%reactionCount
            stop 'Reaction count does not match particle count'
         end if
      end do
      !if (all(wallBCs%wallType==PERIODIC)) write (*,*) 'Particle count matches reaction count.'
   end if
   
end subroutine updateDoiBox

!--------------------------------------
! Computational Doi routines (advection, diffusion, collisions, etc.)
!--------------------------------------

! Fill a given sampling (hydrodynamic) cell with particles
subroutine fillSamplingCell(box, iCell, cellNumberDensity)
   type (DoiBox), target, intent (inout) :: box
   integer, intent (in) :: iCell(nMaxDimensions) ! The position of the sampling cell in the grid
   real (wp), intent (in) :: cellNumberDensity(nSpecies) ! Desired values of number densities

   real (wp), dimension (nDimensions) :: random ! Donev: I changed this to nDimensions
   real (wp) :: rng, mean
   integer :: nParticlesInCell, iSpecies, iParticle, iteration
   logical :: isReservoirCell
   
   do iSpecies = 1, nSpecies       
   
      ! Decide how many particles to generate:
      if(addDensityFluctuations) then ! Generate a Poisson number of particles for each species; Kishore: True right now.
         mean = sampleCellVolume * cellNumberDensity(iSpecies)
         if(usePoisson) then
            call PoissonRNG (nParticlesInCell, mean)
         else ! Approximate with a Gaussian
            call NormalRNG(rng)
            nParticlesInCell = nint(mean + sqrt(mean)*rng)
         end if
      else if(all(nSampleCells==1)) then ! Fix the density at a constant (nearest integer)
         nParticlesInCell = nint( sampleCellVolume * cellNumberDensity(iSpecies))
      else ! Smart rounding to get desired density on average but with minimal fluctuations
         call UniformRNG(random(1))
         nParticlesInCell = floor( sampleCellVolume * cellNumberDensity(iSpecies) + random(1) )
      end if
      isReservoirCell = ( any(iCell < 1) .or. any(iCell(1:nDimensions) > nSampleCells(1:nDimensions)) )

      if(.false.) write(*,*) "Inserting ", nParticlesInCell, " particles of species ", iSpecies, &
         " in cell ", iCell, " isR=", isReservoirCell

      do iteration = 1, nParticlesInCell
         iParticle = box%freeParticle   ! Kishore: This is consistently 1, so later when we index by iParticle
         ! What is the purpose of this?

         if(diffuseByHopping>1) then ! Remain strictly on a lattice
            random=0.5_wp
         else   
            call UniformRNGVec (random, size(random))
         end if   
         ! Kishore: Is this what I would have to change when taking in Ondrej's actin?
         ! I think that this is where I would change the initialization with ondrej's dimers.
         box%particle(iParticle)%position(1:nDimensions) = ( iCell(1:nDimensions) - 1 + random ) * sampleCellLength(1:nDimensions)
         box%particle(iParticle)%position(nDimensions+1:nMaxDimensions) = 0.0_wp

         if (isReservoirCell) then ! Negative species for particles in reservoir cells
            box%particle(iParticle)%species = -iSpecies
         else
            box%particle(iParticle)%species = iSpecies
         end if
         
         ! Update both total species count and individual species count
         box%nParticles(0) = box%nParticles(0) + 1
         box%nParticles(iSpecies) = box%nParticles(iSpecies) + 1
         call updateFreeParticle(box)
         
      end do
   end do

end subroutine fillSamplingCell


! Deletes and repacks particles
subroutine cleanParticles (box)
   type (DoiBox), target :: box
   integer :: p, lastParticle
   
   if(reactionScheme==IRDME) return

   !write(*,*) "Starting search from ", box%freeParticle, count(box%particle(:)%species/=0)
   lastParticle = box%freeParticle - 1
   p = lbound(box%particle,1)

   PackParticles: do

      if (box%particle(p)%species > 0) then
         if(p>=lastParticle) exit ! The array has been packed 
      else  

         ! Find the last particle with a positive species
         do
            if(lastParticle<lbound(box%particle,1)) exit PackParticles ! No particles left
            if(box%particle(lastParticle)%species > 0) exit ! Found a particle
            box%particle(lastParticle)%species = 0 ! Delete any reservoir particles
            lastParticle = lastParticle - 1
         end do
         if(p>=lastParticle) exit ! The array has been packed
         box%particle(p) = box%particle(lastParticle)      
         box%particle(lastParticle)%species = 0
         lastParticle = lastParticle - 1

      end if

      p = p + 1
   end do PackParticles
   box%freeParticle = lastParticle + 1
end subroutine cleanParticles


subroutine fillReservoirCells (box) ! Donev: Not really used at present since only periodic BCs work
   type (DoiBox), target :: box
   integer :: i, j, k, iCell(nMaxDimensions)
   
   if (any(box%boundaryType == RESERVOIR)) then
      if (reactionScheme == IRDME) stop "Reservoir cells not supported for IRDME"
      
      do k = box%lbSample(3), box%ubSample(3)
      do j = box%lbSample(2), box%ubSample(2)
      do i = box%lbSample(1), box%ubSample(1)
         iCell = (/ i, j, k/)
         if ( any(iCell < 1) .or. any(iCell(1:nDimensions) > nSampleCells(1:nDimensions)) ) then
            call fillSamplingCell(box, iCell, numberDensity)
         end if
      end do
      end do
      end do
   end if

end subroutine fillReservoirCells

! Move all particles by dt, taking into account boundary conditions
subroutine mover (box,timestep,randomShift)
   type (DoiBox), target :: box
   real (wp), intent(in) :: timestep
   real (wp), optional, intent(in) :: randomShift(nDimensions)
   integer :: iDimension, p, iSpecies
   integer :: initialSpecies
   real (wp) :: dtime
   real (wp), dimension (nDimensions) :: initialPosition
   logical, parameter :: evolve_r_cm = .true., test_actin = .true.
   !integer :: nLeaving, nEntering ! Donev: These are not really useful

   if (IRDME_trace) return
   if(.false.) then ! Debug the reservoir boundaries
      write(*,*) "Initial n_particles: real=", count(box%particle(:)%species>0), &
         " reservoir=", count(box%particle(:)%species<0)
   end if
   dtime = timestep
   ! QY: If the timestep is 0, then we are actually just doing the second half of random shift in the case in the
   ! case where there is no strang splitting. Hence, we don't need to call brownianMover at all.
   if (timestep>0.0_wp) call brownianMover(box,dtime) ! Move according to simple Brownian motion

   ! PARALLEL: This is a loop that can benefit from parallelization
   !    It is a loop over particles, some of which can be no-ops
   !    So it is important here to use cyclic distribution over threads
   !nLeaving=0; nEntering=0
   do p = lbound(box%particle,1), ubound(box%particle,1)
   
      iSpecies = abs(box%particle(p)%species)
      if (iSpecies == 0) cycle
      
      ! QY: Here we do the random shift. The second way of doing random shift, i.e., by shifting the boundaries, is harder to do
      if (present(randomShift)) box%particle(p)%position(1:nDimensions) = &
         box%particle(p)%position(1:nDimensions) + randomShift * DoiCellLength(1:nDimensions)
      initialSpecies = box%particle(p)%species

      ! Note that due to roundoff error it is possible for a particle to end exactly on the boundary after the position is flushed to memory
      ! This is particularly true for single precision
      do iDimension = 1, nDimensions
         if (box%particle(p)%position(iDimension) < box%boundaryLocation(1,iDimension)) then
            select case(box%boundarytype(1,iDimension))

            case(PERIODIC)
               do while (box%particle(p)%position(iDimension) < box%boundaryLocation(1,iDimension))
                  box%particle(p)%position(iDimension) = box%particle(p)%position(iDimension)+domainLength(iDimension)
               end do
            case(RESERVOIR)
               if (initialSpecies > 0) then
                  !nLeaving = nLeaving + 1
                  box%particle(p)%species = -abs(box%particle(p)%species)
               end if

            case default
               stop 'Boundary condition not implemented yet'

            end select
         else if (box%particle(p)%position(iDimension) >= box%boundaryLocation(2,iDimension)) then
            select case(box%boundarytype(2,iDimension))

            case(PERIODIC)
               do while ((box%particle(p)%position(iDimension) >= box%boundaryLocation(2,iDimension)))
                  box%particle(p)%position(iDimension) = box%particle(p)%position(iDimension)-domainLength(iDimension)
               end do

            case(RESERVOIR)
               if (initialSpecies > 0) then
                  !nLeaving = nLeaving + 1
                  box%particle(p)%species = -abs(box%particle(p)%species)
               end if

            case default
               stop 'Boundary condition not implemented yet'

            end select
         end if

         ! Reservoir particle entering box.  
         if (initialSpecies<0) then
            if (all(box%particle(p)%position(1:nDimensions) > box%boundaryLocation(1,:)) .and. &
                all(box%particle(p)%position(1:nDimensions) < box%boundaryLocation(2,:))) then
               box%particle(p)%species = abs(box%particle(p)%species)
            end if
         end if

      end do
      
      if (isAssertsOn) then
         iSpecies=box%particle(p)%species
         if(iSpecies>0) then ! Positive species must be inside the box
            if ( .not. all( (box%particle(p)%position(1:nDimensions) >= box%boundaryLocation(1,:)) .and. &
                          ( box%particle(p)%position(1:nDimensions) <= box%boundaryLocation(2,:)) )) then
               write (0,*) p, box%particle(p)%species, " Particle outside of box! New r=", &
                  box%particle(p)%species, box%particle(p)%position(1:nDimensions), &
                  " prior r=", initialSpecies, initialPosition, " box low=", box%boundaryLocation(1,:), &
                  " box high=", box%boundaryLocation(2,:)
               stop
            end if
         else if(iSpecies<0) then ! Negative species must be outside of the box
            if ( all( (box%particle(p)%position(1:nDimensions) >= box%boundaryLocation(1,:)) .and. &
                          ( box%particle(p)%position(1:nDimensions) <= box%boundaryLocation(2,:)) )) then
               write (0,*) p, iSpecies," Particle inside of box!", &
                  " New r=", box%particle(p)%species, box%particle(p)%position, &
                  " Old r=", initialSpecies, initialPosition, " box low=", box%boundaryLocation(1,:), &
                  " box high=", box%boundaryLocation(2,:)
               stop
            end if
         end if
      end if

   end do
   if(.false.) then ! Debug the reservoir boundaries
      write(*,*) "Final n_particles: real=", count(box%particle(:)%species>0), &
         " reservoir=", count(box%particle(:)%species<0)
   end if


contains

   subroutine brownianMover(box,dtime)
      type (DoiBox), target :: box
      real (wp), intent(in) :: dtime         ! Kishore: Should I be assuming that this dt is already a 
      ! function of tau (eg. 0.01 * tau)?
      ! Donev: This is whatever it is, you don't control it here (it is set in main, and you don't worry about how it relates to tau (user responsibility))
      
      ! Local variables
      real (wp), dimension(nMaxDimensions) :: disp
      integer :: p, dim, side, k
      real (wp) :: D, probabilities(0:2*nDimensions), r, prob, mu_1, mu_2

      ! PARALLEL: This is the main loop that can benefit from parallelization
      !    It is a loop over particles, some of which can be no-ops
      !    So it is important here to use cyclic distribution over threads
      do p =lbound(box%particle,1),ubound(box%particle,1)
         if (box%particle(p)%species > 0) then
            D = speciesDiffusivity(box%particle(p)%species)
            
            disp = 0.0_wp            
            if(diffuseByHopping>=1) then ! Diffuse by random hops on a lattice   
            
               ! Compute the probability of jumping in each direction:
               k=0  
               probabilities(0) = 0.0_wp
               do dim=1, nDimensions
                  do side=1, 2 ! Left/down or right/up hop
                     k=k+1
                     prob = D*dtime/DoiCellLength(dim)**2
                     probabilities(k) = probabilities(k-1) + prob
                     probabilities(0) = probabilities(0) + prob
                  end do
               end do
               
               probabilities(0) = 1.0_wp - probabilities(0) ! Probability of staying in place
               if(probabilities(0)<0.0_wp) stop "Time step too large (violates CFL condition) for diffusion by hopping"
               
               ! Find a direction to jump in, or stay in place (note that disp=0 by default above)
               call UniformRNG(r)
               k=0
               !write(*,*) "particle=", p, " r=", r, " prob_stay=", probabilities(0)
               FindHop: do dim=1, nDimensions
                  do side=1, 2 ! Left/down or right/up hop
                     k=k+1
                     if(r < probabilities(k)) then ! Hop in this direction
                        disp(dim) = (2*side-3)*DoiCellLength(dim)
                        !write(*,*) "hop=", k, dim, (2*side-3)
                        exit FindHop
                     end if      
                  end do
               end do FindHop

               box%particle(p)%position = box%particle(p)%position + disp
                     
            else ! Continuous random walk  
                     
               ! Kishore: If we have encountered a cross-linker (species 1) and we are odd, then we will diffuse the odd-even PAIR together
               ! What this means is that I must do nothing if p is even. 
               ! Case where we have a dimer that has not bound actin

               ! Cross linkers must come in pairs, but one could be attached to actin (seperate species)
               ! In this case, we would still want to simulate as a spring in some cases. 
               ! So what matters is that EITHER the p+1th or the pth species are of the CL species
               ! If only one is, then it is a bound CL. If neither, then it is a free CL. 
               ! Also only look at odd cases, since otherwise we would double count.
               ! Donev: Just to amke it clearer, please use parenthesis when constructing composite logical expressions, as I did in the line that follows:
               if (add_springs .and. (mod(p,2) == 1) .and. ((box%particle(p)%species == 1) .or. (box%particle(p+1)%species == 1))) then

                  ! Case where one end of the dimer (r1) has bound actin. It is no longer species 1 but the even one is.
                  if (box%particle(p)%species /= 1 .and. box%particle(p+1)%species == 1) then 
                     mu_1 = 0.0_wp
                     mu_2 = mu_2_0

                  ! Case where the other end of dimer (r2) has bound actin. 
                  else if (box%particle(p)%species == 1 .and. box%particle(p+1)%species /= 1) then 
                     
                     mu_1 = mu_1_0
                     mu_2 = 0.0_wp
                     
                  ! Case where we have a free-floating CL. Both ends are species 1. 
                  ! This else clause captures the case ( box%particle(p)%species == 1 .and. box%particle(p+1)%species == 1) 
                  else  
                     mu_1 = mu_1_0
                     mu_2 = mu_2_0
               
                  end if
                  
                  ! Donev: I moved this code inside DiffusionCL.f90, to not polute this code with too much CL stuff and simplify it
                  call moveDimer(dtime, 1, mu_1, mu_2, r_1=box%particle(p)%position, r_2=box%particle(p + 1)%position)
                                    
                  ! Kishore: One thing that is interesting is that if I print p, I get only odd numbers but not consecutively odd
                  ! So I might get 211, 213, 217. 
                  ! Doesn't this mean that the cross-linkers (species 1) are not in consecutive order?
                  !print *, p

                  ! Output tests for debugging purposes:
                  !if (p == 3) then
                  !   open(newunit = myunit1, file = 'output.txt')
                  !      write(myunit1,*) box%particle(p)%position, box%particle(p + 1)%position
                  !   close(myunit1)
                  !end if

               ! If no springs (add_springs is false) then just diffuse normally. 
               ! Also doesn't matter the species, if add_springs is false then 
               ! all species diffuse the same. 
               ! Donev: Slip odd particles of species 1 only when there are springs (notice this looks at p-1 not p+1)
               else if (add_springs .and. (mod(p,2) == 1) .and. ((box%particle(p)%species == 1) .or. (box%particle(p-1)%species == 1))) then
                  ! Do nothing, since particle was already moved above when p-1 was moved
               else if(D > 0) then ! Donev: Removed .and. (.not. add_springs) -- this allows to do ordinary SRBD when add_springs is false
                  call NormalRNGVec(numbers=disp, n_numbers=nDimensions) ! Mean zero and variance one
                  disp = sqrt(2*D*dtime)*disp ! Make the variance be 2*D*time
                  box%particle(p)%position = box%particle(p)%position + disp
               end if

               ! In all 'else' cases, we would do nothing. (D < 0 makes no sense, and D = 0 is immobile)
               ! And if add_springs = .true. but we are NOT odd, we shouldn't
               ! do anything as we would have done that in the odd step. 
               ! So there is no else.
               
            end if
               
            
         end if
      end do
   end subroutine brownianMover

   ! Kishore: This is just what we used before, (havent touched this)
   subroutine brownianMoverOld(box,dtime)
      type (DoiBox), target :: box
      real (wp), intent(in) :: dtime
      
      ! Local variables
      real (wp), dimension(nMaxDimensions) :: disp
      integer :: p, dim, side, k
      real (wp) :: D, probabilities(0:2*nDimensions), r, prob
      
      ! PARALLEL: This is the main loop that can benefit from parallelization
      !    It is a loop over particles, some of which can be no-ops
      !    So it is important here to use cyclic distribution over threads
      do p =lbound(box%particle,1),ubound(box%particle,1)
         if (box%particle(p)%species > 0) then
            D = speciesDiffusivity(box%particle(p)%species)
            
            disp = 0.0_wp            
            if(diffuseByHopping>=1) then ! Diffuse by random hops on a lattice   
            
               ! Compute the probability of jumping in each direction:
               k=0  
               probabilities(0) = 0.0_wp
               do dim=1, nDimensions
                  do side=1, 2 ! Left/down or right/up hop
                     k=k+1
                     prob = D*dtime/DoiCellLength(dim)**2
                     probabilities(k) = probabilities(k-1) + prob
                     probabilities(0) = probabilities(0) + prob
                  end do
               end do
               
               probabilities(0) = 1.0_wp - probabilities(0) ! Probability of staying in place
               if(probabilities(0)<0.0_wp) stop "Time step too large (violates CFL condition) for diffusion by hopping"
               
               ! Find a direction to jump in, or stay in place (note that disp=0 by default above)
               call UniformRNG(r)
               k=0
               !write(*,*) "particle=", p, " r=", r, " prob_stay=", probabilities(0)
               FindHop: do dim=1, nDimensions
                  do side=1, 2 ! Left/down or right/up hop
                     k=k+1
                     if(r < probabilities(k)) then ! Hop in this direction
                        disp(dim) = (2*side-3)*DoiCellLength(dim)
                        !write(*,*) "hop=", k, dim, (2*side-3)
                        exit FindHop
                     end if      
                  end do
               end do FindHop
                     
            else ! Continuous random walk  
                      
               call NormalRNGVec(numbers=disp, n_numbers=nDimensions) ! Mean zero and variance one
               disp = sqrt(2*D*dtime)*disp ! Make the variance be 2*D*time
               
            end if
               
            box%particle(p)%position = box%particle(p)%position + disp
            
         end if
      end do
   end subroutine brownianMoverOld
   
end subroutine mover

subroutine reactor(box,timestep)
   type (DoiBox), target :: box
   real(wp), intent(in) :: timestep
   
   if(reactionScheme==RDME) then
      call sorter(box)
      call reaction_RDME(box,timestep)
   else
      call initializeLLCs(box)
      call reaction_IRDME(box,timestep)
   end if      

end subroutine

!========================================================
! RDME (more correctly, S-BD-RME
!========================================================

subroutine sorter (box)
   type (DoiBox), target :: box

   integer :: cellSum, speciesSumTemp
   integer :: cellIndices (nMaxDimensions), i, j, k, iParticle, iSpecies, positionInOrdering
  
   box%numberPerCell = 0
   box%nParticles=0
   do iParticle = lbound(box%particle,1), ubound(box%particle,1)
      iSpecies = box%particle(iParticle)%species
      if (iSpecies <= 0) cycle
      box%nParticles(0) = box%nParticles(0) + 1
      box%nParticles(iSpecies) = box%nParticles(iSpecies) + 1
      call computeCellIndices(box,iParticle,cellIndices)
      i = cellIndices (1)
      j = cellIndices (2)
      k = cellIndices (3)
      box%numberPerCell (iSpecies, i, j, k) = box%numberPerCell(iSpecies, i, j, k) + 1
      box%numberPerCell (0, i, j, k) = box%numberPerCell(0, i, j, k) + 1
      if (isAssertsOn) call assert (all(box%numberPerCell(0, :, :, :) == sum(box%numberPerCell(1:, :, :, :), 1)), &
         'sorter: numberPerCell not consistent')
   end do
   
   cellSum = 1
   box%speciesSum = 0
   do k = box%lb(3), box%ub(3)
   do j = box%lb(2), box%ub(2)
   do i = box%lb(1), box%ub(1)
      speciesSumTemp = 0
      do iSpecies = 1, nSpecies
         speciesSumTemp = speciesSumTemp + box%numberPerCell(iSpecies, i, j, k)
         box%speciesPointers (iSpecies, i, j, k) = cellSum + speciesSumTemp
         box%speciesSum(iSpecies) = box%speciesSum(iSpecies) + box%numberPerCell(iSpecies, i, j, k)
      end do
      cellSum = cellSum + box%numberPerCell(0, i, j, k)
   end do
   end do
   end do
   
   do iSpecies = 1, nSpecies
      box%speciesSum(0)=box%speciesSum(0)+box%speciesSum(iSpecies)
   end do
   
   do iParticle = lbound(box%particle,1), ubound(box%particle,1)
      iSpecies = box%particle(iParticle)%species
      if (iSpecies <= 0) cycle
      call computeCellIndices(box,iParticle,cellIndices)
      i = cellIndices (1)
      j = cellIndices (2)
      k = cellIndices (3)

      positionInOrdering = box%speciesPointers (iSpecies, i, j, k) - 1
      box%particlesSortedByCell (positionInOrdering) = iParticle
      box%speciesPointers (iSpecies, i, j, k) = box%speciesPointers(iSpecies, i, j, k) - 1
   end do
   if (isAssertsOn) then
      do iParticle = 1, box%nParticles(0)
         if (box%particlesSortedByCell(iParticle)==0) stop 'Particle is not sorted correctly'
         if (box%particle(box%particlesSortedByCell(iParticle))%species == 0) stop 'Particle is not sorted correctly'
      end do
   end if
   
end subroutine sorter


subroutine reaction_RDME(box,timestep)
   type (DoiBox), target :: box
   real(wp), intent(in) :: timestep
   
   integer :: i, j, k
   integer, dimension(nMaxDimensions) :: cellIndices
   
   ! PARALLEL in principle: This loop can be parallelized since the reactions are local
   !   but not important since this is not our goal and we can do this in parallel better in other ways
   ! One must watch out for any global operations in reaction_RDME_cell such as calls to updateFreeParticle
   ! A few of the trouble spots are marked with NOT-PARALLEL
   ! This is a loop over cells so easy to partition, but cyclic is usually best as it works in inhomogeneous cases also
   do k = box%lb(3), box%ub(3)
   do j = box%lb(2), box%ub(2)
   do i = box%lb(1), box%ub(1)
      cellIndices = (/i,j,k/)
      call reaction_RDME_cell(box,timestep,cellIndices)
   end do
   end do
   end do

end subroutine reaction_RDME

subroutine reaction_RDME_cell(box,timestep,cellIndices)
   type (DoiBox), target :: box
   real(wp), intent(in) :: timestep
   integer, dimension(nMaxDimensions), intent(in) :: cellIndices

   real (wp) :: currentTime,tau,r(nDimensions),random,tempR
   integer :: reagent, choiceOfReaction, candidateReactionType,candidate,candidate_1,candidate_2,iParticle
   integer :: temp, temp1,temp2,old,new,iSpecies,typeOfChoice, &
              speciesOfReagent_1,speciesOfReagent_2,speciesOfProduct_1,speciesOfProduct_2, &
              speciesOfSamePosition, speciesOfDifferentPosition, newSpecies
   integer :: choiceReactionTable(1:4)
   real(wp), dimension(0:nReactions) :: propensity
   integer :: i, j, k, index_reaction, n ! iterators
   integer, dimension(0:nSpecies) :: numberPerSpeciesInCell
   integer :: recorder(nSpecies, nMaxParticlesPerCell)

   ! Cell location:
   i=cellIndices(1); j=cellIndices(2); k=cellIndices(3);
      
   numberPerSpeciesInCell = box%numberPerCell(:,i,j,k)
   if(numberPerSpeciesInCell(0)>nMaxParticlesPerCell) stop "Increase nMaxParticlesPerCell"

   ! Make a local copy of the list of particles in this cell
   recorder = 0
   do iSpecies = 1, nSpecies
      do n = 0, numberPerSpeciesInCell(iSpecies)-1
         recorder(iSpecies,1+n)=box%particlesSortedByCell(box%speciesPointers(iSpecies,i,j,k)+n)
      end do
   end do

   currentTime=0
   SSALoop: do

      call calculatePropensity_RDME(propensity)
      ! Ensure the cell is skipped if there's no possible reaction candidate
      if (propensity(0)<=0) exit SSALoop

      ! Sample time of next reaction:
      call calculateTau(propensity(0),tau)
      currentTime=currentTime+tau

      ! Process next reaction
      if (currentTime > timestep) exit SSALoop

      ! Choose reaction to happen next:
      call UniformRNG(random)
      choiceOfReaction = nReactions ! Make sure that this always selects at least one of the reactions even in single precision
      tempR = 0.0_wp
      do index_reaction = 1, nReactions
         if ((tempR <= random).and.(random<tempR+propensity(index_reaction)/propensity(0))) then
            choiceOfReaction = index_reaction
            exit
         end if
         tempR = tempR + (propensity(index_reaction)/propensity(0))
      end do
      typeOfChoice = reactionType(choiceOfReaction)
      ! Keep count
      box%reactionCount(0) = box%reactionCount(0) + 1
      box%reactionCount(choiceOfReaction) = box%reactionCount(choiceOfReaction) + 1

      choiceReactionTable = reactionTable(choiceOfReaction,:)
      speciesOfReagent_1 = choiceReactionTable(1)
      speciesOfReagent_2 = choiceReactionTable(2)
      speciesOfProduct_1 = choiceReactionTable(3)
      speciesOfProduct_2 = choiceReactionTable(4)


      ! Create product particles or remove reacted particles, depending on the type of reaction:
      select case(typeOfChoice)

      ! 1, Birth:    (1)  0 --> A (2) 0 --> A + B (3) 0 --> A + A
      case(R_BIRTH)
      ! Here we generate the position of the products uniformly and randomly inside this cell

      n = box%freeParticle         
      call generateParticlePosition_RDME(box%particle(n)%position(1:nDimensions), cellIndices)
      call createParticle_RDME(n, speciesOfProduct_1)
      call updateFreeParticle(box)

      if (speciesOfProduct_2 /= 0) then
         n = box%freeParticle
         call generateParticlePosition_RDME(box%particle(n)%position(1:nDimensions), cellIndices)
         call createParticle_RDME(n, speciesOfProduct_2)
         call updateFreeParticle(box)
      end if



      ! 2, Annihilation: (1) A --> 0 (2) A + A --> 0 (3) A + B --> 0
      case(R_ANNIHILATION)

      call UniformInteger(candidate,numberPerSpeciesInCell(speciesOfReagent_1))
      call deleteParticle_RDME(candidate, speciesOfReagent_1)
      if (speciesOfReagent_2 /=0 ) then
         call UniformInteger(candidate,numberPerSpeciesInCell(speciesOfReagent_2))
         call deleteParticle_RDME(candidate, speciesOfReagent_2)
      end if



      ! 3, Replacement:  (1) A --> B (2) A --> B + C
      case(R_REPLACEMENT)
      ! Here we change the species of the reacting particle and add a new particle uniformly inside the cell if needed

      call UniformInteger(candidate,numberPerSpeciesInCell(speciesOfReagent_1))
      temp = recorder(speciesOfReagent_1,candidate)

      ! We don't call updateFreeParticle here to avoid wasting one space in the particle array        
      ! Instead, we replace the reactant by the first product particle
      if (speciesOfProduct_2==0) then
         call replaceParticle_RDME(candidate,speciesOfReagent_1,speciesOfProduct_1)
      else ! QY: generate one of the product in position of A and the other randomly in cell.
         call UniformRNG(tempR)
         if (tempR<0.5_wp) then
            speciesOfSamePosition = speciesOfProduct_2
            speciesOfDifferentPosition = speciesOfProduct_1
         else
            speciesOfSamePosition = speciesOfProduct_1
            speciesOfDifferentPosition = speciesOfProduct_2
         end if
         call replaceParticle_RDME(candidate,speciesOfReagent_1,speciesOfSamePosition)
         n = box%freeParticle
         call generateParticlePosition_RDME(box%particle(n)%position(1:nDimensions), cellIndices)
         call createParticle_RDME(n, speciesOfDifferentPosition)
         call updateFreeParticle(box)
      end if

      ! 4, Cata_Birth: (1) A --> A + A (2) A --> A + B
      ! To be consistent, we include A-->A+A here also, even though
      ! it can in principle be treated by the code for replacement.
      ! Here, all we need to do is to create a new particle placed randomly in the cell.
      ! Depending on the products, if it is of type (1), then we create A, if it is of type (2), then we create B.
      case(R_CATA_BIRTH)
      if (speciesOfProduct_1==speciesOfProduct_2) then
         newSpecies = speciesOfProduct_1
      else
         if (speciesOfProduct_1/=speciesOfReagent_1) then
            newSpecies = speciesOfProduct_1
         else
            newSpecies = speciesOfProduct_2
         end if
      end if
      n = box%freeParticle         
      call generateParticlePosition_RDME(box%particle(n)%position(1:nDimensions), cellIndices)
      call createParticle_RDME(n, newSpecies)
      call updateFreeParticle(box)



      ! 5, Merge:       (1) A + A --> B (2) A + B --> C
      case(R_MERGE)
      ! Here we delete one of the two reacting particles and change the species of the remaining one

      ! Decide at whose position the new particle is created
      call UniformRNG(random)
      if (random<0.5_wp) then
         call UniformInteger(candidate,numberPerSpeciesInCell(speciesOfReagent_1))
         call replaceParticle_RDME(candidate,speciesOfReagent_1,speciesOfProduct_1)
         call UniformInteger(candidate,numberPerSpeciesInCell(speciesOfReagent_2))
         call deleteParticle_RDME(candidate,speciesOfReagent_2)
      else
         call UniformInteger(candidate,numberPerSpeciesInCell(speciesOfReagent_2))
         call replaceParticle_RDME(candidate,speciesOfReagent_2,speciesOfProduct_1)
         call UniformInteger(candidate,numberPerSpeciesInCell(speciesOfReagent_1))
         call deleteParticle_RDME(candidate,speciesOfReagent_1)
      end if



      ! 6, Cata_Annihilation: (1) A + B --> A (2) A + A --> A 
      ! In this case B disappears and A remains unchanged
      case(R_CATA_ANNIHILATION)

      if (speciesOfReagent_1==speciesOfProduct_1) then
         call UniformInteger(candidate,numberPerSpeciesInCell(speciesOfReagent_2))
         call deleteParticle_RDME(candidate,speciesOfReagent_2)
      else if (speciesOfReagent_2==speciesOfProduct_1) then
         call UniformInteger(candidate,numberPerSpeciesInCell(speciesOfReagent_1))
         call deleteParticle_RDME(candidate,speciesOfReagent_1)
      else
         stop "catalytic annihilation that has no catylist!"
      end if



      ! 7, Transform:    (1) A + A --> B + C (2) A + B --> C + D (3) A + A --> A + B (4) A + B --> A + A 
      case(R_TRANSFORM)
      ! No new particle positions are generated here, they are simply swapped around

      ! Decide at whose position the new particles are created
      call UniformRNG(random)
      if (random<0.5_wp) then
         call UniformInteger(candidate_1,numberPerSpeciesInCell(speciesOfReagent_1))
         call replaceParticle_RDME(candidate_1,speciesOfReagent_1,speciesOfProduct_1)
         call UniformInteger(candidate_2,numberPerSpeciesInCell(speciesOfReagent_2))
         call replaceParticle_RDME(candidate_2,speciesOfReagent_2,speciesOfProduct_2)
      else
         call UniformInteger(candidate_1,numberPerSpeciesInCell(speciesOfReagent_1))
         call replaceParticle_RDME(candidate_1,speciesOfReagent_1,speciesOfProduct_2)
         call UniformInteger(candidate_2,numberPerSpeciesInCell(speciesOfReagent_2))
         call replaceParticle_RDME(candidate_2,speciesOfReagent_2,speciesOfProduct_1)
      end if


      ! 8, Catalyst:     reaction of form A + B --> A + C
      case(R_CATALYST)
      ! In catalystic reaction A + B --> A + C, the position of catalyst A remains unchanged, so we only need to find out what B and C are
      if (speciesOfReagent_1==speciesOfProduct_1) then
         old = speciesOfReagent_2
         new = speciesOfProduct_2
      else if (speciesOfReagent_1==speciesOfProduct_2) then
         old = speciesOfReagent_2
         new = speciesOfProduct_1
      else if (speciesOfReagent_2==speciesOfProduct_1) then
         old = speciesOfReagent_1
         new = speciesOfProduct_2
      else
         old = speciesOfReagent_1
         new = speciesOfProduct_1
      end if

      ! The following deletes one B and creates one C
      call UniformInteger(candidate,numberPerSpeciesInCell(old))
      call replaceParticle_RDME(candidate,old,new)


      case default

         stop 'Unknown reaction type'

      end select

      if(numberPerSpeciesInCell(0)>nMaxParticlesPerCell) stop "Increase nMaxParticlesPerCell"


   end do SSALoop

   ! Checks if the recorder makes some error.
   if (isAssertsOn) then
      do iSpecies = 1, nSpecies
         do iParticle = 1, numberPerSpeciesInCell(iSpecies)
            if (recorder(iSpecies,iParticle)==0) stop "Problem in recorder"
            if (box%particle(recorder(iSpecies,iParticle))%species /= iSpecies) then
               write (*,*) 'Particle species is ', &
                    box%particle(recorder(iSpecies,iParticle))%species, 'Correct species is ', iSpecies
               write (*,*) iParticle,numberPerSpeciesInCell(iSpecies)
               write (*,*) 'recorder', recorder(iSpecies,iParticle)
               write (*,*) box%nParticles
               stop "Problem in recorder: wrong species. After SSA Loop."
            end if
         end do
      end do
   end if

contains


   subroutine calculatePropensity_RDME(propensity)

      real(wp), intent(out), dimension(0:nReactions) :: propensity
      
      integer :: candidateReactionType, reagent_1, reagent_2, reagentCount

      propensity = 0.0_wp
      ! Compute the reaction rate given stochiometric coefficients
      do candidateReactionType = 1, nReactions
         reagentCount = 0
         reagent_1 = reactionTable(candidateReactionType,1)
         reagent_2 = reactionTable(candidateReactionType,2)
         if (reagent_1/=0) reagentCount = reagentCount + 1
         if (reagent_2/=0) reagentCount = reagentCount + 1
         propensity(candidateReactionType)=reactionRate(candidateReactionType)
         if (reagentCount == 0) then
            propensity(candidateReactionType)=propensity(candidateReactionType)*DoiCellVolume
         else if (reagentCount == 1) then
            propensity(candidateReactionType)=propensity(candidateReactionType)*numberPerSpeciesInCell(reagent_1)
         else
            if (reagent_1 == reagent_2) then ! A+A->something, so skip self-interactions
               propensity(candidateReactionType)=propensity(candidateReactionType)*&
                  numberPerSpeciesInCell(reagent_1)*(numberPerSpeciesInCell(reagent_2)-1)/DoiCellVolume
            else ! A+B->something
               propensity(candidateReactionType)=propensity(candidateReactionType)*&
                  numberPerSpeciesInCell(reagent_1)*numberPerSpeciesInCell(reagent_2)/DoiCellVolume
            end if
         end if
         propensity(0)=propensity(0) + propensity(candidateReactionType)
      end do
   end subroutine calculatePropensity_RDME


   ! Note that the deleteParticle_RDME works slightly different from createParticle_RDME.
   ! Instead of working directly on the box's particle array,
   ! deletaParticle accesses the local recorder first. Hence, the input of these two subroutines are different.
   subroutine deleteParticle_RDME(candidate,species)
      integer, intent(in) :: candidate,species

      if(isAssertsOn) then
         if (candidate==0) stop 'Check RNG'
         if (box%particle(recorder(species,candidate))%species==0) stop 'Deleting nonexistent particles'
      end if
      box%particle(recorder(species,candidate))%species = 0
      
      ! Update the local recorder by moving the last particle here
      recorder(species,candidate)=recorder(species,numberPerSpeciesInCell(species))
      recorder(species,numberPerSpeciesInCell(species))=0
      numberPerSpeciesInCell(species)=numberPerSpeciesInCell(species)-1
      numberPerSpeciesInCell(0)=numberPerSpeciesInCell(0)-1
      
      ! Keep track of this throughout
      ! NOT-PARALLEL: This updates global counters:
      box%nParticles(0)=box%nParticles(0)-1
      box%nParticles(species)=box%nParticles(species)-1
      
   end subroutine
   

   subroutine createParticle_RDME(particle,species)
      integer, intent(in) :: particle,species

      box%particle(particle)%species=species

      ! NOT-PARALLEL: This updates global counters:
      box%nParticles(0)=box%nParticles(0) + 1
      box%nParticles(species)=box%nParticles(species) + 1

      ! Add this to the local recorder:
      numberPerSpeciesInCell(species)=numberPerSpeciesInCell(species)+1
      numberPerSpeciesInCell(0)=numberPerSpeciesInCell(0)+1
      recorder(species,numberPerSpeciesInCell(species))=particle
      
   end subroutine

   ! The routine replaceParticle is equivalent to first calling deleteParticle and then createParticle
   ! It is a bit faster and makes the code more readable, hopefully
   subroutine replaceParticle_RDME(candidate,species_1,species_2)
      integer, intent(in) :: candidate,species_1,species_2

      integer:: particle

      particle = recorder(species_1,candidate)
      if (isAssertsOn) then
         if (particle==0) then
            write (*,*) 'species_1', species_1, 'species_2', species_2, 'candidate', candidate, &
                 'num in sp1', numberPerSpeciesInCell(species_1), &
                 'num in sp2', numberPerSpeciesInCell(species_2)
         end if
      end if
      recorder(species_1,candidate)=recorder(species_1,numberPerSpeciesInCell(species_1))
      recorder(species_1,numberPerSpeciesInCell(species_1))=0
      numberPerSpeciesInCell(species_1)=numberPerSpeciesInCell(species_1)-1

      numberPerSpeciesInCell(species_2)=numberPerSpeciesInCell(species_2)+1
      recorder(species_2,numberPerSpeciesInCell(species_2))=particle

      ! NOT-PARALLEL: This updates global counters:
      box%nParticles(species_1) = box%nParticles(species_1) - 1
      box%nParticles(species_2) = box%nParticles(species_2) + 1
      box%particle(particle)%species = species_2

   end subroutine

   subroutine generateParticlePosition_RDME(position, cellIndices)
      real(kind=wp), intent(out) :: position(1:nDimensions)
      integer, intent(in)  :: cellIndices(1:nDimensions)

      real(kind=wp) :: r(nDimensions)

      if(diffuseByHopping>1) then ! Remain strictly on a lattice
         r=0.5_wp
      else   
         call UniformRNGVec(numbers=r, n_numbers=nDimensions)
      end if   
      
      position = DoiCellLength(1:nDimensions) * (r+(cellIndices(1:nDimensions)-1))

   end subroutine
      
end subroutine reaction_RDME_cell

!========================================================
! IRDME (more correctly, SRBD)
!========================================================

subroutine initializeIRDME(box)
   type (DoiBox), target :: box

   integer :: i,j,k,p,i_2,j_2,k_2,ii,jj,kk,iParticle,iSpecies,iDimension,num,cellIndex,iCell,species,iReaction
   integer :: choiceReactionTable(4)
   integer, dimension(nMaxDimensions) :: cellIndices

   if( any(box%extent(1:nDimensions)<3) ) then
      stop "There must be at least 3 cells in each dimension for IRDME"
   end if
   ! Add an additional check to ensure grid size is reasonable
   do iReaction=1, nReactions      
      choiceReactionTable = reactionTable(iReaction,:)
      if(choiceReactionTable(2)==0) cycle ! Only look at binary reactions
      if (any(DoiCellLength(1:nDimensions) < &
         0.5_wp*(speciesDiameter(choiceReactionTable(1))+speciesDiameter(choiceReactionTable(2))))) then
         write(0,*) "ERROR: grid cells are too small to capture all reactions for reaction #", iReaction, &
            " -- need to be longer than ", 0.5_wp*(speciesDiameter(choiceReactionTable(1))+speciesDiameter(choiceReactionTable(2)))
         stop "Grid cells too small for IRDME"
      end if
   end do
   if(any(box%boundarytype(1:2,1:nDimensions)/=PERIODIC)) stop "Must have periodic BCs for IRDME to work"

   ! Build and store a list of neighboring cells for each cell to avoid recomputing this every time step
   do k = box%lb(3), box%ub(3)
   do j = box%lb(2), box%ub(2)
   do i = box%lb(1), box%ub(1)
      num = 0
      call computeCellIndex(box,(/i,j,k/),iCell)
      do kk = k-neighborhoodSize(3), k+neighborhoodSize(3)
      do jj = j-neighborhoodSize(2), j+neighborhoodSize(2)
      do ii = i-neighborhoodSize(1), i+neighborhoodSize(1)
      ! QY: Here we assume periodic B.C.
         num = num + 1
         i_2 = ii
         if (ii<box%lb(1)) i_2 = box%ub(1)
         if (ii>box%ub(1)) i_2 = box%lb(1)
         j_2 = jj
         if (jj<box%lb(2)) j_2 = box%ub(2)
         if (jj>box%ub(2)) j_2 = box%lb(2)
         k_2 = kk 
         if (kk<box%lb(3)) k_2 = box%ub(3)
         if (kk>box%ub(3)) k_2 = box%lb(3)
         call computeCellIndex(box,(/i_2,j_2,k_2/),cellIndex) 
         box%neighboringCells(num,iCell)=cellIndex
      end do
      end do
      end do
   end do
   end do
   end do
   
   if(IRDME_trace.and.(.not.IRDME_test)) stop "If tracing IRDME then testing must aso be enabled"
   if(IRDME_test) then ! Allocate some counters to count number of reactions in each box and particle
      box%reactionCount = 0
      box%IRDME_test_number = 0
      allocate (box%IRDME_test_number_1(1:box%extent(0)))
      box%IRDME_test_number_1 = 0
      allocate (box%IRDME_test_number_2(1:box%nMaxParticles))
      box%IRDME_test_number_2 = 0
   end if


end subroutine

! QY: The following subroutine initializes the linked list cells (LLCs)
subroutine initializeLLCs (box)
   type (DoiBox), target :: box

   integer :: i,j,k,p,i_2,j_2,k_2,ii,jj,kk,iParticle,iSpecies,iDimension,num,cellIndex,iCell,species,free,freeCount
   integer, dimension(nMaxDimensions) :: cellIndices

   box%firstParticleList = 0
   !box%previousParticle = 0

   ! QY: we need to do the sorting here since sorter is not used in IRDME
   ! Almost PARALLEL: This loop can be parallelized but not trivially
   ! It requires splitting the domain into boxes and sorting the particles into boxes first
   box%numberPerCellList = 0
   box%nParticles = 0
   do iParticle = lbound(box%particle,1),ubound(box%particle,1)
      species = box%particle(iParticle)%species
      if (species>0) then
         box%nParticles(species)=box%nParticles(species)+1
         box%nParticles(0)=box%nParticles(0)+1
         call computeCellIndices(box,iParticle,cellIndices)
         call computeCellIndex(box,cellIndices,cellIndex)
         box%numberPerCellList(species,cellIndex)=box%numberPerCellList(species,cellIndex) + 1
         box%numberPerCellList(0,cellIndex)=box%numberPerCellList(0,cellIndex) + 1
         p = box%firstParticleList(species,cellIndex)
         !box%previousParticle(p)=iParticle
         box%nextParticle(iParticle)=p
         !box%previousParticle(iParticle)=0
         box%firstParticleList(species,cellIndex)=iParticle
      end if
   end do
   
   if (isAssertsOn) then
      do iSpecies = 1, nSpecies
         do k = box%lb(3), box%ub(3)
         do j = box%lb(2), box%ub(2)
         do i = box%lb(1), box%ub(1)
            call computeCellIndex(box,(/i,j,k/),cellIndex)
            if (box%firstParticle(iSpecies,i,j,k)/=box%firstParticleList(iSpecies,cellIndex)) then
               write (*,*) i,j,k
               stop 'Problem in first particle list'
            end if
         end do
         end do
         end do
      end do
      do iSpecies = 0, nSpecies
         do k = box%lb(3), box%ub(3)
         do j = box%lb(2), box%ub(2)
         do i = box%lb(1), box%ub(1)
            call computeCellIndex(box,(/i,j,k/),cellIndex)
            if (box%numberPerCell(iSpecies,i,j,k)/=box%numberPerCellList(iSpecies,cellIndex)) then
               write (*,*) i,j,k
               stop 'Problem in number per cell list'
            end if
         end do
         end do
         end do
      end do
      
      free = box%freeParticle
      freeCount = 0
      do while (free/=0)
         freeCount = freeCount + 1
         if (box%particle(free)%species/=0) then
            write (*,*) 'In initializeLLC: free particle problem: particle', free, 'has species', box%particle(free)%species
            stop 'Check InitializeLLC to find the problem'
         end if
         free = box%nextParticle(free)
      end do
      if (freeCount+box%nParticles(0)/=ubound(box%particle,1)-lbound(box%particle,1)+1) then
         write (*,*) 'Free count is', freeCount
         write (*,*) 'Particle count is', box%nParticles(0)
         write (*,*) 'Total count is', ubound(box%particle,1)-lbound(box%particle,1)+1
         stop 'count mismatch'
      end if
   end if
   
end subroutine initializeLLCs

subroutine reaction_IRDME(box,timestep)
   type (DoiBox), target :: box
   real(wp), intent(in) :: timestep
   
   ! Local variables
   real (wp) :: tau,dist,random,tempR
   integer :: reagent, choiceOfReaction, candidateReactionType,candidate,cellIndex,reactionCell_1, reactionCell_2,length
   integer :: iSpecies,new,typeOfChoice,speciesOfReagent_1,speciesOfReagent_2,speciesOfProduct_1,speciesOfProduct_2,iParticle,temp
   integer :: first, second
   integer, dimension(nMaxDimensions) :: cellLocation,tempCell
   integer :: choiceReactionTable(1:4)
   real (wp) :: propensity(0:totalReactions)
   real (wp), dimension(nMaxDimensions) :: p_position,p_1_position,p_2_position, pos
   integer :: i, j, k, i_2, j_2, k_2, ii, jj, kk, index_reaction, n, n_2, iReaction, iCell ! iterators
   integer :: p, p_prev, p_1, p_1_prev, p_2, p_2_prev, dimens, num
   integer :: nNeighbors(nSpecies)
   logical :: changeInCell_1, changeInCell_2
   
   ! QY: temporary stuff for testing
   logical :: testing
   integer :: testThreshold ! How many reactions to do before starting testing

   testing = .false.
   if(IRDME_trace) then
      ! We turn on testing only after testThreshold reactions have been processed
      testThreshold = -1
      testing = .true.
   else if(IRDME_test) then
      testThreshold = -1 ! Immediately start testing
      testing = .true.
   end if

   ! Build the initial heap:
   !--------------------
   call resetHeap(box%heap)
   box%elapsedTime = 0.0_wp

   ! PARALLEL: This loop can be parallelized but one should first compute the times and store them
   ! and *then* build the queue from them in serial
   ! this requires creating an array to store the tau's in before building the queue
   do iCell = 1, box%extent(0)
      call calculatePropensity_IRDME(box,iCell,propensity,nNeighbors)
      if (propensity(0)>0.0_wp) then
         call calculateTau(propensity(0),tau)
         tau = box%elapsedTime+tau ! box%elapsedTime=0 but add it for code clarity
         if( tau <= timestep) call insertInHeap(box%heap,iCell,tau)
      end if
   end do

   if (isAssertsOn) call LLCtest(box)
   
   ! This loop is very difficult to parallelize
   ! Loop over events:
   !--------------------
   SSALoop : do
      if(box%heap%heapSize <= 0) then ! No more events left to process
         box%elapsedTime = timestep ! Not really needed but to be consistent
         exit SSALoop
      end if
      box%elapsedTime = box%heap%priorityQueue(1)%time ! Jump to this point in time
      changeInCell_1 = .false.
      changeInCell_2 = .false.
      
      if (IRDME_trace) then ! Turn testing on after processing some number of reactions
         if((.not.testing).and.(box%reactionCount(0)>testThreshold)) then
            box%heap%priorityQueue(:)%time = box%heap%priorityQueue(:)%time - box%elapsedTime
            box%elapsedTime = 0.0_wp ! Reset the time to zero to start counting now
            testing = .true.
            !if(isAssertsOn) call testArrays(box)
         end if
      end if
      
      ! Reaction
      reactionCell_1 = box%heap%priorityQueue(1)%elementIndex
      reactionCell_2 = 0
      call calculatePropensity_IRDME(box,reactionCell_1,propensity,nNeighbors)     
      
      ! Choose reaction to happen next:
      call UniformRNG(random)
      choiceOfReaction = totalReactions ! Make sure at least one reaction is always selected
      tempR = 0.0_wp
      do index_reaction = 1, totalReactions
         if ((tempR <= random).and.(random<tempR+propensity(index_reaction)/propensity(0))) then
            choiceOfReaction = index_reaction
            exit
         end if
         tempR = tempR + (propensity(index_reaction)/propensity(0))
      end do
      typeOfChoice = reactionType(choiceOfReaction)

      choiceReactionTable = reactionTable(choiceOfReaction,:)
      speciesOfReagent_1 = choiceReactionTable(1); speciesOfReagent_2 = choiceReactionTable(2)
      speciesOfProduct_1 = choiceReactionTable(3); speciesOfProduct_2 = choiceReactionTable(4)

      ! Create product particles or remove reacted particles, depending on the type of reaction:
      select case(typeOfChoice)


      ! 1, Birth:    (1)  0 --> A (2) 0 --> A + B (3) 0 --> A + A
      ! If testing=false it does the real reactions, otherwise just counts
      ! The code was written carefully to work both in testing mode and in "real" mode
      ! This way we can repeat testing later and in the mean time also use the code to do real runs
      case(R_BIRTH)

      if (testing) then ! Just count this reaction
         box%IRDME_test_number_1(reactionCell_1) = box%IRDME_test_number_1(reactionCell_1)+1
         box%IRDME_test_number(1) = box%IRDME_test_number(1) + 1
         if (speciesOfProduct_2/=0) then
            call UniformRNG(tempR)
            if (tempR<0.5_wp) then
               first = speciesOfProduct_1
               second = speciesOfProduct_2
            else
               first = speciesOfProduct_2
               second = speciesOfProduct_1
            end if
            
            ! Create first particle:
            call UniformRNGVec(numbers=pos, n_numbers=nDimensions)
            n = box%freeParticle
            call revertCellIndex(box,reactionCell_1,cellLocation)
            box%particle(n)%position = DoiCellLength * (pos+real(cellLocation-1,wp))
            call updateFreeParticle(box)
            n_2 = box%freeParticle
            ! Particles born within a sphere of radius equal to sum of particle diameters:
            call getNewPosition(box,n_2,n,speciesDiameter(speciesOfProduct_1)+speciesDiameter(speciesOfProduct_2))
            box%freeParticle = n ! We add this line since in tracing test we don't want to create new particles.
            write(31,"(1000G17.9)") box%particle(n_2)%position-box%particle(n)%position
         end if


      else ! Actually process the reaction
         changeInCell_1 = .true.
         if (speciesOfProduct_2==0) then
            call revertCellIndex(box,reactionCell_1,cellLocation)
            n=box%freeParticle
            call UniformRNGVec(numbers=pos, n_numbers=nMaxDimensions)
            pos(nDimensions+1:nMaxDimensions) = 0.0_wp
            box%particle(n)%position = DoiCellLength * (pos+real(cellLocation-1,wp))
            call createParticle_IRDME(n,speciesOfProduct_1,reactionCell_1)
         else ! Create two particles at once within a reactive distance of each other
            ! The position of the first particle is uniformly distributed inside this cell
            ! while the second one is uniformly distributed in a reactive radius
            changeInCell_2 = .true.
            
            call UniformRNG(tempR)
            if (tempR<0.5_wp) then
               first = speciesOfProduct_1
               second = speciesOfProduct_2
            else
               first = speciesOfProduct_2
               second = speciesOfProduct_1
            end if
            
            ! Create first particle:
            call UniformRNGVec(numbers=pos, n_numbers=nMaxDimensions)
            pos(nDimensions+1:nMaxDimensions) = 0.0_wp
            n = box%freeParticle
            call revertCellIndex(box,reactionCell_1,cellLocation)
            box%particle(n)%position = DoiCellLength * (pos+real(cellLocation-1,wp))
            call createParticle_IRDME(n,first,reactionCell_1)
            
            ! Now create second particle
            n_2 = box%freeParticle

            ! Particles born within a sphere of radius equal to sum of particle diameters:
            call getNewPosition(box,n_2,n,speciesDiameter(speciesOfProduct_1)+speciesDiameter(speciesOfProduct_2))
            call computeCellIndices(box,n_2,cellLocation)
            call computeCellIndex(box,cellLocation,reactionCell_2)

            call createParticle_IRDME(n_2,second,reactionCell_2)
            
         end if
      end if

      box%reactionCount(choiceOfReaction)=box%reactionCount(choiceOfReaction)+1
      box%reactionCount(0)=box%reactionCount(0)+1


      ! 2, Annihilation: (1) A --> 0 (2) A + A --> 0 (3) A + B --> 0
      case(R_ANNIHILATION)

      call UniformInteger(candidate, box%numberPerCellList(speciesOfReagent_1,reactionCell_1))
      p_1 = box%firstParticleList(speciesOfReagent_1,reactionCell_1)
      call findCandidateParticle(box,candidate,p_1,p_1_prev)
      if (speciesOfReagent_2==0) then
         ! Delete directly from the box
         if (testing) then
            box%IRDME_test_number_2(p_1)=box%IRDME_test_number_2(p_1)+1
            box%IRDME_test_number(2) = box%IRDME_test_number(2) + 1
         else
            call deleteParticle_IRDME(p_1,p_1_prev,speciesOfReagent_1,reactionCell_1)
         end if
         box%reactionCount(choiceOfReaction)=box%reactionCount(choiceOfReaction)+1
         box%reactionCount(0)=box%reactionCount(0)+1
         changeInCell_1 = .true.
      else
         call searchNeighbor(box,reactionCell_1,nNeighbors,speciesOfReagent_2,p_2,p_2_prev,reactionCell_2)
         ! Test if they are overlapping and if so, delete from the box
         p_1_position = box%particle(p_1)%position
         p_2_position = box%particle(p_2)%position
         if (p_1/=p_2) then
            if(TestOverlapping(p_1_position,p_2_position, &
                  speciesDiameter(speciesOfReagent_1),speciesDiameter(speciesOfReagent_2))) then

               box%reactionCount(choiceOfReaction)=box%reactionCount(choiceOfReaction)+1
               box%reactionCount(0)=box%reactionCount(0)+1
               
               if (testing) then
                  if (speciesOfReagent_1/=speciesOfReagent_2) then ! we record here the reaction count for A+B->0
                     box%IRDME_test_number(3) = box%IRDME_test_number(3) + 1
                     if(IRDME_trace) then ! Record this pair in the history
                        if(box%IRDME_test_number(3) > nMaxTestPairs) stop "Insufficient room: increase nMaxTestPairs"
                        if (p_1>p_2) then
                           temp = p_1
                           p_1 = p_2
                           p_2 = temp
                        end if
                        box%reactionParticle_AB(box%IRDME_test_number(3),1)=p_1
                        box%reactionParticle_AB(box%IRDME_test_number(3),2)=p_2
                        box%timeOfReactions_AB(box%IRDME_test_number(3)) = box%elapsedTime
                     end if
                  else ! we record here the reaction count for A+A->0
                     box%IRDME_test_number(4) = box%IRDME_test_number(4) + 1
                     if(IRDME_trace) then ! Record this pair in the history
                        if(box%IRDME_test_number(4) > nMaxTestPairs) stop "Insufficient room: increase nMaxTestPairs"
                        if (p_1>p_2) then
                           temp = p_1
                           p_1 = p_2
                           p_2 = temp
                        end if
                        box%reactionParticle_AA(box%IRDME_test_number(4),1)=p_1
                        box%reactionParticle_AA(box%IRDME_test_number(4),2)=p_2
                        box%timeOfReactions_AA(box%IRDME_test_number(4)) = box%elapsedTime
                     end if
                  end if
               else ! Actually annihilate the particles
                  if (box%nextParticle(p_1)==p_2) then
                     call deleteParticle_IRDME(p_2,p_2_prev,speciesOfReagent_2,reactionCell_2)
                     call deleteParticle_IRDME(p_1,p_1_prev,speciesOfReagent_1,reactionCell_1)
                  else
                     call deleteParticle_IRDME(p_1,p_1_prev,speciesOfReagent_1,reactionCell_1)
                     call deleteParticle_IRDME(p_2,p_2_prev,speciesOfReagent_2,reactionCell_2)
                  end if
               end if
               changeInCell_1 = .true.
               changeInCell_2 = .true.
            end if
         end if
      end if


      ! 3, Replacement:  (1) A --> B (2) A --> B + C
      case(R_REPLACEMENT)
      
      ! Choosee a random reagent particle.
      call UniformInteger(candidate,box%numberPerCellList(speciesOfReagent_1,reactionCell_1))
      p = box%firstParticleList(speciesOfReagent_1,reactionCell_1)
      call findCandidateParticle(box,candidate,p,p_prev)
      
      changeInCell_1 = .true.
      if (speciesOfProduct_2 == 0) then
         call replaceParticle_IRDME(p,p_prev,speciesOfReagent_1,speciesOfProduct_1,reactionCell_1)
      else
      
         ! Decide which of the two products the particle A will become
         call UniformRNG(tempR)
         if (tempR<0.5_wp) then
            first = speciesOfProduct_1
            second = speciesOfProduct_2
         else
            first = speciesOfProduct_2
            second = speciesOfProduct_1
         end if
         if (testing) then
            n = box%freeParticle
            call getNewPosition(box,n,p,speciesDiameter(speciesOfProduct_1)+speciesDiameter(speciesOfProduct_2))
            box%count(p) = box%count(p)+1
         else
            ! Turn the A into the new species:
            call replaceParticle_IRDME(p,p_prev,speciesOfReagent_1,first,reactionCell_1)

            ! Now create the second particle
            n = box%freeParticle
            call getNewPosition(box,n,p,speciesDiameter(speciesOfProduct_1)+speciesDiameter(speciesOfProduct_2))
            call computeCellIndices(box,n,cellLocation)
            call computeCellIndex(box,cellLocation,reactionCell_2)
            call createParticle_IRDME(n,second,reactionCell_2)
            changeInCell_2 = .true.
         end if
      end if
      
      box%reactionCount(choiceOfReaction)=box%reactionCount(choiceOfReaction)+1
      box%reactionCount(0)=box%reactionCount(0)+1

      ! 4, Catalitic Birth: (1) A -> A + A (2) A -> A + B
      case(R_CATA_BIRTH)
      ! Choosee a random reagent particle.
      if (speciesOfProduct_1==speciesOfReagent_1) then
         new = speciesOfProduct_2
      else
         new = speciesOfProduct_1
      end if
      call UniformInteger(candidate,box%numberPerCellList(speciesOfReagent_1,reactionCell_1))
      p_1 = box%firstParticleList(speciesOfReagent_1,reactionCell_1)
      call findCandidateParticle(box,candidate,p_1,p_1_prev)
      n = box%freeParticle
      call getNewPosition(box,n,p_1,speciesDiameter(speciesOfProduct_1)+speciesDiameter(speciesOfProduct_2))
      call computeCellIndices(box,n,cellLocation)
      call computeCellIndex(box,cellLocation,reactionCell_2)
      call createParticle_IRDME(n,new,reactionCell_2)
      changeInCell_1 = .true.
      changeInCell_2 = .true.
      
      box%reactionCount(choiceOfReaction)=box%reactionCount(choiceOfReaction)+1
      box%reactionCount(0)=box%reactionCount(0)+1


      ! 5, Merge:       (1) A + A --> B (2) A + B --> C
      case(R_MERGE)

      ! Choosee a random reagent particle.
      call UniformInteger(candidate,box%numberPerCellList(speciesOfReagent_1,reactionCell_1))
      p_1 = box%firstParticleList(speciesOfReagent_1,reactionCell_1)
      call findCandidateParticle(box,candidate,p_1,p_1_prev)
      p_1_position = box%particle(p_1)%position

      ! Choose the second one randomly from any of all the neighboring cells.
      call searchNeighbor(box,reactionCell_1,nNeighbors,speciesOfReagent_2,p_2,p_2_prev,reactionCell_2)
      p_2_position = box%particle(p_2)%position
      if (p_1/=p_2) then
         if (TestOverlapping(p_1_position,p_2_position, &
                  speciesDiameter(speciesOfReagent_1),speciesDiameter(speciesOfReagent_2))) then
            if (box%nextParticle(p_1)==p_2) then
               call deleteParticle_IRDME(p_2,p_2_prev,speciesOfReagent_2,reactionCell_2)
               call deleteParticle_IRDME(p_1,p_1_prev,speciesOfReagent_1,reactionCell_1)
            else
               call deleteParticle_IRDME(p_1,p_1_prev,speciesOfReagent_1,reactionCell_1)
               call deleteParticle_IRDME(p_2,p_2_prev,speciesOfReagent_2,reactionCell_2)
            end if
            ! Decide at whose position the new particle is created
            call UniformRNG(random)
            if (random<0.5_wp) then
               p_1 = box%freeParticle
               call createParticle_IRDME(p_1,speciesOfProduct_1,reactionCell_1)
            else
               p_2 = box%freeParticle
               call createParticle_IRDME(p_2,speciesOfProduct_1,reactionCell_2)
            end if
            changeInCell_1 = .true.
            changeInCell_2 = .true.
            box%reactionCount(choiceOfReaction)=box%reactionCount(choiceOfReaction)+1
            box%reactionCount(0)=box%reactionCount(0)+1
         end if
      end if


      ! 6, Cata_Annihilation: (1) A + B --> A (2) A + A --> A 
      ! In this case B disappears and A remains unchanged
      case(R_CATA_ANNIHILATION)

      ! Choosee a random reagent particle.
      call UniformInteger(candidate,box%numberPerCellList(speciesOfReagent_1,reactionCell_1))
      p_1 = box%firstParticleList(speciesOfReagent_1,reactionCell_1)
      call findCandidateParticle(box,candidate,p_1,p_1_prev)
      p_1_position = box%particle(p_1)%position

      ! Choose the second one randomly from any of all the neighboring cells.
      call searchNeighbor(box,reactionCell_1,nNeighbors,speciesOfReagent_2,p_2,p_2_prev,reactionCell_2)
      p_2_position = box%particle(p_2)%position
      if (p_1/=p_2) then
         if (TestOverlapping(p_1_position,p_2_position, &
                  speciesDiameter(speciesOfReagent_1),speciesDiameter(speciesOfReagent_2))) then

            if (speciesOfReagent_1==speciesOfReagent_2) then
               call UniformRNG(random)
               if (random<0.5_wp) then
                  call deleteParticle_IRDME(p_1,p_1_prev,speciesOfReagent_1,reactionCell_1)
               else
                  call deleteParticle_IRDME(p_2,p_2_prev,speciesOfReagent_2,reactionCell_2)
               end if
            else if (speciesOfReagent_1==speciesOfProduct_1) then
               call deleteParticle_IRDME(p_2,p_2_prev,speciesOfReagent_2,reactionCell_2)
            else
               call deleteParticle_IRDME(p_1,p_1_prev,speciesOfReagent_1,reactionCell_1)
            end if

            changeInCell_1 = .true.
            changeInCell_2 = .true.
            box%reactionCount(choiceOfReaction)=box%reactionCount(choiceOfReaction)+1
            box%reactionCount(0)=box%reactionCount(0)+1
         end if
      end if


      ! 7, Transform:    (1) A + A --> B + C (2) A + B --> C + D (3) A + A --> A + B (4) A + B --> A + A 
      case(R_TRANSFORM)
      ! No new particle positions are generated here, they are simply swapped around

      ! Choosee a random reagent particle.
      call UniformInteger(candidate,box%numberPerCellList(speciesOfReagent_1,reactionCell_1))
      p_1 = box%firstParticleList(speciesOfReagent_1,reactionCell_1)
      call findCandidateParticle(box,candidate,p_1,p_1_prev)
      p_1_position = box%particle(p_1)%position

      ! Choose the second one randomly from any of all the neighboring cells.
      call searchNeighbor(box,reactionCell_1,nNeighbors,speciesOfReagent_2,p_2,p_2_prev,reactionCell_2)
      p_2_position = box%particle(p_2)%position
      if (p_1/=p_2) then
         if (TestOverlapping(p_1_position,p_2_position, &
                  speciesDiameter(speciesOfReagent_1),speciesDiameter(speciesOfReagent_2))) then

            call UniformRNG(random)
            if (random<0.5_wp) then
               if (box%nextParticle(p_1)==p_2) then
                  call replaceParticle_IRDME(p_2,p_2_prev,speciesOfReagent_2,speciesOfProduct_2,reactionCell_2)
                  call replaceParticle_IRDME(p_1,p_1_prev,speciesOfReagent_1,speciesOfProduct_1,reactionCell_1)
               else
                  call replaceParticle_IRDME(p_1,p_1_prev,speciesOfReagent_1,speciesOfProduct_1,reactionCell_1)
                  call replaceParticle_IRDME(p_2,p_2_prev,speciesOfReagent_2,speciesOfProduct_2,reactionCell_2)
               end if
            else
               if (box%nextParticle(p_1)==p_2) then
                  call replaceParticle_IRDME(p_2,p_2_prev,speciesOfReagent_2,speciesOfProduct_1,reactionCell_2)
                  call replaceParticle_IRDME(p_1,p_1_prev,speciesOfReagent_1,speciesOfProduct_2,reactionCell_1)
               else
                  call replaceParticle_IRDME(p_1,p_1_prev,speciesOfReagent_1,speciesOfProduct_2,reactionCell_1)
                  call replaceParticle_IRDME(p_2,p_2_prev,speciesOfReagent_2,speciesOfProduct_1,reactionCell_2)
               end if
            end if

            changeInCell_1 = .true.
            changeInCell_2 = .true.
            box%reactionCount(choiceOfReaction)=box%reactionCount(choiceOfReaction)+1
            box%reactionCount(0)=box%reactionCount(0)+1
         end if
      end if


      ! 8, Catalyst:     reaction of form A + B --> A + C
      case(R_CATALYST)
      ! In catalystic reaction A + B --> A + C, the position of catalyst A remains unchanged, so we only need to find out what B and C are

      ! Choosee a random reagent particle.
      call UniformInteger(candidate,box%numberPerCellList(speciesOfReagent_1,reactionCell_1))
      p_1 = box%firstParticleList(speciesOfReagent_1,reactionCell_1)
      call findCandidateParticle(box,candidate,p_1,p_1_prev)
      p_1_position = box%particle(p_1)%position

      ! Choose the second one randomly from any of all the neighboring cells.
      call searchNeighbor(box,reactionCell_1,nNeighbors,speciesOfReagent_2,p_2,p_2_prev,reactionCell_2)
      p_2_position = box%particle(p_2)%position
      if (p_1/=p_2) then
         if (TestOverlapping(p_1_position,p_2_position, &
                  speciesDiameter(speciesOfReagent_1),speciesDiameter(speciesOfReagent_2))) then

            if (speciesOfReagent_1==speciesOfProduct_1) then
               call replaceParticle_IRDME(p_2,p_2_prev,speciesOfReagent_2,speciesOfProduct_2,reactionCell_2)
            else if (speciesOfReagent_1==speciesOfProduct_2) then
               call replaceParticle_IRDME(p_2,p_2_prev,speciesOfReagent_2,speciesOfProduct_1,reactionCell_2)
            else if (speciesOfReagent_2==speciesOfProduct_1) then
               call replaceParticle_IRDME(p_1,p_1_prev,speciesOfReagent_1,speciesOfProduct_2,reactionCell_1)
            else
               call replaceParticle_IRDME(p_1,p_1_prev,speciesOfReagent_1,speciesOfProduct_1,reactionCell_1)
            end if

            changeInCell_1 = .true.
            changeInCell_2 = .true.
            box%reactionCount(choiceOfReaction)=box%reactionCount(choiceOfReaction)+1
            box%reactionCount(0)=box%reactionCount(0)+1
         end if
      end if



      case default
         stop 'Unknown reaction type'

      end select

      ! Now we need to update all entries in the event queue that got changed:
      if (reactionCell_1==reactionCell_2) then
         reactionCell_2 = 0
         changeInCell_1 = changeInCell_1 .or. changeInCell_2
         changeInCell_2 = .false.
      end if
      call updateQueue(box,reactionCell_1,reactionCell_2,changeInCell_1,changeInCell_2)
      if (isAssertsOn) call LLCtest(box) 

   end do SSALoop

contains

   subroutine calculatePropensity_IRDME(box,cellIndex,propensity,nNeighbors)
      type (DoiBox), target :: box
      integer, intent(in) :: cellIndex
      real(wp), intent(out) :: propensity(0:totalReactions)
      integer, dimension(nSpecies), intent(out) :: nNeighbors

      integer :: iReaction, pairs, neighborParticles, reagent_1, reagent_2, i,j,k, i_2, j_2, k_2, reagentCount, iCell
      integer, dimension(nMaxDimensions) :: tempCell

      do iSpecies=1,nSpecies
         nNeighbors(iSpecies) = 0
      end do
      do reagent = 1, nSpecies
         do iCell = 1, neighborhoodSize(0)
            nNeighbors(reagent) = nNeighbors(reagent)+box%numberPerCellList(reagent,box%neighboringCells(iCell,cellIndex))
         end do
      end do

      propensity = 0.0_wp
      do iReaction = 1, totalReactions
         reagent_1 = reactionTable(iReaction,1)
         reagent_2 = reactionTable(iReaction,2)
         reagentCount = 0
         if (reagent_1 /= 0) reagentCount = reagentCount +1
         if (reagent_2 /= 0) reagentCount = reagentCount +1

         select case (reagentCount)

         case (0)
            propensity(iReaction)=IRDME_reactionRate(iReaction)*DoiCellVolume
         case (1)
            propensity(iReaction)=box%numberPerCellList(reagent_1,cellIndex)*IRDME_reactionRate(iReaction)
         case (2)
            pairs = box%numberPerCellList(reagent_1,cellIndex) * nNeighbors(reagent_2)
            ! Assuming four cases:
            ! For A+B->?, which is duplicated as A+B-> or B+A->
            ! 1) A in center cell, B in center cell; Scheduled either as A+B or as B+A pair twice
            ! 2) A in center cell, B in neighbor center (not center); Scheduled for center cell and for neighbor cell so twice
            ! Or for A+A->?
            ! 3) p1 A in center cell, p2 A in center cell; Scheduled twice as p1-p2 or p2-p1 reaction
            ! 4) A in center cell, A in neighbor cell (not center); Scheduled once for center cell and once for neighbor cell, so twice
            ! In each case, the reaction is scheduled exactly twice, so we need to multiply by 1/2 for binary reactions
            propensity(iReaction)=IRDME_reactionRate(iReaction)*pairs*0.5_wp
         case default
            stop 'Wrong reagent count.'
         end select

         propensity(0)=propensity(0) + propensity(iReaction)
      end do

   end subroutine calculatePropensity_IRDME
   
   
   ! We always need to update the event prediction for reactionCell_1, but sometimes we also need to update many more neighboring cells
   ! If the reaction actually happens and some particle is created or deleted in reactionCell_1,
   ! then we set change_1 to be true and otherwise we set it to be false.
   ! If change_1 is true then we need to update all its neighbors. Otherwise, we just update the cell reactionCell_1.
   ! If it's a unary or birth reaction, we set reactionCell_2 to be 0. In this case, change_2 is also set to false.
   subroutine updateQueue(box,reactionCell_1,reactionCell_2,change_1,change_2)
      type (DoiBox), target :: box
      integer, intent(in) :: reactionCell_1, reactionCell_2
      logical, intent(in) :: change_1, change_2 ! Has the number of particles actually changed?

      real (wp) :: propensity(0:totalReactions)
      integer :: iCell,listIndex,iDimension,newCount
      integer, dimension(2*(3**nDimensions)) :: updateList
      integer, dimension(nMaxDimensions) :: tempCell, cellLocation
      logical :: test, newCell

      updateList = 0 ! List of all cells that need to be updated in the queue
      listIndex = 0
      if (change_1) then ! Add all neighbors (this includes reactionCell_1 itself)
         do iCell = 1, neighborhoodSize(0)
            listIndex = listIndex + 1
            updateList(listIndex)=box%neighboringCells(iCell,reactionCell_1)
         end do
      else if (.not. change_2) then ! Add just reactionCell_1 itself
         ! QY: This *if* check is used to avoid double counting reactinoCell_1.
         ! Without the check, here we include reactionCell_1 and then later we include it again if change_2==.true.
         listIndex = listIndex + 1
         updateList(listIndex) = reactionCell_1
      end if
      
      if (reactionCell_2/=0) then
         if (change_2) then
            if (change_1) then
            
               newCount = 0
               call revertCellIndex(box,reactionCell_1,cellLocation)
               do iCell = 1, neighborhoodSize(0)
                  call revertCellIndex(box,box%neighboringCells(iCell,reactionCell_2),tempCell)
                  ! QY: The code below is equivalent to minimal image convention.
                  newCell = .false.
                  do iDimension = 1, nDimensions
                     if ((abs(cellLocation(iDimension)-tempCell(iDimension))>1).and.&
                        (abs(cellLocation(iDimension)-tempCell(iDimension))/=box%extent(iDimension)-1)) newCell = .true.
                  end do
                  if (newCell) then
                     listIndex = listIndex + 1
                     updateList(listIndex) = box%neighboringCells(iCell,reactionCell_2)
                     newCount = newCount + 1
                  end if
               end do
               
               if (isAssertsOn) then
                  if (newCount>19) then ! QY: Just to be safe
                     write (*,*) 'newCount is', newCount
                     call revertCellIndex(box,reactionCell_1,tempCell)
                     write (*,*) 'ReactionCell_1 is', tempCell, 'its neighbors are'
                     do iCell = 1, neighborhoodSize(0)
                        call revertCellIndex(box,box%neighboringCells(iCell,reactionCell_1),tempCell)
                        write (*,*) tempCell
                     end do
                     call revertCellIndex(box,reactionCell_2,tempCell)
                     write (*,*) 'ReactionCell_2 is', tempCell, 'its neighbors are'

                     do iCell = 1, neighborhoodSize(0)
                        call revertCellIndex(box,box%neighboringCells(iCell,reactionCell_2),tempCell)
                        write (*,*) tempCell
                     end do

                     stop 'How can there be so many new counts?'
                  end if
               end if
               
            else
               do iCell = 1, neighborhoodSize(0)
                  listIndex = listIndex + 1
                  updateList(listIndex) = box%neighboringCells(iCell,reactionCell_2)
               end do
            end if
         end if
      end if
      
      ! Now actually update the event predictions for all the tagged cells
      do iCell = 1,listIndex
         cellIndex = updateList(iCell)
      
         call calculatePropensity_IRDME(box,cellIndex,propensity,nNeighbors)
         if (propensity(0)<=0.0_wp) then
            call deleteFromHeap(box%heap,cellIndex) ! Safe to call even if deleted already
         else
            call calculateTau(propensity(0),tau)
            ! Insert the newly computed time to the heap
            tau = box%elapsedTime+tau
            if(tau > timestep) then ! No need to insert this into the heap at all
               call deleteFromHeap(box%heap,cellIndex) ! Safe to call even if deleted already
            else
               call insertInHeap(box%heap,cellIndex,tau) ! This will call updateHeap if already in the heap
            end if   
         end if
      end do
            
      if (isAssertsOn) then
         call testHeap(box%heap)
         call updateListCheck(box,reactionCell_1,reactionCell_2,change_1,change_2,updateList,listIndex)
      end if
      
   end subroutine updateQueue

   subroutine deleteParticle_IRDME(p,p_prev,species,cellIndex)
      integer, intent(in) :: p, p_prev, species, cellIndex

      integer :: p_next

      box%nParticles(0)=box%nParticles(0)-1
      box%nParticles(species)=box%nParticles(species)-1
      box%particle(p)%species = 0
      
      ! Remove the particle from the LLCs
      box%numberPerCellList(species,cellIndex)=box%numberPerCellList(species,cellIndex)-1
      box%numberPerCellList(0,cellIndex)=box%numberPerCellList(0,cellIndex)-1
      p_next = box%nextParticle(p)
      if (p_prev /= 0) box%nextParticle(p_prev)=p_next
      !if (p_next /= 0) box%previousParticle(p_next)=p_prev
      if (p == box%firstParticleList(species,cellIndex)) then
         if (box%numberPerCellList(species,cellIndex)>0) then
            box%firstParticleList(species,cellIndex)=p_next
         else
            box%firstParticleList(species,cellIndex)=0
         end if
      end if
      
      ! Add this particle to the list of free particles
      box%nextParticle(p)=box%freeParticle
      box%freeParticle=p
     
   end subroutine deleteParticle_IRDME

   subroutine createParticle_IRDME(p,species,cellIndex)
      integer, intent(in) :: p, species, cellIndex

      integer :: p_temp
      box%nParticles(0)=box%nParticles(0) + 1
      box%nParticles(species)=box%nParticles(species) + 1
      box%particle(p)%species = species
      
      ! Add this particle to the LLCs:
      box%numberPerCellList(species,cellIndex)=box%numberPerCellList(species,cellIndex)+1
      box%numberPerCellList(0,cellIndex)=box%numberPerCellList(0,cellIndex)+1
      p_temp = box%firstParticleList(species,cellIndex)
      box%firstParticleList(species,cellIndex) = p
      !if (p_temp/=0) box%previousParticle(p_temp) = p
      call updateFreeParticle(box)
      box%nextParticle(p) = p_temp
   end subroutine createParticle_IRDME


   subroutine replaceParticle_IRDME(p,p_prev,species_old,species_new,cellIndex)
      integer, intent(in) :: p, p_prev, species_old, species_new, cellIndex

      integer :: p_next, p_temp

      if (species_old==species_new) return
      box%nParticles(species_old)=box%nParticles(species_old)-1
      box%nParticles(species_new)=box%nParticles(species_new) + 1
      box%particle(p)%species = species_new
      
      ! Remove the particle from the LLCs
      box%numberPerCellList(species_old,cellIndex)=box%numberPerCellList(species_old,cellIndex)-1
      p_next = box%nextParticle(p)
      if (p_prev /= 0) box%nextParticle(p_prev)=p_next
      !if (p_next /= 0) box%previousParticle(p_next)=p_prev
      if (p == box%firstParticleList(species_old,cellIndex)) then
         if (box%numberPerCellList(species_old,cellIndex)>0) then
            box%firstParticleList(species_old,cellIndex)=p_next
         else
            box%firstParticleList(species_old,cellIndex)=0
         end if
      end if

      box%numberPerCellList(species_new,cellIndex)=box%numberPerCellList(species_new,cellIndex)+1
      p_temp = box%firstParticleList(species_new,cellIndex)
      box%firstParticleList(species_new,cellIndex) = p
      !if (p_temp/=0) box%previousParticle(p_temp) = p
      box%nextParticle(p) = p_temp
   end subroutine replaceParticle_IRDME

   ! For debugging only:
   subroutine updateListCheck(box,reactionCell_1,reactionCell_2,change_1,change_2,updateList,listIndex)
      type (DoiBox), target :: box
      integer, intent(in) :: reactionCell_1, reactionCell_2
      logical, intent(in) :: change_1, change_2 ! Has the number of particles actually changed?
      integer, intent(in), dimension(2*(3**nDimensions)) :: updateList
      integer, intent(in) :: listIndex

      integer :: iCell, jCell, iDimension, count, realCount
      integer, dimension(nMaxDimensions) :: cellLocation_1, cellLocation_2,cellLocation

      count = 0
      do iCell = 1,size(updateList)
         if (updateList(iCell)==0) then
            count = iCell -1
            exit
         end if
      end do
      if (count==0) count = size(updateList)
      if (count/=listIndex) then
         write (*,*) 'Count is', count, 'listIndex is', listIndex
         stop 'Wrong count'
      end if
      do iCell = 1, count
         if (updateList(iCell)==0) stop 'Wrong 0 in update list'
      end do
      do iCell = count+1, size(updateList)
         if (updateList(iCell)/=0) stop 'Wrong non-zero in update list'
      end do
      do iCell = 1,count
         do jCell = iCell+1, count
            if (updateList(iCell)==updateList(jCell)) then
               call revertCellIndex(box,reactionCell_1,cellLocation)
               write (*,*) 'reactionCell_1:', cellLocation
               call revertCellIndex(box,reactionCell_2,cellLocation)
               write (*,*) 'reactionCell_2:', 'index', reactionCell_2, 'location', cellLocation
               call revertCellIndex(box,updateList(iCell),cellLocation)
               write (*,*) 'redundant one:', cellLocation
               stop 'Redundance in update list'
            end if
         end do
      end do
      call revertCellIndex(box,reactionCell_1,cellLocation_1)
      call revertCellIndex(box,reactionCell_2,cellLocation_2)
      dist = 0
      do iDimension = 1, nDimensions
         if (cellLocation_1(iDimension)/=cellLocation_2(iDimension) & 
            .and.(box%extent(iDimension)/=3)) dist = dist +1   
      end do
      realCount = 0
      if ((change_1).and.(change_2)) then
         select case (nDimensions)

         case(1)
         if (dist==0) realCount = 3
         if (dist==1) realCount = 4
         
         case(2)
         if (dist==0) realCount = 9
         if (dist==1) realCount = 12
         if (dist==2) realCount = 14

         case(3)
         if (dist==0) realCount = 27
         if (dist==1) realCount = 36
         if (dist==2) realCount = 42
         if (dist==3) realCount = 46

         case default
         stop 'Problem in update list real count'

         end select
      else if ((change_1).and.(.not. change_2)) then
         select case (nDimensions)
         case (1)
         realCount = 3
         case (2)
         realCount = 9
         case (3)
         realCount = 27
         case default
         stop 'Wrong dimensions'
         end select
      else if ((.not.change_1).and.(change_2)) then
         stop 'Wrong change_1 and change_2'
      else
         realCount = 1
      end if
      if (realCount/=count) then
         write (*,*) 'real count is ', realCount, 'count is ', count
         write (*,*) 'reaction cell 1: ', cellLocation_1
         write (*,*) 'reaction cell 2: ', cellLocation_2
         stop 'Wrong count in update list'
      end if

   end subroutine

end subroutine reaction_IRDME

!--------------------------------
! Utility routines
!--------------------------------

! Called when we create a new particle at location box%freeParticle
! to update the data structure keeping track of the next free particle
subroutine updateFreeParticle(box)
   type (DoiBox), target :: box
   
   integer :: next, temp
   
   if(reactionScheme==RDME) then
      ! Here we just add particles at the end of the array
      box%freeParticle = box%freeParticle + 1
      if (box%freeParticle>ubound(box%particle,1)) then
         write (*,*) 'Allocated space for particle is:', ubound(box%particle,1)
         !temp=-1
         !write(temp,*) "BACKTRACE"
         stop 'Ran out of particle storage, add more'
      end if
   else ! Here freeParticle is the head of a list of free particles
      temp = box%freeParticle
      next = box%nextParticle(temp)
      if (next==0) then
         write (*,*) 'Currently there are ', box%nParticles(0), 'particles'
         write (*,*) 'The current free particle is', box%freeParticle
         stop 'Run out of particle storage in IRDME'
      end if
      box%freeParticle = next
   end if
      
   if(isAssertsOn) then
      if(box%particle(box%freeParticle)%species /= 0) stop 'Free particle has nonzero species'
   end if

end subroutine

subroutine computeCellIndices(box,iParticle,cellIndices)
   type (DoiBox), target :: box
   integer, intent(in) :: iParticle
   integer, dimension(nMaxDimensions), intent(out) :: cellIndices

   ! Donev: UNFINISHED. This may not work for reservoir cells since they have negative indices
   ! but maybe this is only called for particles inside the box
   ! Should we maybe use ceiling here?
   cellIndices(1:nDimensions) = min(box%extent(1:nDimensions), &
               1 + int(box%particle(iParticle)%position(1:nDimensions) / DoiCellLength(1:nDimensions)) )
   cellIndices(nDimensions+1:nMaxDimensions)=1 ! This is faster since it avoids doing computations in the trivial dimensions

end subroutine


subroutine getNewPosition(box,n,p,length)
   type (DoiBox), target :: box
   real (wp), intent(in):: length ! Diameter of sphere
   integer, intent(in) :: n,p ! Put particle n within a distance length/2 from particle p

   real (wp), dimension(nDimensions) :: r, temp, reactiveRadius
   integer :: iDimension

   reactiveRadius = 0.5_wp * length
   ! To generate a point uniformly inside a sphere of unit diameter we use rejection from a point inside a unit square
   GeneratingPosition: do
      call UniformRNGVec(r,nDimensions)
      r = (r-0.5_wp)*2.0_wp
      if (sum(r**2)<1) then ! Accept the point
         temp = r * reactiveRadius
         box%particle(n)%position(1:nDimensions) = box%particle(p)%position(1:nDimensions) + temp
         
         ! Correct for periodic BCs:
         do iDimension = 1, nDimensions
            if (box%particle(n)%position(iDimension) < box%boundaryLocation(1,iDimension)) then
               do while (box%particle(n)%position(iDimension) < box%boundaryLocation(1,iDimension))
                  box%particle(n)%position(iDimension) = box%particle(n)%position(iDimension)+domainLength(iDimension)
               end do
            else if (box%particle(n)%position(iDimension) > box%boundaryLocation(2,iDimension)) then
               do while ((box%particle(n)%position(iDimension) > box%boundaryLocation(2,iDimension)))
                  box%particle(n)%position(iDimension) = box%particle(n)%position(iDimension)-domainLength(iDimension)
               end do
            end if
         end do
         exit GeneratingPosition
      end if
   end do GeneratingPosition
end subroutine

! Select a particle uniformly and randomly from one of the neighboring cells
subroutine searchNeighbor(box,reactionCell_1,nNeighbors,species,p_2,p_2_prev,reactionCell_2)
   type (DoiBox), target :: box
   integer, intent(in) :: reactionCell_1,species
   integer, dimension(nSpecies), intent(in):: nNeighbors
   integer, intent(out) :: p_2, p_2_prev, reactionCell_2

   integer :: tempCount, candidateParticleCount, iCell, neighborCount,candidate,p

   reactionCell_2 = 0
   call UniformInteger(candidate,nNeighbors(species))
   tempCount = 0

   do iCell = 1, neighborhoodSize(0)
      neighborCount = box%numberPerCellList(species,box%neighboringCells(iCell,reactionCell_1))
      if ((candidate>tempCount).and.(candidate <= tempCount+neighborCount)) then
         reactionCell_2 = box%neighboringCells(iCell,reactionCell_1)
         candidate = candidate - tempCount
         exit
      end if
      tempCount = tempCount + neighborCount
   end do
   p_2 = box%firstParticleList(species,reactionCell_2)
   call findCandidateParticle(box,candidate,p_2,p_2_prev)
end subroutine searchNeighbor


! Search through the linked list for the candidate-th particle and its previous particle
subroutine findCandidateParticle(box,candidate,p,p_prev)
   type (DoiBox), target :: box
   integer, intent(in) :: candidate
   integer, intent(inout) :: p
   integer, intent(out):: p_prev

   integer :: num

   p_prev = 0
   do num = 1, candidate-2
      p = box%nextParticle(p)
   end do
   if (candidate>1) then
      p_prev = p
      p = box%nextParticle(p)
   end if
end subroutine findCandidateParticle


! This tests whether two particles actually overlap or not
function TestOverlapping(p_1_position,p_2_position,p_1_diameter,p_2_diameter) result(test)
   real (wp), dimension(nDimensions), intent(in) :: p_1_position, p_2_position
   real (wp), intent(in) :: p_1_diameter, p_2_diameter
   logical :: test
   integer :: iDimension
   real (wp) :: dist(nDimensions)

   do iDimension = 1, nDimensions
      dist(iDimension) = abs(p_1_position(iDimension)-p_2_position(iDimension))
      ! Use the minimum-image convention assuming periodic BCs
      if (dist(iDimension)>DomainLength(iDimension)/2) dist(iDimension)=DomainLength(iDimension)-dist(iDimension)
   end do
   test = (sum(dist**2) < (0.5_wp*(p_1_diameter+p_2_diameter))**2)

end function


! Computes the cell index from index in each dimension
! These are written in 3D always but they will work in 1D and 2D also
subroutine computeCellIndex(box,cellLocation,cellIndex)
   type (DoiBox), target :: box
   integer, dimension(nMaxDimensions), intent(in) :: cellLocation
   integer, intent(out) :: cellIndex

   cellIndex = (cellLocation(3)-1)*(box%extent(1))*(box%extent(2)) + (cellLocation(2)-1) * (box%extent(1)) + cellLocation(1)

end subroutine

! This could be optimized further in 1D and 2D
subroutine revertCellIndex(box,cellIndex,cellLocation)
   type (DoiBox), target :: box
   integer, intent(in) :: cellIndex
   integer, dimension(nMaxDimensions), intent(out) :: cellLocation

   integer :: temp(nMaxDimensions), index

   temp(3) = (cellIndex-1)/((box%extent(1))*(box%extent(2)))
   cellLocation(3) = box%lb(3)+temp(3)
   index = (cellIndex-1) - temp(3)*((box%extent(1))*(box%extent(2)))
   temp(2) = index/box%extent(1)
   cellLocation(2) = box%lb(2)+temp(2)
   cellLocation(1) = index + 1 - temp(2)*box%extent(1)

end subroutine

subroutine calculateTau(totalPropensity,tau)
   real (wp), intent (in) :: totalPropensity
   real (wp), intent (out):: tau

   real (wp) :: random

   call UniformRNG(random)
   ! Note that we use 1-random here to make sure the log argument is never 0 (but can be 1)
   tau = -log(1.0_wp-random)/totalPropensity
   
end subroutine calculateTau

! Since the C rotuine UniformInteger is showing some issues with gfortran (probably C binding problem), write our own here
! This samples an integer uniformly distributed in [1,range]
subroutine myUniformInteger(candidate,range)
   integer, intent(out) :: candidate
   integer, intent(in) :: range

   real (wp) :: random

   call UniformRNG(random)
   random = 1-random ! In (0,1]
   ! This reles on the fact that the floating-point RNG returns in [0,1) strictly
   candidate = ceiling(random*range)
   candidate = min(candidate, range)
   
end subroutine

!--------------------------------
! Test routines (for debugging during development only)
!--------------------------------

subroutine assert (condition, string)
   ! Stops program if specified conditions aren't met
   character (len=*), intent (in) :: string
   logical, intent (in) :: condition
   if ( .not. condition) then
      write(*,*) 'Error: assertion failed with this tag: ', string
      stop 'Program terminated by assert'
   end if
end subroutine assert




! Full LLC test
subroutine LLCtest(box)
   type (DoiBox) :: box

   integer :: n, i, iCell, iSpecies, temp,temp_2,i_2,freeCount,free, alterFreeCount


   do iCell = 1, box%extent(0)
      do iSpecies = 1, nSpecies
         n = box%numberPerCellList(iSpecies,iCell)
         temp = box%firstParticleList(iSpecies,iCell)
         do i = 1,n
            if (temp==0) then
               write (*,*) 'Problem in LLC test'
               write (*,*) 'iSpecies is', iSpecies, 'iCell is', iCell
               write (*,*) 'n is', n, 'number is', i
               write (*,*) 'temp is', temp
               temp_2 = box%firstParticleList(iSpecies,iCell)
               do i_2 = 1, n
                  write (*,*) 'number', i_2, 'is', temp_2
                  temp_2 = box%nextParticle(temp_2)
               end do
            end if
            temp = box%nextParticle(temp)
         end do
      end do
   end do
   free = box%freeParticle
   freeCount = 0
   do while (free/=0)
      if (box%particle(free)%species/=0) stop 'Free particle has species different from 0'
      free = box%nextParticle(free)
      freeCount = freeCount + 1
   end do
   alterFreeCount = 0
   do i = lbound(box%particle,1), ubound(box%particle,1)
      if (box%particle(i)%species == 0) then
         alterFreeCount = alterFreeCount + 1
      end if
   end do
   if (alterFreeCount/=FreeCount) then
      write (*,*) 'Alternative Free Count is', alterFreeCount, 'Free Count is', freeCount
      write (*,*) 'reaction count is', box%reactionCount
      stop 'Two counts unequal'
   end if
   if (freeCount+box%nParticles(0)/=ubound(box%particle,1)-lbound(box%particle,1)+1) then
      if ((.not.IRDME_trace).or.(freeCount+box%nParticles(0)/=ubound(box%particle,1)-lbound(box%particle,1))) then
         write (*,*) 'freeCount is', freeCount
         write (*,*) 'nonFreeCount is', box%nParticles(0)
         write (*,*) 'Total count should be', ubound(box%particle,1)-lbound(box%particle,1)+1
         stop 'Mismatch of free and nonfree particles'
      end if
   end if
end subroutine

! The idea of the full tracing test test is that we first execute one real run,
! where around 1000 reactions happen and then we switch to the fake run and do the test.
! What we want to do in some cases is to run the test while the full code (diffusion+reaction) is running
! The only thing that makes this "testing" is that no reactions will actually be processed
subroutine completeTest(box,nSteps)
   type (DoiBox) :: box
   integer, intent(in) :: nSteps

   integer :: listIndex, iParticle, jParticle, iSpecies, jSpecies, iStep, i, listLen
   
   if(IRDME_trace) then ! We only need this if tracing
      allocate (box%timeOfReactions_AA(nMaxTestPairs))
      allocate (box%timeOfReactions_AB(nMaxTestPairs))
      allocate (box%reactionParticle_AA(nMaxTestPairs,2)) ! Allocate list of reacting particles
      allocate (box%reactionParticle_AB(nMaxTestPairs,2))
      allocate (box%count(ubound(box%particle,1)))
      box%count = 0
   end if

   do iStep = 1, nSteps
      call updateDoiBox(box,inputTimestep)
      write (33,"(1000G17.9)") box%IRDME_test_number(3)
      write (34,"(1000G17.9)") box%IRDME_test_number(4)
   end do
   call testIRDME(box,nSteps)


end subroutine


! Test if the box%nextParticle maintains free particles correctly.
subroutine testNextParticle(box)
   type (DoiBox), target :: box

   integer :: free

   free = box%freeParticle
   do while (free/=0)
      if (box%particle(free)%species/=0) then
         write (*,*) 'free particle problem: particle', free, 'has species', box%particle(free)%species
         stop 'wrong'
      end if
      free = box%nextParticle(free)
   end do
end subroutine


subroutine testIRDME(box,nSteps)
   type (DoiBox), target :: box
   integer, intent(in) :: nSteps

   integer :: iParticle,jParticle,iStep,i,j,iReaction,temp,iCell,iSpecies,jSpecies,iHundred
   integer :: listIndex_AA, listIndex_AB, listLen_AA, listLen_AB
   integer :: iDimension
   integer, dimension(nDimensions) :: iCellIndices, jCellIndices
   real (wp) :: meanReactionPerPair, std, testSum, distance
   real(wp), dimension(nMaxDimensions) :: p_1_position, p_2_position, dist
   real (wp) :: V_AB, V_AA
   integer, allocatable :: numberInEachPair_AA (:),numberInEachPair_AB (:)
   real(wp) :: lastReactionTime_AA, lastReactionTime_AB
   integer :: overlappingCount, particleCount
   integer, allocatable :: reactionList_AA (:,:), reactionList_AB (:,:)
   logical :: found

   ! In addition to automated tests also output these to a file so you can plot them and look at the histogram
   ! and compare it to that of the Poisson distribution. This is the safest test one can do
   ! One can add error bars on a histogram: If the number of counts in a given histogram bin is N
   ! the expected standard deviation is sqrt(N) so the bin count is really N +/- sqrt(N)
   ! By putting error bars on the histogram you can see if each pair is really selected like a Poisson process
   do iCell = 1, box%extent(0)
      write(10,*) iCell, box%IRDME_test_number_1(iCell)/box%elapsedTime/sampleCellVolume, &
         2*sqrt(real(box%IRDME_test_number_1(iCell)))/box%elapsedTime/sampleCellVolume
      write(11,*) box%IRDME_test_number_1(iCell)
   end do
   write(*,*) "Empirical average single-particle birth rate per unit time per unit volume is: ", &
      sum(real(box%IRDME_test_number_1)) / box%elapsedTime / domainVolume
   
   particleCount=0
   do iParticle = 1, size(box%IRDME_test_number_2)
      if(box%particle(iParticle)%species==1) then ! Add some output to file here to examine histogram of cells as well 
         particleCount=particleCount+1
         write(12,*) particleCount, box%IRDME_test_number_2(iParticle)/box%elapsedTime, &
            2*sqrt(real(box%IRDME_test_number_2(iParticle)))/box%elapsedTime

      end if
   end do
   write(*,*) "Empirical average particle death rate per unit time is", &
        real(sum(box%IRDME_test_number_2),wp) / (box%nParticles(1)*box%elapsedTime)
   
   overlappingCount = box%nMaxParticles * nMaxParticlesPerCell * (3**nDimensions) ! Upper bound on number of reactive pairs
   allocate (reactionList_AA(overlappingCount,1:2))
   allocate (reactionList_AB(overlappingCount,1:2))

   ! Now test that pairs of particles are selected uniformly
   !---------------
   ! First make a list of all overlapping pairs of particles
   listIndex_AB = 0
   listIndex_AA = 0
   write (*,*) 'Particles count:', box%nParticles
   do iParticle = lbound(box%particle,1),ubound(box%particle,1)
      iSpecies = box%particle(iParticle)%species;
      if ((iSpecies>0).and.(iSpecies/=3)) then
         do jParticle = iParticle+1, ubound(box%particle,1)
            jSpecies = box%particle(jParticle)%species
            if ((jSpecies>0).and.(jSpecies/=3)) then
               if (iSpecies/=jSpecies) then ! Count A-B pairs
                  p_1_position = box%particle(iParticle)%position; p_2_position = box%particle(jParticle)%position
                  if (TestOverlapping(p_1_position,p_2_position,speciesDiameter(iSpecies),speciesDiameter(jSpecies))) then
                     listIndex_AB = listIndex_AB + 1
                     if(listIndex_AB>ubound(reactionList_AB,1)) stop "Insufficient room in reactionList"
                     reactionList_AB(listIndex_AB,1)=iParticle; reactionList_AB(listIndex_AB,2)=jParticle
                  end if
               else if (iSpecies==1) then ! count A-A pairs
                  p_1_position = box%particle(iParticle)%position; p_2_position = box%particle(jParticle)%position
                  if (TestOverlapping(p_1_position,p_2_position,speciesDiameter(iSpecies),speciesDiameter(jSpecies))) then
                     listIndex_AA = listIndex_AA + 1
                     if(listIndex_AA>ubound(reactionList_AA,1)) stop "Insufficient room in reactionList"
                     reactionList_AA(listIndex_AA,1)=iParticle; reactionList_AA(listIndex_AA,2)=jParticle
                  end if
               end if
            end if
         end do
      end if
   end do
   listLen_AA = listIndex_AA
   listLen_AB = listIndex_AB
   ! Now count how many times each of these pairs was selected by using the record in box%reactionParticle
   allocate (numberInEachPair_AA(1:listLen_AA)) ! How many times the pair was selected
   allocate (numberInEachPair_AB(1:listLen_AB)) ! How many times the pair was selected

   numberInEachPair_AB=0
   numberInEachPair_AA=0
   lastReactionTime_AB=0.0_wp
   lastReactionTime_AA=0.0_wp
   write (*,*) 'Number of A-B reaction pair is', listLen_AB
   write (*,*) 'Number of A-A reaction pair is', listLen_AA
   if (reactionScheme==IRDME) then
      do iReaction = 1, box%IRDME_test_number(3)
         iParticle = box%reactionParticle_AB(iReaction,1)
         jParticle = box%reactionParticle_AB(iReaction,2)
         if (jParticle<iParticle) then
            stop 'jParticle<iParticle'
         end if
         found=.false.
         do i = 1, listLen_AB
            if ((reactionList_AB(i,1)==iParticle).and.(reactionList_AB(i,2)==jParticle)) then
               numberInEachPair_AB(i) = numberInEachPair_AB(i) + 1
               found=.true.
               exit
            end if
         end do
         if(.not.found) then
            write(*,*) "ERROR: Reacting pair not found in list ", iParticle, jParticle, &
               " species=", box%particle(iParticle)%species, box%particle(jParticle)%species
            stop "Reacting pair not found in list of overlapping particles"
         end if
         write(23,"(1000G17.9)") box%timeOfReactions_AB(iReaction)-lastReactionTime_AB
         lastReactionTime_AB=box%timeOfReactions_AB(iReaction)
      end do

      do iReaction = 1, box%IRDME_test_number(4)
         iParticle = box%reactionParticle_AA(iReaction,1)
         jParticle = box%reactionParticle_AA(iReaction,2)
         if (jParticle<iParticle) then
            stop 'jParticle<iParticle'
         end if
         found=.false.
         do i = 1, listLen_AA
            if ((reactionList_AA(i,1)==iParticle).and.(reactionList_AA(i,2)==jParticle)) then
               numberInEachPair_AA(i) = numberInEachPair_AA(i) + 1
               found=.true.
               exit
            end if
         end do
         if(.not.found) then
            write(*,*) "ERROR: Reacting pair not found in list ", iParticle, jParticle, &
               " species=", box%particle(iParticle)%species, box%particle(jParticle)%species
            stop "Reacting pair not found in list of overlapping particles"
         end if
         write(24,"(1000G17.9)") box%timeOfReactions_AA(iReaction)-lastReactionTime_AA
         lastReactionTime_AA=box%timeOfReactions_AA(iReaction)
      end do

      do i = 1, listLen_AB
         write(13,"(1000G17.9)") i, numberInEachPair_AB(i)/box%elapsedTime, 2*sqrt(real(numberInEachPair_AB(i)))/box%elapsedTime
         iParticle = reactionList_AB(i,1)
         jParticle = reactionList_AB(i,2)
         if(numberInEachPair_AB(i)==0) then
            write(*,*) "Pair ", reactionList_AB(i,1:2), " never selected, species=", box%particle(reactionList_AB(i,1:2))%species
            ! Write these to a file so one can visualize them (probably best done in 2D)
         end if
      end do

      do i = 1, listLen_AA
         write(14,"(1000G17.9)") i, numberInEachPair_AA(i)/box%elapsedTime, 2*sqrt(real(numberInEachPair_AA(i)))/box%elapsedTime
         iParticle = reactionList_AA(i,1)
         jParticle = reactionList_AA(i,2)
         if(numberInEachPair_AA(i)==0) then
            write(*,*) "Pair ", reactionList_AA(i,1:2), " never selected, species=", box%particle(reactionList_AA(i,1:2))%species
            ! Write these to a file so one can visualize them (probably best done in 2D)
         end if
      end do

   end if
   write(*,*) "Empirical average binary reaction (A+B->0) rate is ", &
        real(sum(numberInEachPair_AB(1:listLen_AB)),wp) / (listLen_AB * box%elapsedTime)
   write(*,*) "Empirical average binary reaction (A+A->0) rate is ", &
        real(sum(numberInEachPair_AA(1:listLen_AA)),wp) / (listLen_AA * box%elapsedTime)

   if (sum(numberInEachPair_AB)/=box%IRDME_test_number(3)) stop 'Wrong count in A+B->0'
   if (sum(numberInEachPair_AA)/=box%IRDME_test_number(4)) stop 'Wrong count in A+A->0'

   write (*,*) '-------------------'
   write (*,*) 'Elapsed time is', box%elapsedTime
   write (*,*) 'Number of particles at test time is', box%nParticles(0), "=", box%nParticles(1:nSpecies)
   write (*,*) 'Number of all reaction is', box%reactionCount(0), "=", box%reactionCount(1:totalReactions)
   write (*,*) 'Number of fake reaction is', box%IRDME_test_number
   write (*,*) 'Number of candidates with equal probability for each test is', &
        box%extent(0), box%nParticles(1), listLen_AB, listLen_AA
   
   if (nDimensions==2) then
      V_AB = pi*(0.5_wp*(speciesDiameter(1)+speciesDiameter(2)))**2
      V_AA = pi*(speciesDiameter(1))**2
   else
      V_AB = 4.0_wp/3.0_wp*pi*(0.5_wp*(speciesDiameter(1)+speciesDiameter(2)))**3
      V_AA = 4.0_wp/3.0_wp*pi*(speciesDiameter(1))**3
   end if
   write (*,*) 'For selection, theoretical means are', reactionRate(1)*domainVolume*inputTimestep*nSteps/box%extent(0), &
   reactionRate(4)*inputTimestep*nSteps, reactionRate(3)*inputTimestep*nSteps, &
   reactionRate(2)*inputTimestep*nSteps
   
end subroutine

!--------------------------------------
! Sampling cell routines, for hybrid or analysis of hydrodynamics
!--------------------------------------

! Obtain the hydrodynamic state in the *interior* of the box only
! Store the result in the sampling arrays: box%numberDensity
subroutine getState (box)
   type (DoiBox), target :: box

   integer :: iCell (nMaxDimensions), i, j, k, iParticle, iSpecies

   do k = box%lbSample(3), box%ubSample(3)
   do j = box%lbSample(2), box%ubSample(2)
   do i = box%lbSample(1), box%ubSample(1)
      iCell = (/i, j, k/)
      if ( any(iCell < 1) .or. any(iCell > nSampleCells(1:nMaxDimensions)) ) cycle
      box%numberDensity(i,j,k,:) = 0
   end do
   end do
   end do

   do iParticle = lbound(box%particle,1), ubound(box%particle,1)
      iSpecies = box%particle(iParticle)%species
      if (iSpecies <= 0) cycle

      ! We must use ceiling here and not int since there are negative positions for reservoir particles!
      ! Donev: UNFINISHED: Make this consistent with computeCellIndices
      iCell(1:nDimensions) = max(box%lbSample(1:nDimensions), &
            ceiling (box%particle(iParticle)%position(1:nDimensions)/sampleCellLength(1:nDimensions)))
      iCell(nDimensions+1:nMaxDimensions)=1
      i = iCell (1)
      j = iCell (2)
      k = iCell (3)
      if ( any(iCell < box%lbSample) .or. any(iCell > box%ubSample(1:nMaxDimensions)) ) cycle

      ! Count this particle now
      box%numberDensity (i, j, k, iSpecies) = box%numberDensity(i, j, k, iSpecies) + 1.0_wp/sampleCellVolume
      box%numberDensity (i, j, k, 0) = box%numberDensity(i, j, k, 0) + 1.0_wp/sampleCellVolume
   end do
   
end subroutine getState

end module DoiBoxModule
