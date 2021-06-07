program ParticleSpring
    use Precision
    use BoxLibRNGs
    implicit none 

    ! Look at time scales LARGER than tau. 
    ! Donev: Information like the one below should change from run to run:
    ! the best thing is to use an input file to enter these values
    ! In Fortran this is best done using namelists, and this is what the SRBD code uses
    ! It will be worth for you to learn how to use namelists
    ! For this case, I will sample every 100th step where dtau is 0.01*tau. So I am sampling each tau. 

    integer, parameter                  :: wp = r_sp, dim = 3, nsteps = 100
    ! Donev: The values below except for PI should NOT be parameters and should be read from input file
    ! Parameter means "compile-time constant", that is, a value that is known at compile time and never changed
    ! The compiler will typically replace every occurence of these with the actual hardware bit pattern to optimize the most
    real(wp), parameter                 :: k = 10.0_wp, l0 = 1.0_wp, a = 0.2_wp, pi = 4.0_wp*ATAN(1.0_wp), visc = 1.0_wp
    real(wp), dimension(dim)            :: r1, r2, rcm
    integer :: n, i ! Donev: Even though you can sometimes omit :: please always use it for clarity of code
    real(wp), allocatable               :: pos1(:,:), pos2(:,:), pos_cm(:,:) ! Donev: Shorter as dimension(:,:) to not repeat same thing 3 times
    real(wp)                            :: tau, dtau
    character(len=15)                       :: filename = 'data3.txt' ! Donev: Wrote as len=15 is not enough for a file name. It doesn't hurt to make bigger like 128

    ! Donev: The value 0.01 (dimensionless time step size) should be read from input
    tau = (6*pi*a*visc) / (k)
    dtau = tau * 0.01_wp   ! Arbitrarily chosen
    n = 10000 ! Also from input file

    ! Donev: RNGs require a "seed" to get started. It is a good idea to use a seed so that one can reproduce runs: 
    ! if you give the same seed it will generate the same sequence of "random" numbers
    ! call SeedRNG(seed) ! UNCOMMENT later -- seed is an integer read from input file or as command-line argument
    ! In Linux you can use $RANDOM to generate a *truly* random initial seed
    ! So for example run the code like
    ! code.exe $RANDOM < input_file.dat
    ! This way if you run 100s of runs you can make sure that they all use a different seed
    ! Real codes often save an output or log file in addition to input file
    ! and it is common to write the seed there (or to screen) so one can later reproduce the run if there is a bug/problem/etc

    ! Initialize starting positions.
    ! Donev: These are missing _wp ;-)
    r1 = 0.0
    ! Suggest initializing at a distance of l0 to start
    ! An interesting exercise is to figure out how to generate a random vector of length l0...we can discuss
    r2 = 1.0
    rcm = 0.5 * (r1 + r2)

    ! Allocate Memory for all vectors ! Donev: these are "arrays" since rank-2
    ! Donev: Fortran allows you to begin indexing at any value you want, including zero
    ! Since the initial configuration is special, sort of time step 0, it makes more sense to me (but this is a matter of choice so you choose)
    ! to allocate these as pos1(dim,0:n) instead of 1:n+1 as you are doing
    allocate(pos1(dim,n + 1))
    allocate(pos2(dim,n + 1))
    allocate(pos_cm(dim,n + 1))

    ! We will continue to update these position arrays and then will print out only the COM one. The others are for tracking data
    pos1(:,1) = r1
    pos2(:,1) = r2
    ! Donev: The variable rcm seems to be useless, see below. Either use it inside the do loop (suggested) or remove it
    pos_cm(:,1) = rcm
        
    do i = 1, n
        ! Apply one diffusive step to both r1 and r2 simultaneously.
        ! Donev: If you move stores like pos1(:, i + 1) = r1 here you can instead make the loop be i=0:n
        ! and then do pos1(:, i) = r1
        ! and delete the lines 56-59 above and avoid code duplication
        ! Donev: EM does not need to know about viscosity and a
        ! Instead, it only cares about diffusion coefficient
        ! This allows to separate applied math pieces (solving an ODE) from physics pieces better
        call Euler_Maruyama(dtau, nsteps, a, visc, k, l0, r1, r2, dim)        ! Just one tau step each iterate.
        pos1(:, i + 1) = r1
        pos2(:, i + 1) = r2
        pos_cm(:, i + 1) = 0.5 * (pos1(:,i+1) + pos2(:, i + 1))

    end do

    ! Donev: Why store pos1 and pos2 if they are never used?
    ! Sometimes we do runs that go on for billions of steps and don't want to store huge arrays
    ! So one can instead not store r1 and r2 and instead store r_cm but also r_diff=r1-r2 which are sufficient to reconstruct r1 and r2
    ! Often we don't even store these arrays and instead do write(file,fmt) r_com INSIDE the for loop so we only write to disk but not memory
    ! Since in your case you are analyzing the data after the fact, this is the right thing to do -- you don't need to store any arrays to memory
    ! the process is Markov which means the next point only depends on the current one so we don't need to store a history
    call writeToFile(pos_cm, n, dim, filename)

    deallocate(pos1)
    deallocate(pos2)
    deallocate(pos_cm)
        


    contains

        ! r1 is the position of the bead which we are diffusing by one time step. 
        subroutine Euler_Maruyama(dt,nsteps,a,visc,k,l0,r1,r2, dim)
            ! Donev: Replace KB, a, visc with mu as input (compute mu above in the main program)
            real(wp), intent(in) :: dt, a, visc, k, l0
            integer, intent(in) :: nsteps, dim            
            real(wp), dimension(dim), intent(inout) :: r1, r2
            
            ! Donev: Local variables start here not below            
            real(wp), dimension(dim) :: disp1, disp2, temp
            real(wp), parameter :: KB = 1.38065E-23_wp, T = 300.0_wp, pi = 4.0_wp*ATAN(1.0_wp)      ! These are typically constants of problem so not passed as parameters.
            real(wp) mu, l12


            ! Local Variables ! Donev: disp1 etc are also local
            integer :: i

            ! Initialization of constants
            ! Donev: Move to main program before the loop -- this is physics and should be separate from the SODE solver (math)
            mu = ( 1.0 / (6*pi*visc*a) ) 

            ! Brownian Motion with Deterministic Drift Realization. 
            ! Donev: The comment below does not seem to be correct to me -- RNGVec should work for vectors of length 1 (any length in fact)
            ! Note that if dim = 1, disp 1 will need to become disp1(1) and same for disp2. 
            do i = 1, nsteps
                call NormalRNGVec(numbers = disp1, n_numbers = dim) ! Mean zero and variance one
                call NormalRNGVec(numbers = disp2, n_numbers = dim) ! Mean zero and variance one

                l12 = norm2(r1-r2) ! Donev: Function norm2 is UNDEFINED (this won't compile)?

                temp = r1 - r2
                ! Donev: Observe that it is better to define
                ! temp = (mu * k * (l12 - l0) * (temp) / l12)
                ! here instead of repeating this twice below
                ! Even though compilers are good enough to recognize the repetition and create a temporary themselves
                ! it would make your code clearer if you do it, and also avoid repeating the same thing twice which will avoid errors
                ! Same comment goes for sqrt(2*KB*T*mu*dt)
                ! Donev: Names like "temp" are cryptic and not a good idea, try to use something more decriptive like vel (this has units of velocity)

                ! Subtle comment: Avoid unnecessary parenthesis
                ! (a*b)*c
                ! forces the compiler to multiply a*b and THEN multiply by c
                ! a*b*c allows the compiler freedom to multiply in any order and thus allows better optimization
                ! In fact, these optimization rules in Fortran are what make it generally easier to optimize for numerical computing than C
                ! The name Fortran means formula translation...and this is what it is designed for
                r1 = r1 + (mu * k * (l12 - l0) * (-temp) / l12) * dt + sqrt(2*KB*T*mu*dt)*disp1 ! Apply one Euler-Maruyama Step   
                r2 = r2 + (mu * k * (l12 - l0) * (temp) / l12) * dt + sqrt(2*KB*T*mu*dt)*disp2

            end do

        end subroutine      
                


        ! Write to File
        ! Donev: See comment above about often writing files INSIDE the DO loop
        subroutine writeToFile(pos,n,dim, filename)
            integer, intent(in)                     :: n, dim
            real(wp), dimension(3,n+1), intent(in)  :: pos
            character(15)                           :: filename ! Donev: see comment about length above -- here you should use character(len=*)
            ! and the length will be "assumed" from the actual argument
            integer j, i

            open(12, file = filename) ! Donev: Hard-wiring unit numbers like 12 is generally a bad idea
            ! In a recent revision of Fortran (supported by gfortran) we added an intrinsic routine called new_unit or something like that
            ! google it and if not successful ask me
            do j = 1, n + 1
                write(12, *) (pos(i,j), i = 1, dim) ! Donev: You can just do pos(:,j) here to accomplish the same
            end do

            close(12)
        

        end subroutine
    
end program
