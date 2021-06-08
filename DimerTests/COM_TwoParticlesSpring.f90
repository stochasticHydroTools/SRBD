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

    integer, parameter                                  :: wp = r_sp, dim = 3
    ! Donev: The values below except for PI should NOT be parameters and should be read from input file
    ! Parameter means "compile-time constant", that is, a value that is known at compile time and never changed
    ! The compiler will typically replace every occurence of these with the actual hardware bit pattern to optimize the most
    real(wp), parameter                                 :: pi = 4.0_wp*ATAN(1.0_wp), KB = 1.38065E-23_wp, T = 300.0_wp
    real(wp), dimension(dim)                            :: r1, r2
    integer                                             :: n, i, nsteps, rc, iounit ! Donev: Even though you can sometimes omit :: please always use it for clarity of code
    real(wp), dimension(:,:), allocatable               :: pos_cm
    real(wp), dimension(:), allocatable                 :: pos_diff
    real(wp)                                            :: tau, dtau, k, l0, a, visc, dimensionlessTSSize, mu, D
    character(len=128)                                  :: filenameDiff = 'diffusion.txt', filenameRel = 'relDist.txt'

    ! Not sure why the below does not work, so I will manually initialize for now. When the below code is run, 
    !namelist / diffusiveSpringParam / n, k, l0, a, visc, dimensionlessTSSize, nsteps
    !inquire (file="diffusiveSpringParam.nml", iostat=rc)
    !open(newunit = iounit, file="diffusiveSpringParam.nml", iostat = rc, action = "read")
    !read(unit = iounit, nml = diffusiveSpringParam, iostat = rc)
    !close(iounit)

    n = 100000
    k = 0.029_wp
    l0 = 1.0E-8_wp
    a = 5.29E-11_wp
    visc = 8.9E-4_wp
    dimensionlessTSSize = 0.001_wp
    nsteps = 1000

    ! Donev: The value 0.01 (dimensionless time step size) should be read from input
    mu = 1.0 / (6*pi*a*visc)
    tau = 1.0 / (mu * k)
    D = KB*T*mu ! 9.0_wp ! Uncomment to left, but to calculate diffusive coefficients it is easier to have a nice value
    dtau = tau * dimensionlessTSSize   ! Arbitrarily chosen


    ! Donev: RNGs require a "seed" to get started. It is a good idea to use a seed so that one can reproduce runs: 
    ! if you give the same seed it will generate the same sequence of "random" numbers
    ! call SeedRNG(5) ! UNCOMMENT later -- seed is an integer read from input file or as command-line argument
    ! In Linux you can use $RANDOM to generate a *truly* random initial seed
    ! So for example run the code like
    ! code.exe $RANDOM < input_file.dat
    ! This way if you run 100s of runs you can make sure that they all use a different seed
    ! Real codes often save an output or log file in addition to input file
    ! and it is common to write the seed there (or to screen) so one can later reproduce the run if there is a bug/problem/etc

    ! Initialize starting positions.
    r1 = 0.0_wp
    ! Suggest initializing at a distance of l0 to start
    ! An interesting exercise is to figure out how to generate a random vector of length l0...we can discuss
    r2 = 1.20*l0       ! Initial distance is l0, 


    ! Allocate Memory for all arrays 
    allocate(pos_cm(dim,0:n))
    allocate(pos_diff(0:n))
        
    do i = 0, n
        pos_cm(:,i) = 0.5 * (r1 + r2)       ! No longer collecting pos1, pos2.
        pos_diff(i) = norm2(r1 - r2)        ! Can be easily modified to collect the norm as well (if desired)

        ! Donev: EM does not need to know about viscosity and a
        ! Instead, it only cares about diffusion coefficient
        ! This allows to separate applied math pieces (solving an ODE) from physics pieces better
        call Euler_Maruyama(dtau, nsteps, mu, k, D, l0, r1, r2, dim)        ! Just one tau step each iterate.
       
    end do

    pos_cm(:,n) = 0.5 * (r1 + r2)       
    pos_diff(n) = norm2(r1 - r2) 

    ! Donev: Why store pos1 and pos2 if they are never used?
    ! Sometimes we do runs that go on for billions of steps and don't want to store huge arrays
    ! So one can instead not store r1 and r2 and instead store r_cm but also r_diff=r1-r2 which are sufficient to reconstruct r1 and r2
    ! Often we don't even store these arrays and instead do write(file,fmt) r_com INSIDE the for loop so we only write to disk but not memory
    ! Since in your case you are analyzing the data after the fact, this is the right thing to do -- you don't need to store any arrays to memory
    ! the process is Markov which means the next point only depends on the current one so we don't need to store a history
    call writeToFile(pos_cm, n, dim, filenameDiff)

    ! Write the relative differences to a file. 
    call writeToFile(pos_diff,n,1,filenameRel)


    deallocate(pos_cm)
    deallocate(pos_diff)
        


    contains

        ! r1 is the position of the bead which we are diffusing by one time step. 
        subroutine Euler_Maruyama(dt, nsteps, mu, k, D, l0, r1, r2, dim)
            ! Donev: Replace KB, a, visc with mu as input (compute mu above in the main program)
            real(wp), intent(in)                        :: dt, mu, k, l0, D
            integer, intent(in)                         :: nsteps, dim            
            real(wp), dimension(dim), intent(inout)     :: r1, r2
            
            ! Local variables           
            real(wp), dimension(dim)                    :: disp1, disp2, vel
            real(wp)                                    :: l12, sdev
            integer                                     :: i


            ! Brownian Motion with Deterministic Drift Realization. 
            ! Donev: The comment below does not seem to be correct to me -- RNGVec should work for vectors of length 1 (any length in fact)
            ! Note that if dim = 1, disp 1 will need to become disp1(1) and same for disp2. 
            do i = 1, nsteps
                call NormalRNGVec(numbers = disp1, n_numbers = dim) ! Mean zero and variance one
                call NormalRNGVec(numbers = disp2, n_numbers = dim) ! Mean zero and variance one

                l12 = norm2(r1-r2) 

                vel = r1 - r2
                vel = mu * k * (l12 - l0) * vel / l12
                sdev = sqrt(2 * D * dt)

                r1 = r1 - vel * dt + sdev*disp1 ! Apply one Euler-Maruyama Step to both r1, r2.
                r2 = r2 + vel * dt + sdev*disp2

            end do

        end subroutine      
                


        ! Write to File
        ! Donev: See comment above about often writing files INSIDE the DO loop
        subroutine writeToFile(pos,n,dim, filename)
            integer, intent(in)                             :: n, dim
            real(wp), dimension(dim,0:n), intent(in)        :: pos
            character(len = *)                              :: filename 
            integer                                         :: j, myunit

            open(newunit = myunit, file = filename) ! Donev: Hard-wiring unit numbers like 12 is generally a bad idea
            ! In a recent revision of Fortran (supported by gfortran) we added an intrinsic routine called new_unit or something like that
            ! google it and if not successful ask me
            do j = 0, n
                write(myunit, *) pos(:,j) ! Donev: You can just do pos(:,j) here to accomplish the same
            end do

            close(myunit)
        

        end subroutine
    
end program