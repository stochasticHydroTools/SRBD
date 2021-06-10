program ParticleSpring
    use Precision
    use BoxLibRNGs
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit       ! This is so I can track errors.
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
    ! DONEV: For testing, one can set mu=1, D=1, instead of using real values
    ! Alternatively, and even better,
    ! all output should be normalized, e.g., length should be normalized by the standard deviation predicted by Gibbs-Boltzmann distribution
    ! time by tau, etc.
    real(wp), parameter                                 :: pi = 4.0_wp*ATAN(1.0_wp), KB = 1.38065E-23_wp, T = 300.0_wp
    real(wp), dimension(dim)                            :: r1, r2, rcm
    integer                                             :: n, i, nsteps, myunit1, myunit2 ! Donev: Even though you can sometimes omit :: please always use it for clarity of code
    real(wp)                                            :: tau, dtau, k, l0, a, visc, dimensionlessTSSize, mu, D, r_rel
    character(len=128)                                  :: filenameDiff = 'diffusion.txt', filenameRel = 'relDist.txt', nml_file = "diffusiveSpringParam.nml"


    call read_namelist(nml_file, n, k, l0, a, visc, dimensionlessTSSize, nsteps)

    !Below are Parameters which work, to read from namelist.
    ! DONEV: Add integer seed to this and uncomment seeding below
    !n = 100000
    !k = 0.029_wp
    !l0 = 1.0E-8_wp
    !a = 5.29E-11_wp
    !visc = 8.9E-4_wp
    !dimensionlessTSSize = 0.001_wp
    !nsteps = 1000

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
    r2 = sqrt(l0**2 / dim)       ! Initial distance is l0, 

    open(newunit = myunit1, file = filenameDiff)
    open(newunit = myunit2, file = filenameRel) 
         
    do i = 0, n
        rcm = 0.5 * (r1 + r2)
        ! DONEV: Suggest writing to file r1-r2 (a vector) and then Matlab script can compute norms instead
        ! This is more flexible and does not loose information (e.g., the sign of r1-r2 in 1D)
        r_rel = norm2(r1-r2)   !r_rel must be a array of dimension 1 because writeToFile only accepts arrays, will crash if sent a scalar.

        !pos_cm(:,i) = rcm       ! No longer collecting pos1, pos2.
        !pos_diff(i) = r_rel       ! Can be easily modified to collect the displacement as well (if desired)

        write(myunit1,*) rcm(:) ! DONEV: The (:) is unnecessary 
        write(myunit2,*) r_rel

        ! Donev: EM does not need to know about viscosity and a
        ! Instead, it only cares about diffusion coefficient
        ! This allows to separate applied math pieces (solving an ODE) from physics pieces better
        call Euler_Maruyama(dtau, nsteps, mu, k, D, l0, r1, r2, dim)        ! Just one tau step each iterate.
       
    end do

    ! DONEV: It is OK to just skip this for the very last step to avoid duplication
    ! One can also use a do while loop instead of a do loop with a "break"
    rcm = 0.5 * (r1 + r2)
    r_rel = norm2(r1-r2)
    write(myunit1,*) rcm(:)
    write(myunit2,*) r_rel

    close(myunit1)
    close(myunit2)


    !pos_cm(:,n) = 0.5 * (r1 + r2)       
    !pos_diff(n) = norm2(r1 - r2) 
    ! Write the relative differences to a file. 
    !call writeToFile(pos_diff,n,1,filenameRel)


        


    contains

        ! r1 is the position of the bead which we are diffusing by one time step. 
        subroutine Euler_Maruyama(dt, nsteps, mu, k, D, l0, r1, r2, dim)
            ! DONEV: You should not pass dim here as an argument, since this undoes the "parameter" statement above
            ! Fixing dim=3 at compile time is better for execution speed since it means that the compiler knows the size is 3 and can optimize for that
            ! So just remove dim as an argument and don't pass it in (it will be taken from the variable above)
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
                

        ! Subroutine to read in the data from a namelist file. 
        subroutine read_namelist(file_path, n, k, l0, a, visc, dimensionlessTSSize, nsteps)
            !! Reads Namelist from given file.
            ! DONEV: It is not necessary to pass all these parameters. Variables in the parent scoping unit,
            ! which is the program here, are visible here inside the routine. In fact,
            ! you can treat these inlined/contained routines simply as some text that you can cut and paste above where the "call" is
            ! For routines where it is important to clarify what is input/output etc like the ones above it is a good idea to use arguments
            ! but once you have to pass long lists of arguments it is easier to just use scoping
            character(len=*),  intent(in)  :: file_path
            integer,           intent(out) :: n, nsteps
            real,              intent(out) :: k, l0, a, visc, dimensionlessTSSize
            integer                        :: unit, check
    
            ! Namelist definition.
            namelist /diffusiveSpringParam/ n, k, l0, a, visc, dimensionlessTSSize, nsteps
    
            ! Check whether file exists.
            inquire (file=file_path, iostat=check)
    
            ! Here we have some checks, this just makes sure that the file is around. 
            if (check /= 0) then
                write (stderr, '(3a)') 'Error: The file "', trim(file_path), '" does not exist.'
                return
            end if
    
            ! Open and read Namelist file.
            open (action='read', file=file_path, iostat=check, newunit=unit)
            read (nml=diffusiveSpringParam, iostat=check, unit=unit)
    
            ! This is to keep people aware of cases like : End of File runtime errors.
            if (check /= 0) then
                write (stderr, '(a)') 'Error: invalid Namelist format.'
            end if
    
            close (unit)
        end subroutine read_namelist
    
    
end program