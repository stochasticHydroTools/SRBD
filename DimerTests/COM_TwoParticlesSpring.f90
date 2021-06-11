program ParticleSpring
    use Precision
    use BoxLibRNGs
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit       ! This is so I can track errors.
    implicit none 

    ! Look at time scales LARGER than tau. 

    integer, parameter                                  :: wp = r_sp, dim = 1
    ! DONEV: For testing, one can set mu=1, D=1, instead of using real values
    ! Alternatively, and even better,
    ! all output should be normalized, e.g., length should be normalized by the standard deviation predicted by Gibbs-Boltzmann distribution
    ! time by tau, etc.
    real(wp), parameter                                 :: pi = 4.0_wp*ATAN(1.0_wp), KB = 1.38065E-23_wp, T = 300.0_wp
    real(wp), dimension(dim)                            :: r1, r2, rcm, r_rel
    integer                                             :: n, i, nsteps, myunit1, myunit2 
    real(wp)                                            :: tau, dt, k, l0, a, visc, dimensionlessTSSize, mu, D
    character(len=128)                                  :: filenameDiff = 'diffusion.txt', filenameRel = 'relDist.txt', nml_file = "diffusiveSpringParam.nml"


    !call read_namelist(nml_file)

    !Below are Parameters which work, to read from namelist.
    integer                                             :: seed = 5
    n = 100000
    k = 0.029_wp
    l0 = 0.0_wp
    a = 5.29E-11_wp
    visc = 8.9E-4_wp
    dimensionlessTSSize = 0.01_wp
    nsteps = 100

    mu = 1.0_wp !/ (6*pi*a*visc)
    tau = 1.0 / (mu * k)
    D = 1.0_wp ! kb*T*mu ! Uncomment to left, but to calculate diffusive coefficients it is easier to have a nice value
    dt = tau * dimensionlessTSSize   ! Arbitrarily chosen (this is delta t, which is a fraction of tau)


    call SeedRNG(seed) 

    ! Initialize starting positions.
    r1 = -2E-4_wp
    ! Suggest initializing at a distance of l0 to start
    ! An interesting exercise is to figure out how to generate a random vector of length l0...we can discuss
    r2 = 2E-4_wp ! sqrt(l0**2 / dim)       ! Initial distance is l0, 

    open(newunit = myunit1, file = filenameDiff)
    open(newunit = myunit2, file = filenameRel) 
         
    do i = 0, n
        rcm = 0.5 * (r1 + r2)
        r_rel = r1-r2   

        write(myunit1,*) rcm 
        write(myunit2,*) r_rel

        call Euler_Maruyama(dt, nsteps, mu, k, D, l0, r1, r2)        ! Just one tau step each iterate. So dtau*nsteps = tau
        
    end do

    close(myunit1)
    close(myunit2)   

    contains

        subroutine Euler_Maruyama(dt, nsteps, mu, k, D, l0, r1, r2)
            real(wp), intent(in)                        :: dt, mu, k, l0, D
            integer, intent(in)                         :: nsteps          
            real(wp), dimension(dim), intent(inout)     :: r1, r2
            
            ! Local variables           
            real(wp), dimension(dim)                    :: disp1, disp2, vel
            real(wp)                                    :: l12, sdev
            integer                                     :: i

            ! Brownian Motion with Deterministic Drift Realization. 
            do i = 1, nsteps
                call NormalRNGVec(numbers = disp1, n_numbers = dim) ! Mean zero and variance one
                call NormalRNGVec(numbers = disp2, n_numbers = dim) ! Mean zero and variance one

                l12 = norm2(r1-r2) 

                sdev = sqrt(2 * D * dt)


                ! If > 1 dimensions:
                if (dim /= 1) then 
                    vel = r1 - r2
                    vel = mu * k * (l12 - l0) * vel / l12

                    r1 = r1 - vel * dt + sdev*disp1 ! Apply one Euler-Maruyama Step to both r1, r2.
                    r2 = r2 + vel * dt + sdev*disp2 

                ! If 1 Dimension
                else 
                    r1 = r1 + mu * k * (r2 - r1 - l0) * dt + sdev*disp1
                    r2 = r2 + mu * k * (r1 - r2 - l0) * dt + sdev*disp2
                end if

            end do

        end subroutine      
                

        ! Subroutine to read in the data from a namelist file. 
        subroutine read_namelist(file_path)
            character(len=*),  intent(in)  :: file_path
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