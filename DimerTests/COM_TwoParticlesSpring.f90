program ParticleSpring
    use Precision
    use BoxLibRNGs
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit       ! This is so I can track errors.
    implicit none 

    ! Look at time scales LARGER than tau. 

    integer, parameter                                  :: wp = r_sp, dim = 1
    real(wp), parameter                                 :: pi = 4.0_wp*ATAN(1.0_wp), KB = 1.38065E-23_wp, T = 300.0_wp
    real(wp), dimension(dim)                            :: r1, r2, r_cm, r_rel
    integer                                             :: n, i, nsteps, myunit1, myunit2 
    real(wp)                                            :: tau, dt, k, l0, a, visc, dimensionlessTSSize, mu, D, mu_eff, D_cm, D_d, mu1, mu2
    character(len=128)                                  :: filenameDiff = 'diffusion.txt', filenameRel = 'relDist.txt', nml_file = "diffusiveSpringParam.nml"


    !call read_namelist(nml_file)

    !Below are Parameters which work, to read from namelist.
    integer                                             :: seed = 6  !Started at 5
    n = 100000 !500000
    k = 1.0_wp
    l0 = 0.0_wp  ! / sqrt(kb * T /  k) or 1E-4
    a = 5.29E-6_wp
    visc = 8.9E-4_wp
    dimensionlessTSSize = 0.01
    nsteps = 100     ! Used to be 100

    mu = 1.0_wp !/ (6*pi*visc*a)
    tau = 1.0 / (mu * k)
    D = 1.0_wp ! kb*T*mu ! Uncomment to left, but to calculate diffusive coefficients it is easier to have a nice value
    dt = tau * dimensionlessTSSize   ! Arbitrarily chosen (this is delta t, which is a fraction of tau)

    ! Optional Parameters
    mu1 = 1.0_wp
    mu2 = 1.0_wp

    mu_eff = mu1 + mu2
    D_d = KB * T * mu_eff
    D_cm = KB * T * mu1 * mu2 / (mu1 + mu2)


    ! For simplicity, I overwrite the diffusion coefficients calculated above to make it nicer, but for a real simulation would delete the below two lines
    D_d = 2.0_wp
    D_cm = 0.5_wp

    call SeedRNG(seed) 

    ! Initialize starting positions.
    r1 = -2E-4_wp
    ! Suggest initializing at a distance of l0 to start
    ! An interesting exercise is to figure out how to generate a random vector of length l0...we can discuss
    r2 = 2E-4_wp ! sqrt(l0**2 / dim)       ! Initial distance is l0, 
  
    open(newunit = myunit1, file = filenameDiff)
    open(newunit = myunit2, file = filenameRel) 

    r_cm = mu2 * r1 / mu_eff + mu1 * r2 / mu_eff
    r_rel = r1-r2  
         
    do i = 0, n
        !rcm = 0.5 * (r1 + r2) (for euler uncomment)
        !r_rel = r1-r2   

        write(myunit1,*) r_cm 
        write(myunit2,*) r_rel

        !call Euler_Maruyama(dt, nsteps, mu, k, D, l0, r1, r2)        ! Just one tau step each iterate. So dtau*nsteps = tau
        !call explicitMidpoint(dt, nsteps, mu_eff, k, D_cm, D_d, l0, r_cm, r_rel)
        call implicitTrapezoidal(dt, nsteps, mu_eff, k, D_cm, D_d, l0, r_cm, r_rel)

        
    end do

    close(myunit1)
    close(myunit2)   

    contains

        ! Implicit trapezoidal integrator as detailed in "Multiscale Temporal Integrators for Fluctuating Hydrodynamics" Delong et. al. 
        ! Implements the scheme found in equation (33), that is where L = mu_eff * k * (l0-l12)/l12,
        ! x^{p,n+1} = x^n + dt/2 * L * (x^n + x^{p,n+1}) + sqrt(2 D dt) N_1(0,1)
        ! x^{n+1} = x^n + dt/2 * (L(x^n)x^n + L(x^n+1)x^{n+1}) + sqrt(dt 2D)N_1(0,1)
        ! Note that these are sampled from the same Normal distribution (they will change b/w r_cm and r_d but the predictor and corrector step seem to use the same W increment in the paper)
        subroutine implicitTrapezoidal(dt, nsteps, mu_eff, k, D_cm, D_d, l0, r_cm, r_rel)
            real(wp), intent(in)                        :: dt, mu_eff, k, l0, D_cm, D_d
            integer, intent(in)                         :: nsteps
            real(wp), dimension(dim), intent(inout)     :: r_cm, r_rel

            ! Local Variables
            real(wp), dimension(dim)                    :: disp1, disp2, r_cm_pred, r_rel_pred
            real(wp)                                    :: l12, sdev_cm, sdev_d, L2_n, L_n
            integer                                     :: i


            ! Implicit Midpoint loop
            do i = 1, nsteps
                call NormalRNGVec(numbers = disp1, n_numbers = dim) ! Mean zero and variance one
                call NormalRNGVec(numbers = disp2, n_numbers = dim) ! Mean zero and variance one

                sdev_cm = sqrt(2 * D_cm * dt)
                sdev_d = sqrt(2 * D_d * dt)

                ! Will want to apply one Implicit Midpoint step on BOTH rcm and rd (rcm can be commented out)

                ! COM STEP
                ! Predictor Step
                r_cm_pred = r_cm + sdev_cm * disp1
                !   Corrector Step
                r_cm = r_cm + sdev_cm * disp1

                ! R DIFFERENCE STEP. Note I think that for both predicting and correcting step in the implicit trap method is the SAME Wiener increment, as in the paper there is no subscript. (Not true for explicit midpoint)
                l12 = norm2(r_rel) ! Evaluate l12 when we are at x_n. Note that this is DEPENDANT on where rd is so l12 = l12(rd). Make sure to update appropriately
                L_n = mu_eff * k * (l0 - l12) / l12

                !Predictor Step
                r_rel_pred = r_rel + (dt / 2) * L_n * r_rel + sdev_d * disp2
                r_rel_pred = r_rel_pred / (1 - dt * L_n / 2)

                ! Corrector L has now changed, as l12 evaluated at the predictor stage is now different. 

                r_rel = r_rel + (dt/2) * L_n * r_rel + sdev_d * disp2

                l12 = norm2(r_rel_pred)  ! For evaluation of L_n+1
                L2_n = mu_eff * k * (l0 - l12) / l12

                r_rel = r_rel / (1 - dt * L2_n / 2)
            end do



        end subroutine
    

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
        
        ! Explicit midpoint integrator as detailed in "Multiscale Temporal Integrators for Fluctuating Hydrodynamics" Delong et. al. 
        ! Implements the scheme found in equation (31), that is, where L = mu_eff * k * (l0-l12)/l12,
        ! x^{p,n+1/2} = x^n + dt/2 * L(x^n) * (x^n) + sqrt(D dt) N_1(0,1)
        ! x^{n+1} = x^n + dt * (L(x^{n+1/2})x^{p,n+1/2} + sqrt(dt D) (N_1(0,1) + N_2(0,1))
        subroutine explicitMidpoint(dt, nsteps, mu_eff, k, D_cm, D_d, l0, r_cm, r_rel)
            real(wp), intent(in)                        :: dt, mu_eff, k, l0, D_cm, D_d
            integer, intent(in)                         :: nsteps
            real(wp), dimension(dim), intent(inout)     :: r_cm, r_rel

            ! Local Variables
            real(wp), dimension(dim)                    :: disp1, disp2, r_cm_pred, r_rel_pred, disp3, disp4
            real(wp)                                    :: l12, sdev_cm, sdev_d
            integer                                     :: i


            ! Explicit Midpoint loop
            do i = 1, nsteps
                call NormalRNGVec(numbers = disp1, n_numbers = dim) ! Mean zero and variance one
                call NormalRNGVec(numbers = disp2, n_numbers = dim) ! Mean zero and variance one
                call NormalRNGVec(numbers = disp3, n_numbers = dim) ! Mean zero and variance one
                call NormalRNGVec(numbers = disp4, n_numbers = dim) ! Mean zero and variance one

                l12 = norm2(r_rel)
                sdev_cm = sqrt(2 * D_cm * dt)
                sdev_d = sqrt(2 * D_d * dt)

                ! Will want to apply one Explicit Midpoint on ONLY rd (the rcm code is commented, but I kept it just so I could make sure it was working correctly)

                ! COM STEP
                ! Predictor Step
                r_cm_pred = r_cm + sqrt(1/2.0) * sdev_cm * disp1
                ! Corrector Step
                r_cm = r_cm + sqrt(1/2.0) * sdev_cm * (disp1 + disp2)

                ! R DIFFERENCE STEP
                l12 = norm2(r_rel) ! Evaluate l12 when we are at x_n. Note that this is DEPENDANT on where rd is so l12 = l12(rd). Make sure to update appropriately
                !Predictor Step
                r_rel_pred = r_rel + dt * mu_eff * k * r_rel * (l0 - l12) / (l12 * 2) + sqrt(1 / 2.0) * sdev_d * disp3
                ! Corrector L has now changed, as l12 evaluated at the predictor stage is now different. 
                l12 =norm2(r_rel_pred) ! L evaluated at n + 1/2
                r_rel = r_rel + dt * mu_eff * k * r_rel_pred * (l0 - l12) / (l12) + sqrt(1 / 2.0) * sdev_d * (disp3 + disp4)
                

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