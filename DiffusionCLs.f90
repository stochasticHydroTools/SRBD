module DiffusionCLs
    use Precision
    use BoxLibRNGs
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit       ! This is so I can track errors.
    implicit none 
    save

    integer, parameter, private                     :: pr = r_sp, dim = 3       ! MUST BE CONSISTENT WITH MAIN CODE!!!!
    real(pr)                                        :: kbT = 4.0E-3_pr ! Boltzmann constant in units of ???
    real(pr)                                        :: a_1, a_2, visc, k_s, mu_1_0, mu_2_0, tau, l0
    integer                                         :: sde_integrator_enum      != 1: Euler-Maruyama, 2: Explicit Midpoint, 3: Implicit Trapezoidal
    character(len = 128)                            :: nml_file = "diffCLs.nml" ! Can be changed to read namelist from different file
    logical                                         :: evolve_r_cm, add_springs     ! initialized in namelist, must be public, not parameter cause initialized in nml
    real(pr), parameter, private                    :: pi = 4.0_pr*ATAN(1.0_pr)


    private                                         :: eulerMaruyama, implicitTrapezoidal, explicitMidpoint

    contains 

        ! Donev: I made it so that this does not read the namelist, but only initializes, to be consistent with what is done in DoiBoxModule.f90
        ! Reads the namelist and open files and allocate arrays (we may need extra arrays to
        ! keep track of reactions ourselves)
        subroutine initializeCLs(nml_unit)
            integer, intent(in), optional :: nml_unit ! Read namelist from open file if passed in

            ! Will read namelist and put in all values.
            ! Other things, like filenames to write to are not included, because this is oly needed in
            ! COM_TwoParticlesSpring.f90 and may not be used in all cases when I use this module.            
            mu_1_0 = 1.0_pr / (6 * pi * visc * a_1)
            mu_2_0 = 1.0_pr / (6 * pi * visc * a_2)
            
            tau = 1.0_pr / (k_s * (mu_1_0 + mu_2_0))

        end subroutine
        

        ! Subroutine to read in the data from a namelist file. 
        subroutine readCLParameters(nml_unit)
            integer,  intent(in), optional  :: nml_unit
            integer                         :: unit, check

            ! Namelist definition.
            namelist /diffCLs/ k_s, l0, a_1, a_2, visc, sde_integrator_enum, evolve_r_cm, add_springs
            
            if(present(nml_unit)) then
               unit=nml_unit
            else ! Open the file to read namelist from

               ! Check whether file exists.
               inquire (file=nml_file, iostat=check)

               ! Here we have some checks, this just makes sure that the file is around. 
               if (check /= 0) then
                   write (stderr, '(3a)') 'Error: The file "', trim(nml_file), '" does not exist.'
                   return
               end if

               ! Open and read Namelist file.
               open (action='read', file=nml_file, iostat=check, newunit=unit)
            end if
            
            ! Donev: I am not sure using iostat is the best option -- you want the code to abort with a useful error message
            ! I would delete iostat=check here   
            read (nml=diffCLs, iostat=check, unit=unit)

            ! This is to keep people aware of cases like : End of File runtime errors.
            if (check /= 0) then
                write (stderr, '(a)') 'Error: invalid Namelist format.'
                print *, check
            end if

            close (unit)
        end subroutine


        subroutine moveDimer(dt, nsteps, mu_1, mu_2, r_1, r_2)
            real(pr), intent(in)                        :: dt, mu_1, mu_2
            integer, intent(in)                         :: nsteps          
            real(pr), dimension(dim), intent(inout)     :: r_1, r_2

            real(pr), dimension(dim) :: r_cm, r_rel
            
            r_cm = mu_2 * r_1 / (mu_1 + mu_2) + mu_1 * r_2 / (mu_1 + mu_2)
            r_rel = r_1 - r_2

            select case(sde_integrator_enum)
            case(1)
               call eulerMaruyama(dt, nsteps, mu_1, mu_2, r_cm, r_rel)
            case(2)
               call explicitMidpoint(dt, nsteps, mu_1, mu_2, r_cm, r_rel)
            case(3)
               call implicitTrapezoidal(dt, nsteps, mu_1, mu_2, r_cm, r_rel) 
            end select

            ! Use r_CM and r_rel to reconstruct the displacement of each cross-linker
            r_1 = r_rel + (r_cm * (mu_1 + mu_2) - mu_2 * r_rel) / (mu_1+ mu_2)
            r_2 = (r_cm * (mu_1 + mu_2) - mu_2 * r_rel) / (mu_1+ mu_2)
        
        end subroutine
         
        subroutine eulerMaruyama(dt, nsteps, mu_1, mu_2, r_cm, r_rel)
            real(pr), intent(in)                        :: dt, mu_1, mu_2
            integer, intent(in)                         :: nsteps          
            real(pr), dimension(dim), intent(inout)     :: r_cm, r_rel
            
            ! Local variables           
            real(pr), dimension(dim)                    :: disp1, disp2
            real(pr)                                    :: l12, sdev_cm, sdev_rel, D_rel, D_cm, mu_eff
            integer                                     :: i

            ! Since we pass mu_1, mu_2, we must now calculate these quantities here.
            mu_eff = mu_1 + mu_2                            ! Effective Mobility for r_rel
            D_rel = kbT * mu_eff                            ! Diffusive coeff for r_rel
            D_cm = kbT * mu_1 * mu_2 / (mu_eff)        ! Diff coeff for r_cm


            ! Brownian Motion with Deterministic Drift Realization. 
            do i = 1, nsteps
                call NormalRNGVec(numbers = disp1, n_numbers = dim) ! Mean zero and variance one
                call NormalRNGVec(numbers = disp2, n_numbers = dim) ! Mean zero and variance one

                l12 = norm2(r_rel) 

                sdev_cm = sqrt(2 * D_cm * dt)
                sdev_rel = sqrt(2 * D_rel * dt)

                if (evolve_r_cm) then
                    r_cm = r_cm + sdev_cm * disp1                 
                end if

                r_rel = r_rel + mu_eff * k_s * dt * r_rel * (l0 - l12) / l12 + sdev_rel * disp2

            end do

        end subroutine   

        ! Implicit trapezoidal integrator as detailed in "Multiscale Temporal Integrators for Fluctuating Hydrodynamics" 
        ! Delong et. al. 
        ! Implements the scheme found in equation (33), that is where L = mu_eff * k_s * (l0-l12)/l12,
        ! x^{p,n+1} = x^n + dt/2 * L * (x^n + x^{p,n+1}) + sqrt(2 D dt) N_1(0,1)
        ! x^{n+1} = x^n + dt/2 * (L(x^n)x^n + L(x^n+1)x^{n+1}) + sqrt(dt 2D)N_1(0,1)
        subroutine implicitTrapezoidal(dt, nsteps, mu_1, mu_2, r_cm, r_rel)
            real(pr), intent(in)                        :: dt, mu_1, mu_2
            integer, intent(in)                         :: nsteps
            real(pr), dimension(dim), intent(inout)     :: r_cm, r_rel

            ! Local Variables
            real(pr), dimension(dim)                    :: disp1, disp2, r_rel_pred
            real(pr)                                    :: l12, sdev_cm, sdev_d, L2_n, L_n, D_cm, D_rel, mu_eff
            integer                                     :: i

            mu_eff = mu_1 + mu_2            ! Effective Mobility for r_rel
            D_rel = kbT * mu_eff
            D_cm = kbT * mu_1 * mu_2 / (mu_eff)

            ! Implicit Trapezoidal loop
            do i = 1, nsteps
                call NormalRNGVec(numbers = disp2, n_numbers = dim) ! Mean zero and variance one

                sdev_cm = sqrt(2 * D_cm * dt)
                sdev_d = sqrt(2 * D_rel * dt)

                ! Will want to apply one Implicit Trapezoidal step

                ! COM STEP using Euler-Maruyama
                if (evolve_r_cm) then
                    call NormalRNGVec(numbers = disp1, n_numbers = dim) ! Mean zero and variance one
                    r_cm = r_cm + sdev_cm * disp1                       ! This is really EM, since this is exact
                end if

                ! Evaluate l12 when we are at x_n. Note that this is DEPENDANT on where rd is so l12 = l12(rd). 
                ! Make sure to update appropriately
                l12 = norm2(r_rel) 
                L_n = mu_eff * k_s * (l0 - l12) / l12

                !Predictor Step
                r_rel_pred = r_rel + (dt / 2) * L_n * r_rel + sdev_d * disp2
                r_rel_pred = r_rel_pred / (1 - dt * L_n / 2)

                ! Corrector L has now changed, as l12 evaluated at the predictor stage is now different. 

                r_rel = r_rel + (dt/2) * L_n * r_rel + sdev_d * disp2

                l12 = norm2(r_rel_pred)  ! For evaluation of L_n+1
                L2_n = mu_eff * k_s * (l0 - l12) / l12

                r_rel = r_rel / (1 - dt * L2_n / 2)


            end do

        end subroutine
    
        ! Explicit midpoint integrator as detailed in "Multiscale Temporal Integrators 
        ! for Fluctuating Hydrodynamics" Delong et. al. 
        ! Implements the scheme found in equation (31), that is, where L = mu_eff * k_s * (l0-l12)/l12,
        ! x^{p,n+1/2} = x^n + dt/2 * L(x^n) * (x^n) + sqrt(D dt) N_1(0,1)
        ! x^{n+1} = x^n + dt * (L(x^{n+1/2})x^{p,n+1/2} + sqrt(dt D) (N_1(0,1) + N_2(0,1))
        subroutine explicitMidpoint(dt, nsteps, mu_1, mu_2, r_cm, r_rel)
            real(pr), intent(in)                        :: dt, mu_1, mu_2
            integer, intent(in)                         :: nsteps
            real(pr), dimension(dim), intent(inout)     :: r_cm, r_rel

            ! Local Variables
            real(pr), dimension(dim)                    :: disp1, disp2, disp3, r_rel_pred
            real(pr)                                    :: l12, sdev_cm, sdev_d, D_rel, D_cm, mu_eff
            integer                                     :: i

            mu_eff = mu_1 + mu_2                            ! Effective Mobility for r_rel
            D_rel = kbT * mu_eff                            ! Diffusive coeff for r_rel
            D_cm = kbT * mu_1 * mu_2 / (mu_eff)        ! Diffusive coeff for r_cm


            ! Explicit Midpoint loop
            do i = 1, nsteps

                l12 = norm2(r_rel)
                sdev_cm = sqrt(2 * D_cm * dt)
                sdev_d = sqrt(2 * D_rel * dt)


                ! Euler-Maruyama for r_cm
                if (evolve_r_cm) then
                    call NormalRNGVec(numbers = disp1, n_numbers = dim) ! Mean zero and variance one
                    r_cm = r_cm + sdev_cm * disp1                       
                end if

                ! R DIFFERENCE STEP
                l12 = norm2(r_rel) ! Evaluate l12 when we are at x_n
                call NormalRNGVec(numbers = disp2, n_numbers = dim) ! Mean zero and variance one
                call NormalRNGVec(numbers = disp3, n_numbers = dim) ! Mean zero and variance one
                !Predictor Step
                r_rel_pred = r_rel + dt * mu_eff * k_s * r_rel * (l0 - l12) / (l12 * 2) + sqrt(0.5_pr) * sdev_d * disp2
                ! Corrector L has now changed, as l12 evaluated at the predictor stage is now different. 
                l12 =norm2(r_rel_pred) ! L evaluated at n + 1/2
                r_rel = r_rel + dt * mu_eff * k_s * r_rel_pred * (l0 - l12) / (l12) + sqrt(0.5_pr) * sdev_d * (disp2 + disp3)
                

            end do

        end subroutine                


end module
