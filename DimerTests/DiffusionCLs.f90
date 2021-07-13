module DiffusionCLs
    use Precision
    use BoxLibRNGs
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit       ! This is so I can track errors.
    implicit none 
    save

    integer, parameter, private                     :: pr = r_sp, dim = 3       ! MUST BE CONSISTENT WITH MAIN CODE!!!!
    real(pr)                                        :: kbT = 4.0E-3_pr ! Boltzmann constant in units of pN * um 
    real(pr)                                        :: a_1=1.0_pr, a_2=1.0_pr, visc=1.0_pr, k_s=1.0_pr, l0=1.0_pr
    real(pr)                                        :: mu_1_0, mu_2_0, tau_s, tau_r ! Computed values
    integer                                         :: sde_integrator_enum=3, nsteps_CLs = 1 !sde_integrator_enum = 
    !1 : Euler-Maruyama, 2: Explicit Midpoint, 3: Implicit Trapezoidal
    character(len = 128)                            :: nml_file = "diffCLs.nml" ! Can be changed to read namelist from different file
    logical                                         :: evolve_r_cm = .true.
    real(pr), parameter, private                    :: pi = 4.0_pr*ATAN(1.0_pr) !1.0_pr/6.0_pr ! Donev: TEMPORARY, set 6*pi=1, 4.0_pr*ATAN(1.0_pr)
    
    ! These parameters are in this module since they are required for SRBD, however, they are not actually used in this module
    ! This is done so that minimal changes are made to DoiBoxModule and most changes are here
    integer :: n_dimers = 10 ! How many dimers to start with; if add_springs=.false. this is number of monomers
    integer :: n_fiber_blobs = 10 ! How many particles composing the fibers
    logical :: add_springs=.true.     ! whether to do translational diffusion and add springs

    private                                         :: eulerMaruyama, implicitTrapezoidal, explicitMidpoint

    contains 

        ! Reads the namelist and open files and allocate arrays (we may need extra arrays to
        ! keep track of reactions ourselves)
        subroutine initializeCLs(nml_unit)
            integer, intent(in), optional :: nml_unit ! Read namelist from open file if passed in

            ! Will read namelist and put in all values.
            ! Other things, like filenames to write to are not included, because this is oly needed in
            ! COM_TwoParticlesSpring.f90 and may not be used in all cases when I use this module.            
            mu_1_0 = 1.0_pr / (6 * pi * visc * a_1)
            mu_2_0 = 1.0_pr / (6 * pi * visc * a_2)
            
            tau_s = 1.0_pr / (k_s * (mu_1_0 + mu_2_0))       ! Spring time scale
            tau_r = l0 * l0 / (2 * kbT * (mu_1_0 + mu_2_0))  ! Rotational time scale

        end subroutine

        ! Assuming 3D. 
        subroutine outputCLs(particle, specie, x, y, z)
            real(pr), intent(in)                    :: x, y, z 
            integer, intent(in)                     :: particle, specie

            write(77,*) particle, specie, x, y, z 
            
        end subroutine
        
        ! Subroutine to read in the data from a namelist file. 
        subroutine readCLParameters(nml_unit)
            integer,  intent(in), optional  :: nml_unit
            integer                         :: unit, check

            ! Namelist definition.
            namelist /diffCLs/ k_s, l0, a_1, a_2, visc, sde_integrator_enum, evolve_r_cm, &
                               add_springs, n_dimers, n_fiber_blobs, nsteps_CLs
            
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
            case(4)
                call seperateTimeScales(dt, nsteps, mu_1, mu_2, r_cm, r_rel)
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

        ! Apply the implicit trapezoidal method to the dl equation and
        ! the Euler-Lie integrator for rotation on the unit sphere, du
        ! Probably best to think of a better subroutine name
        subroutine seperateTimeScales(dt, nsteps, mu_1, mu_2, r_cm, r_rel)
            real(pr), intent(in)                        :: dt, mu_1, mu_2
            integer, intent(in)                         :: nsteps
            real(pr), dimension(dim), intent(inout)     :: r_cm, r_rel

            ! Local Variables
            real(pr), dimension(dim)                    :: u, disp1, disp3, N_hat
            real(pr)                                    :: L_n, mu_eff, D_rel, G, G_pred, theta, sdev_cm, D_cm
            integer                                     :: i
            real(pr), dimension(1)                      :: disp2, l, l_pred

            ! Declare variables
            mu_eff = mu_1 + mu_2 
            D_rel = kbT * mu_eff
            D_cm = kbT * mu_1 * mu_2 / (mu_eff)        ! Diffusive coeff for r_cm

            ! Change coordinates
            l = norm2(r_rel)
            u = r_rel / l(1)

            do i = 1, nsteps

                ! Euler-Lie Integrator for du (must do this first bc uses l)
                ! Rodrigues Rotation Formula
                call NormalRNGVec(numbers = disp1, n_numbers = dim) ! Mean zero and variance one
                ! Rotating angle
                theta = sqrt(2 * D_rel * dt) * norm2(disp1) / l(1) 
                
                N_hat = disp1 / norm2(disp1)
                call rotate(u, theta, N_hat)

                ! Implicit Trapezoidal method for dl. 
                L_n = - mu_eff * k_s 
                G = 2 * kbT * mu_eff / l(1) + l0 * mu_eff * k_s ! Must be scalar so choose 1st component

                ! Random variate for the scalar equation dl
                call NormalRNGVec(numbers = disp2, n_numbers = 1) ! Always a scalar
                disp2 = sqrt(2 * D_rel * dt) * disp2

                ! Predictor
                l_pred = l + 0.5_pr * dt * L_n * l + dt * G + disp2
                l_pred = l_pred / (1 - dt * L_n / 2)

                ! Corrector
                G_pred = 2 * kbT * mu_eff / l_pred(1) + l0 * mu_eff * k_s
                l = l + 0.5_pr * dt * (L_n * l) + 0.5_pr * dt * (G + G_pred) + disp2
                l = l / (1 - dt * L_n / 2)

                ! Euler-Maruyama for r_cm 
                if (evolve_r_cm) then
                    sdev_cm = sqrt(2 * D_cm * dt)
                    call NormalRNGVec(numbers = disp3, n_numbers = dim) ! Mean zero and variance one
                    r_cm = r_cm + sdev_cm * disp3                       
                end if

            end do

            r_rel = u * l(1) 
            
        end subroutine

        subroutine rotate(u, theta, N)
            real(pr), dimension(dim), intent(inout)     :: u
            real(pr), intent(in)                        :: theta 
            real(pr), dimension(dim), intent(in)        :: N

            ! Assumed that N is a unit vector
            u = u * cos(theta) + cross(N,u)*sin(theta) + N * dot_product(N,u) * (1-cos(theta))

        end subroutine

        ! Cross product is well-defined in the 3D case so built own function assuming 3D. 
        ! No intrinsic fortran routine I could find
        function cross(a, b)
            real(pr), dimension(3) :: cross
            real(pr), dimension(3), intent(in) :: a, b
          
            cross(1) = a(2) * b(3) - a(3) * b(2)
            cross(2) = a(3) * b(1) - a(1) * b(3)
            cross(3) = a(1) * b(2) - a(2) * b(1)
        end function cross
          


end module