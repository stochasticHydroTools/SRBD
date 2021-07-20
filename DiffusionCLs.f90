module DiffusionCLs
    use Precision
    use MinHeapModule
    use BoxLibRNGs
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit       ! This is so I can track errors.
    implicit none 
    save

    integer, parameter, private                     :: dim = 3  ! Dimension must be consistent with main code
    real(wp), parameter, private                    :: pi = 4.0_wp*ATAN(1.0_wp)
    integer, private :: outputUnit

    real(wp)                                        :: kbT = 4.0E-3_wp ! Boltzmann constant in units of pN * um 
    real(wp)                                        :: a_1=1.0_wp, a_2=1.0_wp, visc=1.0_wp, k_s=1.0_wp, l0=1.0_wp
    real(wp)                                        :: mu_1_0, mu_2_0, tau_s, tau_r, tau_d ! Computed values
    integer                                         :: sde_integrator_enum=4, nsteps_CLs = 1 !sde_integrator_enum = 
    !1 : Euler-Maruyama, 2: Explicit Midpoint, 3: Implicit Trapezoidal, 4: rotation-vibration integrator w/  
    ! 0.5_wp * dt * (G + G_wped), 5: rotationVibration with g(1/2*(x^n+x_wped))
    character(len = 128)                            :: nml_file = "diffCLs.nml", blobInitializer = "FibersDiameterApart.txt", outputFile = "DimerSimulation.txt" 
    logical                                         :: evolve_r_cm = .true.

    
    ! These parameters are in this module since they are required for SRBD, however, they are not actually used in this module
    ! This is done so that minimal changes are made to DoiBoxModule and most changes are here
    integer :: n_dimers = 10 ! How many dimers to start with; if add_springs=.false. this is number of monomers
    integer :: n_fiber_blobs = 10 ! How many particles composing the fibers
    logical :: add_springs=.true.     ! whether to do translational diffusion and add springs
    integer :: nOutputCLsStep=0 ! How often to output positions of CLs to the file outputFile
    logical :: debug_CLs = .false. ! Whether to print out some more info to screen for testing purposes

    private                                         :: eulerMaruyama, implicitTrapezoidal, explicitMidpoint, rotationVibration

    contains 

        ! Reads the namelist and open files and allocate arrays (we may need extra arrays to
        ! keep track of reactions ourselves)
        subroutine initializeCLs(nml_unit)
            integer, intent(in), optional :: nml_unit ! Read namelist from open file if passed in
            
            real(wp)        ::  lA, D_cm
            
            ! Will read namelist and put in all values.
            ! Other things, like filenames to write to are not included, because this is oly needed in
            ! COM_TwoParticlesSpring.f90 and may not be used in all cases when I use this module.            
            mu_1_0 = 1.0_wp / (6 * pi * visc * a_1)
            mu_2_0 = 1.0_wp / (6 * pi * visc * a_2)
            
            D_cm = kbT * mu_1_0 * mu_2_0 / (mu_1_0 + mu_2_0)
            write(*,*) "CLs diffusing with D=", D_cm
            
            tau_s = 1.0_wp / (k_s * (mu_1_0 + mu_2_0))       ! Spring time scale
            !lA = ((l0**4)*k_s*k_s + 6*l0*l0*kbT*k_s+3*kbT*kbT) / (k_s*(k_s*l0*l0 + kbT)) Rotational time scale for non stiff springs (replace l0*l0
            ! with lA)
            tau_r = l0*l0 / (2 * kbT * (mu_1_0 + mu_2_0))  ! Rotational time scale for stiff springs
            tau_d = l0*l0 / D_cm
            
            write(*,*) "CLs: tau_s=", tau_s, " tau_r=", tau_r, " tau_d=", tau_d

            if (nOutputCLsStep>0) open(newunit = outputUnit, file = outputFile)

        end subroutine

        subroutine destroyCLs()
            if (nOutputCLsStep>0) close(outputUnit)
        end subroutine

        ! Donev: This can output a separate file for each step (now a required argument)
        ! Or write to one file but in this case you want to somehow know which step/time a given position refers to
        subroutine outputCLs(step, time, particle, specie, position)
            integer, intent(in)                               :: step
            real(wp), intent(in), optional                    :: time
            ! The rest of these must all be present *together*:
            integer, intent(in), optional                     :: particle, specie
            real(wp), dimension(dim), intent(in), optional    :: position 
            
            if(nOutputCLsStep<=0) return
            
            if(present(particle)) then

               if(specie /= 2) write(outputUnit,*) particle, specie, position
               
            else
            
               if(step<=0) then
                  write(outputUnit,*) ! newline
               else if(present(time)) then
                  write(outputUnit,*) "# step=", step, " t=", time
               else
                  write(outputUnit,*) "# step=", step
               end if     
            
            end if    
            
        end subroutine
        
        ! Subroutine to read in the data from a namelist file. 
        subroutine readCLParameters(nml_unit)
            integer,  intent(in), optional  :: nml_unit
            integer                         :: unit, check

            ! Namelist definition.
            namelist /diffCLs/ k_s, l0, a_1, a_2, visc, sde_integrator_enum, evolve_r_cm, &
                               add_springs, debug_CLs, n_dimers, n_fiber_blobs, nsteps_CLs, blobInitializer, &
                               outputFile, nOutputCLsStep
            
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
            
            read (nml=diffCLs, unit=unit)

            close (unit)
        end subroutine


        subroutine moveDimer(dt, nsteps, mu_1, mu_2, r_1, r_2)
            real(wp), intent(in)                        :: dt, mu_1, mu_2
            integer, intent(in)                         :: nsteps       
            real(wp), dimension(dim), intent(inout)     :: r_1, r_2
            real(wp), dimension(dim) :: r_cm, r_rel

            
            r_cm = mu_2 * r_1 / (mu_1 + mu_2) + mu_1 * r_2 / (mu_1 + mu_2)
            r_rel = r_1 - r_2

            select case(sde_integrator_enum)
            case(1)
               call eulerMaruyama(dt, nsteps, mu_1, mu_2, r_cm, r_rel)
            case(2)
               call explicitMidpoint(dt, nsteps, mu_1, mu_2, r_cm, r_rel)
            case(3)
               call implicitTrapezoidal(dt, nsteps, mu_1, mu_2, r_cm, r_rel) 
            case(4, 5)
                call rotationVibration(dt, nsteps, mu_1, mu_2, r_cm, r_rel)
            end select

            ! Use r_CM and r_rel to reconstruct the displacement of each cross-linker
            r_1 = r_rel + (r_cm * (mu_1 + mu_2) - mu_2 * r_rel) / (mu_1+ mu_2)
            r_2 = (r_cm * (mu_1 + mu_2) - mu_2 * r_rel) / (mu_1+ mu_2)
        
        end subroutine
         
        subroutine eulerMaruyama(dt, nsteps, mu_1, mu_2, r_cm, r_rel)
            real(wp), intent(in)                        :: dt, mu_1, mu_2
            integer, intent(in)                         :: nsteps          
            real(wp), dimension(dim), intent(inout)     :: r_cm, r_rel
            
            ! Local variables           
            real(wp), dimension(dim)                    :: disp1, disp2
            real(wp)                                    :: l12, sdev_cm, sdev_rel, D_rel, D_cm, mu_eff
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
            real(wp), intent(in)                        :: dt, mu_1, mu_2
            integer, intent(in)                         :: nsteps
            real(wp), dimension(dim), intent(inout)     :: r_cm, r_rel

            ! Local Variables
            real(wp), dimension(dim)                    :: disp1, disp2, r_rel_wped
            real(wp)                                    :: l12, sdev_cm, sdev_d, L2_n, L_n, D_cm, D_rel, mu_eff
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
                r_rel_wped = r_rel + (dt / 2) * L_n * r_rel + sdev_d * disp2
                r_rel_wped = r_rel_wped / (1 - dt * L_n / 2)

                ! Corrector L has now changed, as l12 evaluated at the predictor stage is now different. 

                r_rel = r_rel + (dt/2) * L_n * r_rel + sdev_d * disp2

                l12 = norm2(r_rel_wped)  ! For evaluation of L_n+1
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
            real(wp), intent(in)                        :: dt, mu_1, mu_2
            integer, intent(in)                         :: nsteps
            real(wp), dimension(dim), intent(inout)     :: r_cm, r_rel

            ! Local Variables
            real(wp), dimension(dim)                    :: disp1, disp2, disp3, r_rel_wped
            real(wp)                                    :: l12, sdev_cm, sdev_d, D_rel, D_cm, mu_eff
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
                r_rel_wped = r_rel + dt * mu_eff * k_s * r_rel * (l0 - l12) / (l12 * 2) + sqrt(0.5_wp) * sdev_d * disp2
                ! Corrector L has now changed, as l12 evaluated at the predictor stage is now different. 
                l12 =norm2(r_rel_wped) ! L evaluated at n + 1/2
                r_rel = r_rel + dt * mu_eff * k_s * r_rel_wped * (l0 - l12) / (l12) + sqrt(0.5_wp) * sdev_d * (disp2 + disp3)
                

            end do

        end subroutine                

        ! Apply the implicit trapezoidal method to the dl equation and
        ! the Euler-Lie integrator for rotation on the unit sphere, du
        subroutine rotationVibration(dt, nsteps, mu_1, mu_2, r_cm, r_rel)
            real(wp), intent(in)                        :: dt, mu_1, mu_2
            integer, intent(in)                         :: nsteps
            real(wp), dimension(dim), intent(inout)     :: r_cm, r_rel

            ! Local Variables
            real(wp), dimension(dim)                    :: u, disp1, disp3, N_hat
            real(wp)                                    :: L_n, mu_eff, D_rel, G, G_wped, theta, sdev_cm, D_cm
            integer                                     :: i
            real(wp)                                    :: disp2, l, l_wped, explicit

            ! Declare variables
            mu_eff = mu_1 + mu_2 
            D_rel = kbT * mu_eff
            D_cm = kbT * mu_1 * mu_2 / (mu_eff)        ! Diffusive coeff for r_cm

            ! Change coordinates
            l = norm2(r_rel)
            u = r_rel / l

            do i = 1, nsteps

                ! Euler-Lie Integrator for du (must do this first bc uses l)
                ! Rodrigues Rotation Formula
                call NormalRNGVec(numbers = disp1, n_numbers = dim) ! Mean zero and variance one
                ! Rotating angle
                theta = sqrt(2 * D_rel * dt) * norm2(disp1) / l
                
                N_hat = disp1 / norm2(disp1)
                call rotate(u, theta, N_hat)

                ! Implicit Trapezoidal method for dl. 
                L_n = - mu_eff * k_s 
                !G = 2 * kbT * mu_eff / l + l0 * mu_eff * k_s ! Must be scalar so choose 1st component

                ! Random variate for the scalar equation dl
                call NormalRNG(disp2) ! Always a scalar 
                disp2 = sqrt(2 * D_rel * dt) * disp2

                ! Predictor
                l_wped = l + 0.5_wp * dt * L_n * l + dt * explicitTerm(l) + disp2
                l_wped = l_wped / (1 - dt * L_n / 2)

                ! Corrector
                select case(sde_integrator_enum)
                case(4) ! 0.5_wp * dt * (G + G_wped)
                    explicit = 0.5_wp * dt * (explicitTerm(l) + explicitTerm(l_wped)) 
                case(5) ! g(1/2*(x^n+x_wped))
                    explicit = explicitTerm(0.5_wp * (l + l_wped))
                end select

                l = l + 0.5_wp * dt * (L_n * l) + explicit + disp2
                l = l / (1 - dt * L_n / 2)

                ! Euler-Maruyama for r_cm 
                if (evolve_r_cm) then
                    sdev_cm = sqrt(2 * D_cm * dt)
                    call NormalRNGVec(numbers = disp3, n_numbers = dim) ! Mean zero and variance one
                    r_cm = r_cm + sdev_cm * disp3                       
                end if

            end do

            r_rel = u * l
            
            contains 

            function explicitTerm(l)
                real(wp), intent(in)        :: l
                real(wp)                    :: explicitTerm

                explicitTerm = 2 * kbT * mu_eff / l + l0 * mu_eff * k_s

            end function explicitTerm

        end subroutine

        subroutine rotate(u, theta, N)
            real(wp), dimension(dim), intent(inout)     :: u
            real(wp), intent(in)                        :: theta 
            real(wp), dimension(dim), intent(in)        :: N

            ! Assumed that N is a unit vector
            u = u * cos(theta) + cross(N,u)*sin(theta) + N * dot_product(N,u) * (1-cos(theta))

        end subroutine

        function cross(a, b)
            real(wp), dimension(3) :: cross
            real(wp), dimension(3), intent(in) :: a, b
          
            cross(1) = a(2) * b(3) - a(3) * b(2)
            cross(2) = a(3) * b(1) - a(1) * b(3)
            cross(3) = a(1) * b(2) - a(2) * b(1)
        end function cross

        ! Donev: Changed this to a subroutine since it has side effects (calls RNG)
        subroutine randomDimer(randomDisp)
            real(wp), dimension(dim), intent (out) :: randomDisp
            
            real(wp), dimension(dim)                :: orientationVector
            real(wp)                                :: lengthRand

            call NormalRNGVec(numbers = orientationVector, n_numbers = dim) ! Mean zero and variance one
            orientationVector = orientationVector / norm2(orientationVector)  
   
            ! Random length
            call NormalRNG(lengthRand) ! Mean Zero, Variance 1 
            lengthRand = l0 * l0 * sqrt(kbT / k_s) * lengthRand + l0   ! Mean l0, variance of gibbs * l0^2

            randomDisp = orientationVector * lengthRand
   
        end subroutine randomDimer
          
end module
