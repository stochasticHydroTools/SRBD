module DiffusionCLs
    use Precision
    use MinHeapModule
    use BoxLibRNGs
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit       ! This is so I can track errors.
    implicit none 
    save

    integer, parameter, private                     :: dim = 3!, wp = r_sp   ! Dimension must be consistent with main code
    real(wp), parameter, private                    :: pi = 4.0_wp*ATAN(1.0_wp)
    integer, private :: outputUnit

    real(wp)                                        :: kbT = 4.0E-3_wp ! Boltzmann constant in units of pN * um 
    
    real(wp)                                        :: a_1=1.0_wp, a_2=1.0_wp, visc=1.0_wp, k_s=1.0_wp, l0=1.0_wp    
    ! For isotropic diffusion, these determine the mobility/diffusion of the two monomers in each dimer
    ! For anisotropic diffusion, a_1 and mu_1_0=kbT/(6*pi*eta*a_1) determine the translational diffusion coefficient
    ! while a_2 and mu_2_0=kbT/(6*pi*eta*a_2) determine the rotational diffusion coefficient
    ! Note that this refers to the diffusion of the WHOLE dimer as a body, not to individual monomers
    ! Specifically, for anisotropic diffusion
    ! D_cm = 1/2 * kbT * mu_1_0 (translational diffusion coefficient)
    ! D_r = 2 * kbT * mu_2_0 * <l^(-2)> (rotational diffusion coefficient)
    ! This convention ensures that if a_1=a_2 there is no difference between anisotropicDiffusion=T or F
    real(wp)                                        :: mu_1_0, mu_2_0, mu_par, mu_perp, tau_s, tau_r, tau_d ! Computed values

    integer                                         :: sde_integrator_enum=4, nsteps_CLs = 1 !sde_integrator_enum = 
    !1 : Euler-Maruyama, 2: Explicit Midpoint, 3: Implicit Trapezoidal, 4: rotation-vibration integrator w/  
    ! 0.5_wp * dt * (G + G_pred), 5: rotationVibration with g(1/2*(x^n+x_pred))

    character(len = 128)                            :: nml_file = "diffCLs.nml", blobInitializer = "FibersDiameterApart.txt", outputFile = "DimerSimulation.txt" 
    logical                                         :: evolve_r_cm = .true.

    
    ! These parameters are in this module since they are required for SRBD, however, they are not actually used in this module
    ! This is done so that minimal changes are made to DoiBoxModule and most changes are here
    integer :: n_dimers = 10 ! How many dimers to start with; if add_springs=.false. this is number of monomers
    integer :: n_fiber_blobs = 10 ! How many particles composing the fibers
    logical :: add_springs=.true.     ! whether to do translational diffusion and add springs
    integer :: nOutputCLsStep=0 ! How often to output positions of CLs to the file outputFile
    logical :: debug_CLs = .false. ! Whether to print out some more info to screen for testing purposes
    logical :: anisotropicMobility  = .false.   ! Whether or not parallel mobility is seperate from perpendicular mobility

    private                                         :: eulerMaruyama, implicitTrapezoidal, explicitMidpoint, rotationVibration

    contains 

        ! Reads the namelist and open files and allocate arrays (we may need extra arrays to
        ! keep track of reactions ourselves)
        subroutine initializeCLs(nml_unit)
            integer, intent(in), optional :: nml_unit ! Read namelist from open file if passed in
            
            real(wp)        ::  linv, D_cm, dl
            
            ! Will read namelist and put in all values.
            ! Other things, like filenames to write to are not included, because this is oly needed in
            ! COM_TwoParticlesSpring.f90 and may not be used in all cases when I use this module.            
            mu_1_0 = 1.0_wp / (6 * pi * visc * a_1)
            mu_2_0 = 1.0_wp / (6 * pi * visc * a_2)

            ! Gibbs-Boltzmann is independent of mobility:
            dl = sqrt(kbT/k_s)
            linv = sqrt(2.0_wp) * sqrt(pi) * (erf(l0 * sqrt(2.0_wp) / dl / 2.0_wp) + 1.0_wp) / (2.0_wp * exp(-l0 ** 2 / dl ** 2 / 2.0_wp) * &
                   dl * l0 + sqrt(2.0_wp) * sqrt(pi) * (dl ** 2 + l0 ** 2) * (erf(l0 * sqrt(2.0_wp) / dl / 2.0_wp) + 1.0_wp))
            
            
            ! IF we are in an anisotropic case, then our convention is as follows:
            ! mu_1_0 refers to translational mobility of a monomer. 
            ! mu_2_0 refers to the perpendicular/rotational mobility, mu_perp
            if(anisotropicMobility) then
               
               D_cm = 0.5_wp * kbT * mu_1_0  ! Here mu_1_0 refers to mobility of one "monomer" (D_cm=0.5*(2/3*D_perp+1/3*D_par))
               
               mu_par = 3.0_wp * mu_1_0 - 2.0_wp * mu_2_0 !Compute mu_par from mu_1_0=mu_cm and mu_2_0=mu_r
               if (mu_par < 0 ) stop "Cannot have negative parallel mobility"
               mu_perp = mu_2_0 ! Here mu_2_0 refers to rotation = perp
               
               tau_s = 1.0_wp / (k_s * 2 * mu_par)      ! Spring time scale involves only parallel mobility 
               tau_r = 1.0_wp / (2 * kbT * mu_perp  * linv)  ! Rotational time scale involves only perpendicular mobility
            else   

               D_cm = kbT * mu_1_0 * mu_2_0 / (mu_1_0 + mu_2_0) ! Center of mass diffusion coefficient
               tau_s = 1.0_wp / (k_s * (mu_1_0 + mu_2_0))       ! Spring time scale
               tau_r = 1.0_wp / (kbT * (mu_1_0 + mu_2_0) * linv)  ! Rotational time scale

            end if
            write(*,*) "CLs diffusing with D=", D_cm
            tau_d = l0*l0 / D_cm
            
            write(*,*) "CLs: tau_s=", tau_s, " tau_r=", tau_r, " tau_d=", tau_d

            if (nOutputCLsStep>0) open(newunit = outputUnit, file = outputFile)

        end subroutine

        subroutine destroyCLs()
            if (nOutputCLsStep>0) close(outputUnit)
        end subroutine

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
                               outputFile, nOutputCLsStep, anisotropicMobility
            
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


        subroutine moveDimer(dt, nsteps, immobile, r_1, r_2)
            real(wp), intent(in)                        :: dt
            integer, intent(in)                         :: immobile ! <=0 if both mobile, or 1 or 2 to indicate immobile particle
            integer, intent(in)                         :: nsteps       
            real(wp), dimension(dim), intent(inout)     :: r_1, r_2
            
            ! Local variables:
            real(wp), dimension(dim) :: r_cm, r_rel 
            real(wp)                 :: mu_1, mu_2

            
            if(immobile>2) return ! Nothing to move

            if(anisotropicMobility) then
            
               if(sde_integrator_enum<4) stop "Only rotationVibration integrator supported for anisotropic mobility"
               
               ! Here mu_1 and mu_2 are only used for computing r_cm and not to actually update the dimers
               ! These are not actual mobilities just set to 1 (mobile) or 0 (immobile)
               select case(immobile)
               case(0)
                  mu_1=1.0_wp
                  mu_2=1.0_wp
               case(1)   
                  mu_1=0.0_wp
                  mu_2=1.0_wp
               case(2)   
                  mu_1=1.0_wp
                  mu_2=0.0_wp               
               end select
               
            else 
                select case(immobile)
                case(0)
                   mu_1= mu_1_0
                   mu_2= mu_2_0
                case(1)   
                   mu_1=0.0_wp
                   mu_2=mu_2_0
                case(2)   
                   mu_1=mu_1_0
                   mu_2=0.0_wp               
                end select  
            
            end if
            
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
               call rotationVibration(dt, nsteps, immobile, r_cm, r_rel)
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

                if(evolve_r_cm .and. (D_cm>0.0_wp)) then
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
            real(wp), dimension(dim)                    :: disp1, disp2, r_rel_pred
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
                if(evolve_r_cm .and. (D_cm>0.0_wp)) then
                    call NormalRNGVec(numbers = disp1, n_numbers = dim) ! Mean zero and variance one
                    r_cm = r_cm + sdev_cm * disp1                       ! This is really EM, since this is exact
                end if

                ! Evaluate l12 when we are at x_n. Note that this is DEPENDENT on where rd is so l12 = l12(rd). 
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
            real(wp), intent(in)                        :: dt, mu_1, mu_2
            integer, intent(in)                         :: nsteps
            real(wp), dimension(dim), intent(inout)     :: r_cm, r_rel

            ! Local Variables
            real(wp), dimension(dim)                    :: disp1, disp2, disp3, r_rel_pred
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
                if(evolve_r_cm .and. (D_cm>0.0_wp)) then
                    call NormalRNGVec(numbers = disp1, n_numbers = dim) ! Mean zero and variance one
                    r_cm = r_cm + sdev_cm * disp1                       
                end if

                ! R DIFFERENCE STEP
                l12 = norm2(r_rel) ! Evaluate l12 when we are at x_n
                call NormalRNGVec(numbers = disp2, n_numbers = dim) ! Mean zero and variance one
                call NormalRNGVec(numbers = disp3, n_numbers = dim) ! Mean zero and variance one
                !Predictor Step
                r_rel_pred = r_rel + dt * mu_eff * k_s * r_rel * (l0 - l12) / (l12 * 2) + sqrt(0.5_wp) * sdev_d * disp2
                ! Corrector L has now changed, as l12 evaluated at the predictor stage is now different. 
                l12 =norm2(r_rel_pred) ! L evaluated at n + 1/2
                r_rel = r_rel + dt * mu_eff * k_s * r_rel_pred * (l0 - l12) / (l12) + sqrt(0.5_wp) * sdev_d * (disp2 + disp3)
                

            end do

        end subroutine                

        ! Apply the implicit trapezoidal method to the dl equation and
        ! the Euler-Lie integrator for rotation on the unit sphere, du
        subroutine rotationVibration(dt, nsteps, immobile, r_cm, r_rel)
            real(wp), intent(in)                        :: dt
            integer, intent(in)                         :: nsteps
            integer, intent(in)                         :: immobile ! <=0 if both mobile, or 1 or 2 to indicate immobile particle            
            real(wp), dimension(dim), intent(inout)     :: r_cm, r_rel

            ! Local Variables
            integer                                     :: i
            real(wp), dimension(dim)                    :: u, u_n, disp1, disp3, N_hat
            real(wp)                                    :: L_n, D_u, D_l, G, G_pred, theta, sdev_cm, D_cm, mu_1, mu_2
            real(wp)                                    :: disp2, l, l_pred, explicit, mu_u, mu_l, sdev_cm_perp, sdev_cm_par, sdev_l
                                    
            if (.not. anisotropicMobility) then
                ! These cases are if we are isotropic, then we handle immobility as follows
                if(immobile==1) then
                   mu_1=0.0_wp
                else
                   mu_1 = mu_1_0
                end if     
                if(immobile==2) then
                   mu_2=0.0_wp
                else
                   mu_2 = mu_2_0
                end if

                mu_l = mu_1 + mu_2 
                mu_u = mu_l
                D_cm = kbT * mu_1 * mu_2 / mu_l        ! Diffusion coeff for r_cm
                
            else if(immobile<=0) then ! Both monomers are mobile

                mu_u = 2.0_wp * mu_perp
                mu_l = 2.0_wp * mu_par 
                D_cm = 0.5_wp * kbT * ((2.0_wp/3.0_wp)*mu_perp + (1.0_wp/3.0_wp)*mu_par)
                
            else ! We are doing anisotropicMobility but one monomer is immobile so rotational mobility is cut in half when partly immobilized
            
                mu_u = mu_perp
                mu_l = mu_par 
                D_cm = 0.0_wp ! No motion of center of mass since partly immobilized
                
            end if
            D_u = kbT * mu_u ! Diffusion coefficient for u
            D_l = kbT * mu_l ! Diffusion coefficient for l

            ! Change coordinates
            l = norm2(r_rel)
            u = r_rel / l
            u_n = u  ! This is the first u before iterations and might be used in computing r_cm

            do i = 1, nsteps

                ! Euler-Lie Integrator for du (must do this first bc uses l)
                ! Rodrigues Rotation Formula
                call NormalRNGVec(numbers = disp1, n_numbers = dim) ! Mean zero and variance one
                ! Rotating angle
                theta = sqrt(2 * D_u * dt) * norm2(disp1) / l
                
                N_hat = disp1 / norm2(disp1)
                ! Evolve Brownian Motion on the unit sphere
                call rotate(u, theta, N_hat)

                ! Implicit Trapezoidal method for dl. 
                L_n = - mu_l * k_s 

                ! Random variate for the scalar equation dl
                call NormalRNG(disp2) ! Always a scalar 
                sdev_l = sqrt(2 * D_l * dt) * disp2

                ! Predictor
                l_pred = l + 0.5_wp * dt * L_n * l + dt * explicitTerm(l) + sdev_l
                l_pred = l_pred / (1 - dt * L_n / 2)

                ! Corrector
                select case(sde_integrator_enum)
                case(4) ! 0.5_wp * dt * (G + G_pred)
                    explicit = 0.5_wp * dt * (explicitTerm(l) + explicitTerm(l_pred)) 
                case(5) ! g(1/2*(x^n+x_pred))
                    explicit = explicitTerm(0.5_wp * (l + l_pred))
                end select

                l = l + 0.5_wp * dt * (L_n * l) + explicit + sdev_l
                l = l / (1 - dt * L_n / 2)


                ! Euler-Maruyama for r_cm if D_cm>0
                if(evolve_r_cm .and. (D_cm>0.0_wp)) then
                  if(anisotropicMobility) then
                    call NormalRNGVec(numbers = disp3, n_numbers = dim) ! Mean zero and variance one 
                    sdev_cm_perp = sqrt(kbT * mu_perp * dt)    
                    sdev_cm_par = sqrt(kbT * mu_par * dt)
                    r_cm = r_cm + sdev_cm_perp * disp3 + &
                     (sdev_cm_par - sdev_cm_perp) * u_n * dot_product(u_n, disp3)   ! Dot product with the initial u vector  
                  else                
                    sdev_cm = sqrt(2 * D_cm * dt)
                    call NormalRNGVec(numbers = disp3, n_numbers = dim) ! Mean zero and variance one
                    r_cm = r_cm + sdev_cm * disp3                       
                  end if
                end if


            end do

            r_rel = u * l
            
            contains 

            function explicitTerm(l)
                real(wp), intent(in)        :: l
                real(wp)                    :: explicitTerm

                explicitTerm = 2 * D_l / l + l0 * mu_l * k_s

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

        ! Simulates r_d and r_cm for non-isotropic mobilities
        ! This is solving the WRONG SDE and missing a drift term, so should NEVER be used. 
        subroutine anisotropicDimer(dt, nsteps, mu_perp, mu_para, r_cm, r_rel)
            real(wp), intent(in)                        :: dt, mu_perp, mu_para
            integer, intent(in)                         :: nsteps
            real(wp), dimension(dim), intent(inout)     :: r_cm, r_rel

            ! Local Variables
            real(wp), dimension(dim)                    :: disp1, disp2, r_rel_pred, u, u_n, L_n, K
            real(wp)                                    :: l12, mu_eff_para, mu_eff_perp, mu_cm_para, mu_cm_perp, &
                                                            sdev_cm_perp, sdev_cm_par, sdev_rel_para, sdev_rel_perp
            integer                                     :: i
            
            l12 = norm2(r_rel)
            u_n = r_rel / l12
            u = u_n

            ! Non-isotropic components for mu, assuming equal mobilities
            mu_eff_para = 2 * mu_para
            mu_eff_perp = 2 * mu_perp
            mu_cm_para = 0.5_wp * mu_para
            mu_cm_perp = 0.5_wp * mu_perp


            do i = 1, nsteps

                ! Euler-Maruyama for r_cm
                if (evolve_r_cm) then
                    call NormalRNGVec(numbers = disp1, n_numbers = dim) ! Mean zero and variance one
                    sdev_cm_perp = sqrt(2 * kbT * mu_cm_perp * dt)
                    sdev_cm_par = sqrt(2 * kbT * mu_cm_para * dt)
                    r_cm = r_cm + sdev_cm_perp * disp1 + (sdev_cm_par - sdev_cm_perp) * u_n * dot_product(u_n, disp1)   ! Dot product with the initial u vector                       
                end if

                ! Implicit Trapezoidal for r_rel
                call NormalRNGVec(numbers = disp2, n_numbers = dim) ! Mean zero and variance one  
                sdev_rel_perp = sqrt(2 * kbT * mu_eff_perp * dt)
                sdev_rel_para = sqrt(2 * kbT * mu_eff_para * dt)

                ! Evaluate l12 when we are at x_n. Note that this is DEPENDENT on where rd is so l12 = l12(rd). 
                ! Make sure to update appropriately
                l12 = norm2(r_rel) 
                u = r_rel / l12
                L_n =  (k_s * (l0 - l12) / l12) * (mu_eff_perp * r_rel + &
                        (mu_eff_para - mu_eff_perp) * u * dot_product(u, r_rel) ) ! Technically, this is L*r_rel not just L, so I avoid matricies
                K = sdev_rel_perp * disp2 + (sdev_rel_para - sdev_rel_perp) * u * dot_product(u, disp2)             ! Technically K*W, not just K


                ! Predictor Step
                r_rel_pred = r_rel + (0.5_wp * dt) * L_n + K
                ! Inversion Predictor step.
                r_rel_pred =  alpha(r_rel)*r_rel + beta(r_rel)*u*dot_product(u,r_rel)

                ! Corrector L has now changed, as l12 evaluated at the predictor stage is now different. 
                r_rel = r_rel + (0.5_wp * dt) * L_n + K 
                ! Inversion Corrector Step
                r_rel = alpha(r_rel_pred)*r_rel + beta(r_rel_pred)*u*dot_product(u, r_rel)

            end do



            contains 

            ! First component alpha of matrix inverse of form alpha *I + beta * u u^T
            function alpha(r_rel)
                real(wp), dimension(dim), intent(in)     :: r_rel
                real(wp)                                 :: alpha, l12

                l12 = norm2(r_rel)
                alpha = 1 / (1 - dt * k_s * (l0 - l12) * mu_perp / (2*l12))

            end function alpha

            ! Second component of matrix inverse of form alpha *I + beta * u u^T
            function beta(r_rel)
                real(wp), dimension(dim), intent(in)     :: r_rel
                real(wp)                                 :: beta, l12, a, b

                l12 = norm2(r_rel)
                a = (1 - dt * k_s * (l0 - l12) * mu_perp / (2*l12))
                b = dt * k_s * (l0 -l12) * (mu_perp - mu_para) / (2*l12)

                beta = -b / (a*(a+b))

            end function beta

        end subroutine anisotropicDimer



          
end module
