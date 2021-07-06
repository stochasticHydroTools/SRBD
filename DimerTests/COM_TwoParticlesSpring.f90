program ParticleSpring
    use Precision
    use BoxLibRNGs
    use DiffusionCLs
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit       ! This is so I can track errors.
    implicit none 


    integer, parameter                                  :: wp = r_sp, dim = 3
    !real(wp), parameter                                :: pi = 4.0_wp*ATAN(1.0_wp), KB = 1.38065E-23_wp, T = 300.0_wp
    real(wp), dimension(dim)                            :: r1, r2, r_cm, r_rel
    integer                                             :: i, myunit3, myunit4, n, nsteps
    character(len=128)                                  :: filenameDiff = 'diffusion.txt', filenameRel = 'relDist.txt'
    logical, parameter                                  :: rotational = .false.
    integer                                             :: seed = 5
    real(wp)                                            :: mu_1, mu_2, tau, dimensionlessTSSize, dt

    call initializeCLs()


    n = 500000          ! How many data points saved (this is not too related to DiffusionCLs, so am putting it here) 
    dimensionlessTSSize = 0.01   
    nsteps = 100 


    call SeedRNG(seed) 

    ! Initialize starting positions.
    r1 = -1.0_wp
    r2 = 1.0_wp 
  
    open(newunit = myunit3, file = filenameDiff)
    open(newunit = myunit4, file = filenameRel) 

    mu_1 = 1.0_wp
    mu_2 = 1.0_wp

    tau = 1 / (k_s * (mu_1 + mu_2))
    dt = dimensionlessTSSize * tau


    r_cm = mu_2 * r1 / (mu_1 + mu_2) + mu_1 * r2 / (mu_1 + mu_2)
    r_rel = r1-r2  

         
    do i = 0, n

        if (evolve_r_cm) write(myunit3,*) r_cm 
        
        if (rotational) then
            write(myunit4,*) r_rel / norm2(r_rel)     

        else
            write(myunit4,*) r_rel       
        end if

        select case (sde_integrator_enum)

        ! For now, will not use the exact sol since this is only useful for l0 = 0
        case(1)
            call eulerMaruyama(dt, nsteps, mu_1, mu_2, r_cm, r_rel)

        case(2)
            call explicitMidpoint(dt, nsteps, mu_1, mu_2, r_cm, r_rel)

        case(3)
            call implicitTrapezoidal(dt, nsteps, mu_1, mu_2, r_cm, r_rel)

        end select
        
    end do

    close(myunit3)
    close(myunit4)       
    
    contains

        ! Subroutine uses the exact solution of the Ornstein-Uhlenbeck process to solve the dr_d SODE (not the dr_cm SODE)
        ! This is detailed in the Ornstein-Uhlenbeck wikipedia page, and it states that
        ! r_d = r0 exp( -theta * t) + sigma / sqrt(2*theta)exp(-theta*t)W*
        ! Where W* is a Wiener increment with variance exp(2*theta*t) - 1, and also theta = - mu_eff *k * (l0-l12)/l12
        subroutine exactSol(dt, nstepsCLs, mu_eff, k, D_rel, l0, r_rel)   
            real(wp), intent(in)                        :: dt, mu_eff, k, l0, D_rel
            integer, intent(in)                         :: nstepsCLs
            real(wp), dimension(dim), intent(inout)     :: r_rel

            ! Local Variables
            real(wp), dimension(dim)                    :: disp1
            real(wp)                                    :: l12, theta
            integer                                     :: i

            do i = 1, nstepsCLs
                call NormalRNGVec(numbers = disp1, n_numbers = dim) ! Mean zero and variance one


                l12 = norm2(r_rel)
                theta = -mu_eff * k * (l0 - l12) / l12
                r_rel = r_rel * exp(-theta * dt) + sqrt( D_rel / theta) * exp(-theta * dt) * sqrt(exp(2*theta*dt) - 1) * disp1

            end do
            

        end subroutine
end program