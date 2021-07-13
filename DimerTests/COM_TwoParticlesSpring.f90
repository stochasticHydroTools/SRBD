program ParticleSpring
    use Precision
    use BoxLibRNGs
    use DiffusionCLs
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit       ! This is so I can track errors.
    implicit none 


    integer, parameter                                  :: wp = r_sp, dim = 3
    !real(wp), parameter                                :: pi = 4.0_wp*ATAN(1.0_wp), KB = 1.38065E-23_wp, T = 300.0_wp
    real(wp), dimension(dim)                            :: r1, r2, r_cm, r_rel
    integer                                             :: i, myunit3, myunit4, n, nsteps, numIterates, m
    character(len=128)                                  :: filenameDiff = 'diffusion.txt', filenameRel = 'relDistTest.txt'
    logical, parameter                                  :: rotational = .false., timeScaleSpring = .true., confidenceInterval = .false.
    integer                                             :: seed = 5
    real(wp)                                            :: mu_1, mu_2, dimensionlessTSSize, dt
    real(wp), parameter                                 :: pi = 4.0_wp*ATAN(1.0_wp)
    character(len = 128), dimension(16)                 :: seedFilenames

    call readCLParameters()
    call initializeCLs()

    n = 100000         ! How many data points saved (this is not too related to DiffusionCLs, so am putting it here) 
    dimensionlessTSSize = 1.5_wp   
    nsteps = 1
    numIterates = 1

    seedFilenames = (/'diffusionCSeed05.txt','diffusionCSeed06.txt','diffusionCSeed07.txt', &
    'diffusionCSeed08.txt', 'diffusionCSeed09.txt', 'diffusionCSeed10.txt', &
    'diffusionCSeed11.txt', 'diffusionCSeed12.txt','diffusionCSeed13.txt', &
    'diffusionCSeed14.txt', 'diffusionCSeed15.txt', 'diffusionCSeed16.txt', &
    'diffusionCSeed17.txt', 'diffusionCSeed18.txt', 'diffusionCSeed19.txt', &
    'diffusionCSeed20.txt' /)

    !seedFilenames = (/'rotationalTimeScale_k_5.0_dtau_s_0.1Seed05.txt','rotationalTimeScale_k_5.0_dtau_s_0.1Seed06.txt','rotationalTimeScale_k_5.0_dtau_s_0.1Seed07.txt', &
    !'rotationalTimeScale_k_5.0_dtau_s_0.1Seed08.txt', 'rotationalTimeScale_k_5.0_dtau_s_0.1Seed09.txt', 'rotationalTimeScale_k_5.0_dtau_s_0.1Seed10.txt', &
    !'rotationalTimeScale_k_5.0_dtau_s_0.1Seed11.txt', 'rotationalTimeScale_k_5.0_dtau_s_0.1Seed12.txt','rotationalTimeScale_k_5.0_dtau_s_0.1Seed13.txt', &
    !'rotationalTimeScale_k_5.0_dtau_s_0.1Seed14.txt', 'rotationalTimeScale_k_5.0_dtau_s_0.1Seed15.txt', 'rotationalTimeScale_k_5.0_dtau_s_0.1Seed16.txt', &
    !'rotationalTimeScale_k_5.0_dtau_s_0.1Seed17.txt', 'rotationalTimeScale_k_5.0_dtau_s_0.1Seed18.txt', 'rotationalTimeScale_k_5.0_dtau_s_0.1Seed19.txt', &
    !'rotationalTimeScale_k_5.0_dtau_s_0.1Seed20.txt' /)

    if (confidenceInterval) numIterates = 16

    do m = 0, numIterates - 1
        seed = seed + m
        call SeedRNG(seed) 

        ! Initialize starting positions.
        r1 = 1.0_wp
        r2 = -1.0_wp 
    
        if(confidenceInterval) filenameDiff = seedFilenames(m + 1)


        open(newunit = myunit3, file = filenameDiff)
        open(newunit = myunit4, file = filenameRel) 

        mu_1 = mu_1_0
        mu_2 = mu_2_0

        if (timeScaleSpring) then 
            dt = dimensionlessTSSize * tau_s        ! For spring time scale
        else
            dt = dimensionlessTSSize * tau_r        ! For rotational time scale
        end if
    
        r_cm = mu_2 * r1 / (mu_1 + mu_2) + mu_1 * r2 / (mu_1 + mu_2)
        r_rel = r1-r2  

            
        do i = 0, n

            if (evolve_r_cm) write(myunit3,*) r_cm 
            
            if (rotational) then
                write(myunit4,*) r_rel / norm2(r_rel)     

            else
                write(myunit4,*) r_rel       
            end if


            r1 = r_rel + (r_cm * (mu_1 + mu_2) - mu_2 * r_rel) / (mu_1 + mu_2)
            r2 = (r_cm * (mu_1 + mu_2) - mu_2 * r_rel) / (mu_1 + mu_2)

            call moveDimer(dt, nsteps, mu_1, mu_2, r1, r2)

            r_cm = mu_2 * r1 / (mu_1 + mu_2) + mu_1 * r2 / (mu_1 + mu_2)
            r_rel = r1-r2  
            
        end do

        close(myunit3)
        close(myunit4)   
    end do    
    
    contains

        ! Subroutine uses the exact solution of the Ornstein-Uhlenbeck process to solve the dr_d SODE (not the dr_cm SODE)
        ! This is detailed in the Ornstein-Uhlenbeck wikipedia page, and it states that
        ! r_d = r0 exp( -theta * t) + sigma / sqrt(2*theta)exp(-theta*t)W*
        ! Where W* is a Wiener increment with variance exp(2*theta*t) - 1, and also theta = - mu_eff *k * (l0-l12)/l12
        subroutine exactSol(dt, nsteps, mu_eff, k, D_rel, l0, r_rel)   
            real(wp), intent(in)                        :: dt, mu_eff, k, l0, D_rel
            integer, intent(in)                         :: nsteps
            real(wp), dimension(dim), intent(inout)     :: r_rel

            ! Local Variables
            real(wp), dimension(dim)                    :: disp1
            real(wp)                                    :: l12, theta
            integer                                     :: i

            do i = 1, nsteps
                call NormalRNGVec(numbers = disp1, n_numbers = dim) ! Mean zero and variance one


                l12 = norm2(r_rel)
                theta = -mu_eff * k * (l0 - l12) / l12
                r_rel = r_rel * exp(-theta * dt) + sqrt( D_rel / theta) * exp(-theta * dt) * sqrt(exp(2*theta*dt) - 1) * disp1

            end do
            

        end subroutine
end program