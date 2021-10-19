! This program will sample different dtau (dimensionless time constant) and gather the data for each one. 
! In matlab I find the SD for each of these cases. 
! In particular, we are only looking at relative distances.
! Uses DiffusionCLs module to reduce duplication
! Different from COM_TwoParticlesSpring because this samples many dt.
program diffTimeSteps
    use Precision
    use BoxLibRNGs
    use DiffusionCLs
    implicit none 

    ! Look at time scales LARGER than tau. 

    integer, parameter                                  :: dim_tau = 9 ! Number of time steps sampled
    real(wp), dimension(dim)                            :: r1, r2, r_rel, r_cm
    integer                                             :: n, i, myunit2, l,n_steps 
    character(len=128), dimension(dim_tau)              :: filenameRels
    real(wp), dimension(dim_tau)                        :: timeStepSize = (/ 2.0, 1.0, 2.0**(-1), 2.0**(-2), 2.0**(-3), 2.0**(-4), &
                                                        2.0**(-5), 2.0**(-6), 2.0**(-7) /) 
    ! Other options I have looked at: (/1.0,0.5,0.25,0.20,0.10,0.05/) 
    !(/2.5, 2.0, 1.5, 1.0, 2.0**(-1), 2.0**(-2), 2.0**(-3), 2.0**(-4), 2.0**(-5), 2.0**(-6), 2.0**(-7) /) 
    integer, dimension(dim_tau)                         :: numSteps = (/1, 1, 2, 4, 8, 16, 32, 64, 128/) ! How many time steps in method
    character(len = 128)                                :: filenameRel
    real(wp)                                            :: dimensionless_TSSize, dtau


    !Here is all of the file names I will write out.
    filenameRels = (/'relDistTest_dt01_seed20.txt', 'relDistTest_dt02_seed20.txt', &
    'relDistTest_dt03_seed20.txt', 'relDistTest_dt04_seed20.txt', &
    'relDistTest_dt05_seed20.txt', 'relDistTest_dt06_seed20.txt', &
    'relDistTest_dt07_seed20.txt','relDistTest_dt08_seed20.txt','relDistTest_dt09_seed20.txt'/) 
    !'relDistTest_dt10_seed20.txt','relDistTest_dt11_seed20.txt'/) 

    call initializeCLs()
    n = 100000

    ! Iterate over dt
    do l = 1, dim_tau
        dimensionless_TSSize = timeStepSize(l)
        n_steps = numSteps(l)     ! Used to be 100

        dtau = tau * dimensionless_TSSize   ! Arbitrarily chosen (this is delta t, which is a fraction of tau)
        filenameRel = filenameRels(l)

        call SeedRNG(seed) 

        ! Initialize starting positions.
        r1 = 2.0E-4_wp
        r2 = -2.0E-4_wp 
    
        open(newunit = myunit2, file = filenameRel) 
        
        r_rel = r1-r2  
        r_cm = mu_2 * r1 / mu_eff + mu_1 * r2 / mu_eff

        do i = 0, n

            write(myunit2,*) r_rel

            select case (enumer)

            ! For now, will not use the exact sol since this is only useful for l0 = 0
            case(1)
                call Euler_Maruyama(dt, nsteps, mu_1, mu_2, k_s, l0, r_cm, r_rel)
        
            case(2)
                call explicitMidpoint(dt, nsteps, mu_1, mu_2, k_s, l0, r_cm, r_rel)

            case(3)
                call implicitTrapezoidal(dt, nsteps, mu_1, mu_2, k_s, l0, r_cm, r_rel)
        
            end select  

        end do

        close(myunit2)   
            
    end do

end program
