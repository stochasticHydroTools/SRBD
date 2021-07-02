program ParticleSpring
    use Precision
    use BoxLibRNGs
    use DiffusionCLs
    !use DiffusionCLs
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit       ! This is so I can track errors.
    implicit none 


    !integer, parameter                                 :: wp = r_sp, dim = 3
    !real(wp), parameter                                :: pi = 4.0_wp*ATAN(1.0_wp), KB = 1.38065E-23_wp, T = 300.0_wp
    real(wp), dimension(dim)                            :: r1, r2, r_cm, r_rel
    integer                                             :: i, myunit3, myunit4, n
    character(len=128)                                  :: filenameDiff = 'diffusion.txt', filenameRel = 'relDist.txt'
    
    ! Uncomment below if want to read from namefile.
    !call read_namelist(nml_file)

    call initializeCLs()

    
    n = 100000

    call SeedRNG(seed) 

    ! Initialize starting positions.
    r1 = 1.0_wp
    r2 = -1.0_wp 
  
    open(newunit = myunit3, file = filenameDiff)
    open(newunit = myunit4, file = filenameRel) 

    r_cm = mu_2 * r1 / mu_eff + mu_1 * r2 / mu_eff
    r_rel = r1-r2  

         
    do i = 0, n

        if (evolve_r_cm) write(myunit1,*) r_cm 
        
        write(myunit4,*) r_rel       ! unit vector if and only if looking at rotational diffusion

        select case (enumer)

        case(1)
            call Euler_Maruyama(dt, nsteps, mu_1, mu_2, k_s, l0, r_cm, r_rel)

        case(2)
            call explicitMidpoint(dt, nsteps, mu_1, mu_2, k_s, l0, r_cm, r_rel)

        case(3)
            call implicitTrapezoidal(dt, nsteps, mu_1, mu_2, k_s, l0, r_cm, r_rel)

        end select
        
    end do

    close(myunit3)
    close(myunit4)       
    
end program