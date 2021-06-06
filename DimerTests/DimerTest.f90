program EM 
    use Precision
    use BoxLibRNGs ! Random number generation
    implicit none   ! Must include

    integer, parameter :: wp = r_sp ! Double precision
    
    real(wp), dimension(3) :: r1, r2
    
    r1 = (/0.0_wp, 0.0_wp, 0.0_wp/)
    r2 = (/0.0_wp, 0.0_wp, 0.5_wp/)

    call Euler_Maruyama(0.05_wp, 100, 1.0_wp, 3.0_wp, 5.0_wp, 1.0_wp, r1, r2)

    print *, r1

    contains 
     
        ! r1 is the position of the bead which we are diffusing by one time step. 
        subroutine Euler_Maruyama(dt,nsteps,a,visc,k,l0,r1,r2)
            implicit none
            real(wp), intent(in) :: dt, a, visc, k, l0
            integer, intent(in) :: nsteps
            real(wp), dimension(3), intent(inout) :: r1, r2
            real(wp), dimension(3) :: disp
            real(wp), parameter :: KB = 1.38065E-23_wp, T = 300.0_wp, pi = 4.0_wp*ATAN(1.0_wp)
            real(wp) mu, l12

            ! Local Variables
            integer :: i

            ! Initialization of constants
            mu = ( 1.0 / (6*pi*visc*a) ) 
                    
            ! Brownian motion (only) realization (code could be modified to be in R3)
            ! do i = 1, nsteps
            !    call NormalRNG(number = disp) ! Mean zero and variance one
            !
                ! Apply one Euler-Maruyama step
            !    x = x + sqrt(2*KB*T*dt)*disp
            !end do   

            ! Brownian Motion with Deterministic Drift Realization
            do i = 1, nsteps
		l12 = NORM2(r1-r2)
                call NormalRNGVec(numbers=disp, n_numbers= 3) ! Mean zero and variance one
                r1 = r1 + (mu * k * (l12 - l0) * (r1-r2) / l12) * dt + sqrt(2*KB*T*mu*dt)*disp ! Apply one Euler-Maruyama Step        
            end do   

        end subroutine      
 
end program