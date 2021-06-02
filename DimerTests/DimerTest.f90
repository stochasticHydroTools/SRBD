program EM 
    use MinHeapModule
    use BoxLibRNGs ! Random number generation
    implicit none
    real, dimension(3) :: r1, r2, x
    r1 = (/0.0, 0.0, 0.0/)
    r2 = (/0.0, 0.0, 0.5/)
    x = (/1.0, 1.0, 1.0/)

    call Euler_Maruyama(0.05, x, 100, 1.0, 3.0, 5.0, 1.0, r1, r2)

    contains 
        subroutine Euler_Maruyama(dt,x,nsteps,a,visc,k,l0,r1,r2)
            real :: dt, disp, a, visc, k, l0, mu, l12
            integer :: nsteps
            real, parameter :: KB = 1.38065*(10**(-23)), T = 300, pi = 4.D0*DATAN(1.D0)
            real, dimension(3), intent(in) :: r1, r2
            real, dimension(3) :: x
        
            ! Local Variables
            integer :: i

            ! Initialization of constants
            mu = ( 1 / (6*pi*visc*a) ) 
            l12 = NORM2(r1-r2)
                    
            ! Brownian motion (only) realization (code could be modified to be in R3)
            ! do i = 1, nsteps
            !    call NormalRNG(number = disp) ! Mean zero and variance one
            !
                ! Apply one Euler-Maruyama step
            !    x = x + sqrt(2*KB*T*dt)*disp
            !end do   

            ! Brownian Motion with Deterministic Drift Realization
            do i = 1, nsteps
                call NormalRNGVec(numbers=disp, n_numbers= 3) ! Mean zero and variance one
                x = x + (mu * k * (l12 - l0) * (r1-r2) / l12) * dt + sqrt(2*KB*T*dt)*disp ! Apply one Euler-Maruyama Step        
            end do   

        end subroutine
 
end program


