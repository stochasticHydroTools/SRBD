program relDiff
    use Precision
    use BoxLibRNGs
    implicit none 

    ! Look at time scales LARGER than tau. 

    integer, parameter                  :: wp = r_sp
    real(wp), parameter                 :: k = 0.029_wp, l0 = 1.0E-8, a = 5.29E-11_wp, pi = 4.0_wp*ATAN(1.0_wp), visc = 8.9E-4_wp  ! Using Bohr radius as estimates
    real(wp), dimension(3)              :: r1, r2
    integer n, i
    real(wp), allocatable               :: pos1(:,:), pos2(:,:), pos_rel(:)
    real(wp)                            :: tau, r_rel
    

    tau = (6*pi*a*visc) / (k)

    n = 100000

    ! Initialize starting positions.
    r1 = 0.0_wp
    r2 = 1.20*l0
    r_rel = norm2(r1-r2)

    ! Allocate Memory for all vectors
    allocate(pos1(3,n + 1))
    allocate(pos2(3,n + 1))
    allocate(pos_rel(n + 1))

    ! We will continue to update these position arrays and then will print out only the COM one. The others are for tracking data
    pos1(:,1) = r1
    pos2(:,1) = r2
    pos_rel(1) = r_rel
    
    do i = 1, n
        ! Apply one diffusive step to both r1 and r2 simultaneously.
        call Euler_Maruyama(tau * 0.01_wp, 100, a, visc, k, l0, r1, r2)        ! Just one tau step each iterate.
        pos1(:, i + 1) = r1
        pos2(:, i + 1) = r2
        pos_rel(i + 1) = norm2(r1-r2)

    end do

    call writeToFile(pos_rel, n)

    deallocate(pos1)
    deallocate(pos2)
    deallocate(pos_rel)
        


    contains

        ! r1 is the position of the bead which we are diffusing by one time step. 
        subroutine Euler_Maruyama(dt,nsteps,a,visc,k,l0,r1,r2)
            implicit none
            real(wp), intent(in) :: dt, a, visc, k, l0
            integer, intent(in) :: nsteps
            real(wp), dimension(3), intent(inout) :: r1, r2
            real(wp), dimension(3) :: disp1, disp2, temp
            real(wp), parameter :: KB = 1.38065E-23_wp, T = 300.0_wp, pi = 4.0_wp*ATAN(1.0_wp)
            real(wp) mu, l12


            ! Local Variables
            integer :: i

            ! Initialization of constants
            mu = ( 1.0 / (6*pi*visc*a) ) 

            ! Brownian Motion with Deterministic Drift Realization
            do i = 1, nsteps
                call NormalRNGVec(numbers=disp1, n_numbers= 3) ! Mean zero and variance one
                call NormalRNGVec(numbers=disp2, n_numbers= 3) ! Mean zero and variance one

                l12 = norm2(r1-r2)

                if (norm2(sqrt(2*KB*T*mu*dt)*disp1) > 1) then 
                    print *, "Large"
                endif

                if (norm2(sqrt(2*KB*T*mu*dt)*disp2) > 1.0_wp) then 
                    print *, "Large2"
                endif

                temp = r1 - r2

                !if (norm2(temp/l12) /= 1.0_wp) then 
                !    print *, norm2(temp/l12)
                !endif

                r1 = r1 + (mu * k * (l12 - l0) * (-temp) / l12) * dt + sqrt(2*KB*T*mu*dt)*disp1 ! Apply one Euler-Maruyama Step   
                r2 = r2 + (mu * k * (l12 - l0) * (temp) / l12) * dt + sqrt(2*KB*T*mu*dt)*disp2

            end do   

        end subroutine      


        ! Write to File
        subroutine writeToFile(pos,n)
            integer, intent(in) :: n
            real(wp), dimension(n+1), intent(in) :: pos
            integer j

            open(12, file = 'data4.txt')
            do j = 1, n + 1
                write(12, *) pos(j)
            end do

            close(12)

        end subroutine

    
end program