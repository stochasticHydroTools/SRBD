program relDiff
    use Precision
    use BoxLibRNGs
    implicit none 

    ! Look at times LARGER than tau (to see BM). Can sample each tau. 

    integer, parameter                  :: wp = r_sp, dim = 3, nsteps = 100
    real(wp), parameter                 :: k = 0.029_wp, l0 = 1.0E-8_wp, a = 5.29E-11_wp, pi = 4.0_wp*ATAN(1.0_wp), visc = 8.9E-4_wp  ! Using Bohr radius as estimates, not if l0 order 1e-11, does not satisfy nicely histogram
    real(wp), dimension(dim)            :: r1, r2
    integer n, i
    real(wp), allocatable               :: pos1(:,:), pos2(:,:), pos_rel(:)
    real(wp)                            :: tau, r_rel, dtau
    character(15)                       :: filename = 'data4.txt'
    

    tau = (6*pi*a*visc) / (k)
    dtau = 0.01_wp * tau ! This is arbitrary, but this is the dimensionless time step.

    n = 100000

    ! Initialize starting positions.
    r1 = 0.0_wp
    r2 = 1.20*l0
    r_rel = norm2(r1-r2)

    ! Allocate Memory for all arrays
    allocate(pos1(dim,n + 1))
    allocate(pos2(dim,n + 1))
    allocate(pos_rel(n + 1))

    ! We will continue to update these position arrays and then will print out only the COM one. The others are for tracking data
    pos1(:,1) = r1
    pos2(:,1) = r2
    pos_rel(1) = r_rel
    
    do i = 1, n
        ! Apply one diffusive step to both r1 and r2 simultaneously.
        call Euler_Maruyama(dtau, nsteps, a, visc, k, l0, r1, r2, dim)        ! Just one tau step each iterate.
        pos1(:, i + 1) = r1
        pos2(:, i + 1) = r2
        pos_rel(i + 1) = norm2(r1-r2)

    end do

    call writeVectorToFile(pos_rel, n, filename)

    deallocate(pos1)
    deallocate(pos2)
    deallocate(pos_rel)
        


    contains

        ! r1 is the position of the bead which we are diffusing by one time step. 
        subroutine Euler_Maruyama(dt,nsteps,a,visc,k,l0,r1,r2,dim)
            implicit none
            real(wp), intent(in) :: dt, a, visc, k, l0
            integer, intent(in) :: nsteps, dim
            real(wp), dimension(dim), intent(inout) :: r1, r2
            real(wp), dimension(dim) :: disp1, disp2, temp
            real(wp), parameter :: KB = 1.38065E-23_wp, T = 300.0_wp, pi = 4.0_wp*ATAN(1.0_wp)
            real(wp) mu, l12


            ! Local Variables
            integer :: i

            ! Initialization of constants
            mu = ( 1.0 / (6*pi*visc*a) ) 

            ! Brownian Motion with Deterministic Drift Realization
            do i = 1, nsteps
                call NormalRNGVec(numbers=disp1, n_numbers= dim) ! Mean zero and variance one
                call NormalRNGVec(numbers=disp2, n_numbers= dim) ! Mean zero and variance one

                l12 = norm2(r1-r2)

                temp = r1 - r2

                r1 = r1 + (mu * k * (l12 - l0) * (-temp) / l12) * dt + sqrt(2*KB*T*mu*dt)*disp1 ! Apply one Euler-Maruyama Step   
                r2 = r2 + (mu * k * (l12 - l0) * (temp) / l12) * dt + sqrt(2*KB*T*mu*dt)*disp2

            end do   

        end subroutine      


        ! Write to File
        subroutine writeVectorToFile(pos,n, filename)
            integer, intent(in)                     :: n       ! n represents how many items to print out. usually same as size, but not always
            real(wp), dimension(n+1), intent(in)    :: pos
            character(15)                           :: filename
            integer j

            open(12, file = filename)
            do j = 1, n + 1
                write(12, *) pos(j)
            end do

            close(12)
        

        end subroutine

    
end program
