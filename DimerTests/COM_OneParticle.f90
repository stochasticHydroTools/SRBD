! File to diffuse and gather data on ONE particle
program OneParticle
    use Precision
    use BoxLibRNGs ! Random number generation
    implicit none

    integer, parameter      :: wp = r_sp, dim = 3
    real(wp), parameter     :: dt = 0.001_wp, D = 9.0_wp
    real(wp), dimension(dim)  :: r = 0.0_wp
    real(wp), allocatable   :: pos(:,:)
    integer                 :: n , i
    character(15)           :: filename = 'data1.txt'

    n = 10000 ! Number of data entries (not including the starting point) -- will be large

    allocate(pos(dim,n + 1))
    pos(:,1) = r

    do i = 1, n
        call brownianStep(r, i, dt,dim, D)
        pos(:,i+1) = r
    end do

    ! Write to File
    call writeToFile(pos,n,dim, filename)

    deallocate(pos)

    contains

        ! Brownian Step
        subroutine brownianStep(r, i, dt,dim, D)
            integer , intent(in)                            :: i, dim
            real(wp), intent(in)                            :: dt, D
            real(wp), dimension(dim)                        :: disp
            real(wp), dimension(dim), intent(inout)     :: r

            call NormalRNGVec(numbers=disp, n_numbers=dim) ! Mean zero and variance one

            ! The below two lines would need to be modified once I consider the cases of a spring, or no spring but two partices. 
            disp = sqrt(2*D*dt)*disp ! Make the variance be 2*D*time
            r = r + disp

        end subroutine
        
        ! Write to File
        subroutine writeToFile(pos,n,dim, filename)
            integer, intent(in)                     :: n, dim
            real(wp), dimension(3,n+1), intent(in)  :: pos
            character(15)                           :: filename
            integer j, i

            open(12, file = filename)
            do j = 1, n + 1
                write(12, *) (pos(i,j), i = 1, dim)
            end do

            close(12)
        

        end subroutine
        
end program
