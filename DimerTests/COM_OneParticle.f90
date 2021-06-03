! File to diffuse and gather data on ONE particle
program test
    use MinHeapModule
    use BoxLibRNGs ! Random number generation
    implicit none

    real :: tau = 1.0/(0.1)
    real, dimension(3) :: r = 0.0
    real, allocatable :: pos(:,:)
    integer :: n , i

    n = 100 ! Number of data entries (not including the starting point)
    allocate(pos(3,n + 1))
    pos(:,1) = r

    do i = 1, n
        call brownianStep(r, pos, i, tau)
    end do

    ! Write to File
    call writeToFile(pos)

    deallocate(pos)

    contains

        ! Brownian Step
        subroutine brownianStep(r, pos, i, tau)
            integer , intent(in) :: i 
            real, intent(in) :: tau
            real disp
            real, dimension(3), intent(inout) :: r 
            real, dimension(3,n+1), intent(inout) :: pos
            real :: D = 9 ! Diffusion coefficient of species A

            call NormalRNGVec(numbers=disp, n_numbers=3) ! Mean zero and variance one
            disp = sqrt(2*D*tau)*disp ! Make the variance be 2*D*time
            pos(:,i + 1) = pos(:,i) + disp

        end subroutine
        
        ! Write to File
        subroutine writeToFile(pos)
            real, dimension(3,n+1), intent(in) :: pos

            open(7,file = 'output.txt')
            write(7,*) pos
            close(7)

        end subroutine
        
end program
