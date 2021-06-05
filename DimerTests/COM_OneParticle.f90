! File to diffuse and gather data on ONE particle
program OneParticle
    use Precision
    use BoxLibRNGs ! Random number generation
    implicit none

    integer, parameter      :: wp = r_sp
    real, parameter         :: pi = 4.0_wp*ATAN(1.0_wp)
    real(wp)                :: tau
    real(wp), dimension(3)  :: r = 0.0_wp
    real(wp), allocatable   :: pos(:,:)
    integer                 :: n , i

    n = 10000 ! Number of data entries (not including the starting point) -- will be large
    tau = 0.001_wp      ! Should be very small. 

    allocate(pos(3,n + 1))
    pos(:,1) = r

    do i = 1, n
        call brownianStep(pos, i, tau)
    end do

    ! Write to File
    call writeToFile(pos,n)

    deallocate(pos)

    contains

        ! Brownian Step
        subroutine brownianStep(pos, i, tau)
            integer , intent(in)                        :: i 
            real(wp), intent(in)                        :: tau
            real(wp), dimension(3)                      :: disp
            real(wp), dimension(3,n+1), intent(inout)   :: pos
            real(wp)                                    :: D = 9 ! Diffusion coefficient of species A

            call NormalRNGVec(numbers=disp, n_numbers=3) ! Mean zero and variance one

            ! The below two lines would need to be modified once I consider the cases of a spring, or no spring but two partices. 
            disp = sqrt(2*D*tau)*disp ! Make the variance be 2*D*time
            pos(:,i + 1) = pos(:,i) + disp

        end subroutine
        
        ! Write to File
        subroutine writeToFile(pos,n)
            integer, intent(in) :: n
            real(wp), dimension(3,n+1), intent(in) :: pos
            integer j

            open(12, file = 'data1.txt')
            do j = 1, n + 1
                write(12, *) pos(1,j), pos(2,j), pos(3,j)
            end do

            close(12)
            

          !  open(unit=2, file='data1.txt', ACTION="write", STATUS="replace")
           ! do i=1,3
           !      write(2, '(1000F14.7)')( real(pos(i,j)) ,j=1,n+1)
           ! end do
        
            !close(2)
        

        end subroutine
        
end program
