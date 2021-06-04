! File to diffuse and gather data on TWO SEPERATELY DIFUSING particles, and calculate the centre of mass 
program TwoParticles
    use Precision
    use BoxLibRNGs ! Random number generation
    implicit none

    integer, parameter      :: wp = r_sp

    real(wp)                :: tau = 1.0_wp/(0.1_wp)        ! Arbitraily chosen tau (time scale)
    real(wp), dimension(3)  :: r1 = 0.0_wp, r2 = 1.0_wp     ! Initialize where r1 and r2 start.
    real(wp), allocatable   :: pos1(:,:), pos2(:,:), pos_cm(:,:)
    integer                 :: n , i
    real(wp), dimension(3)  :: r_cm 


    n = 100 ! Number of data entries (not including the starting point)
    allocate(pos1(3,n + 1))
    allocate(pos2(3,n + 1))
    allocate(pos_cm(3,n + 1))


    r_cm = 0.5 * (r1 + r2)

    pos1(:,1)   = r1
    pos2(:,1)   = r2
    pos_cm(:,1) = r_cm

    do i = 1, n
        call brownianStep(r1, r2, pos1, pos2, pos_cm, i, tau)
    end do

    ! Write to File
    call writeToFile(pos_cm,n)

    deallocate(pos1)
    deallocate(pos2)
    deallocate(pos_cm)

    contains

        ! Brownian Step
        subroutine brownianStep(r1, r2, pos1, pos2, pos_cm, i, tau)
            integer , intent(in)                        :: i 
            real(wp), intent(in)                        :: tau
            real(wp), dimension(3)                      :: disp1, disp2, r_cm
            real(wp), dimension(3), intent(inout)       :: r1, r2
            real(wp), dimension(3,n+1), intent(inout)   :: pos1, pos2, pos_cm
            real(wp)                                    :: D = 9 ! Diffusion coefficient of species A

            call NormalRNGVec(numbers=disp1, n_numbers=3) ! Mean zero and variance one
            call NormalRNGVec(numbers=disp2, n_numbers=3) ! Mean zero and variance one

            ! The below two lines would need to be modified once I consider the cases of a spring, or no spring but two partices. 
            disp1 = sqrt(2*D*tau)*disp1 ! Make the variance be 2*D*time
            disp2 = sqrt(2*D*tau)*disp2 ! Make the variance be 2*D*time

            r1 = pos1(:,i) + disp1
            r2 = pos2(:,i) + disp2

            r_cm = 0.5 * (r1 + r2)

            pos1(:,i + 1)   = r1
            pos2(:, i + 1)  = r2
            pos_cm(:,i + 1) = r_cm

        end subroutine
        
        ! Write to File
        subroutine writeToFile(pos,n)
            integer, intent(in) :: n
            real(wp), dimension(3,n+1), intent(in) :: pos
            integer i, j
            

            open(unit=2, file='data2.txt', ACTION="write", STATUS="replace")
            do i=1,3
                 write(2, '(1000F14.7)')( real(pos(i,j)) ,j=1,n+1)
            end do
        
            close(2)
        

        end subroutine
        
end program
