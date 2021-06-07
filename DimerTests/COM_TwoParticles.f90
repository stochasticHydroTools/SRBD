! File to diffuse and gather data on TWO SEPERATELY DIFUSING particles, and calculate the centre of mass 
program TwoParticles
    use Precision
    use BoxLibRNGs ! Random number generation
    implicit none

    integer, parameter      :: wp = r_sp, dim = 3

    real(wp)                :: dt = 0.001_wp, D1 = 9, D2 = 9        ! Arbitraily chosen dt 
    real(wp), dimension(dim)  :: r1 = 0.0_wp, r2 = 1.0_wp     ! Initialize where r1 and r2 start.
    real(wp), allocatable   :: pos1(:,:), pos2(:,:), pos_cm(:,:)
    integer                 :: n , i
    real(wp), dimension(dim)  :: r_cm 
    character(15)           :: filename = 'data2.txt'


    n = 10000 ! Number of data entries (not including the starting point)
    allocate(pos1(dim,n + 1))
    allocate(pos2(dim,n + 1))
    allocate(pos_cm(dim,n + 1))


    r_cm = 0.5 * (r1 + r2)

    pos1(:,1)   = r1
    pos2(:,1)   = r2
    pos_cm(:,1) = r_cm

    do i = 1, n
        call brownianStep(r1, r2, r_cm, dt, D1, D2, dim)
        pos1(:,i+1) = r1
        pos2(:,i+1) = r2
        pos_cm(:,i+1) = r_cm
    end do

    ! Write to File
    call writeToFile(pos_cm,n,dim, filename)

    deallocate(pos1)
    deallocate(pos2)
    deallocate(pos_cm)

    contains

        ! Brownian Step
        subroutine brownianStep(r1, r2, r_cm, dt, D1, D2, dim)
            integer, intent(in)                             :: dim
            real(wp), intent(in)                            :: dt, D1, D2
            real(wp), dimension(dim)                        :: disp1, disp2
            real(wp), dimension(dim), intent(inout)         :: r1, r2, r_cm

            call NormalRNGVec(numbers=disp1, n_numbers=3) ! Mean zero and variance one
            call NormalRNGVec(numbers=disp2, n_numbers=3) ! Mean zero and variance one

            ! The below two lines would need to be modified once I consider the cases of a spring, or no spring but two partices. 
            disp1 = sqrt(2*D1*dt)*disp1 ! Make the variance be 2*D*time
            disp2 = sqrt(2*D2*dt)*disp2 ! Make the variance be 2*D*time

            r1 = r1 + disp1
            r2 = r2 + disp2
            r_cm = 0.5 * (r1 + r2)

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
