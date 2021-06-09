! File to diffuse and gather data on ONE particle
program OneParticle
    use Precision
    use BoxLibRNGs
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit       ! This is so I can track errors.
    implicit none 

    integer, parameter          :: wp = r_sp, dim = 3
    real(wp)                    :: dt, D
    real(wp), dimension(dim)    :: r = 0.0_wp
    real(wp), allocatable       :: pos(:,:)
    integer                     :: n , i, unit1
    character(len = 128)        :: filename = 'data1.txt', nml_file = 'brownianParticle.nml'

    call read_namelist(nml_file, dt, D, n)

    open(newunit = unit1, file = filename)

    do i = 1, n
        write(unit1,*) r(:)
        call brownianStep(r, dt,dim, D)
    end do

    write(unit1,*) r(:)

    close(unit1)

    contains

        ! Brownian Step
        subroutine brownianStep(r, dt,dim, D)
            integer , intent(in)                            :: dim
            real(wp), intent(in)                            :: dt, D
            real(wp), dimension(dim)                        :: disp
            real(wp), dimension(dim), intent(inout)         :: r

            call NormalRNGVec(numbers=disp, n_numbers=dim) ! Mean zero and variance one

            ! The below two lines would need to be modified once I consider the cases of a spring, or no spring but two partices. 
            disp = sqrt(2*D*dt)*disp ! Make the variance be 2*D*time
            r = r + disp

        end subroutine
        
        ! Subroutine to read in the data from a namelist file. 
        subroutine read_namelist(file_path, dt, D, n)
            ! Reads Namelist from given file.
            character(len=*),  intent(in)  :: file_path
            integer,           intent(out) :: n
            real,              intent(out) :: dt, D
            integer                        :: unit, check
    
            ! Namelist definition.
            namelist /brownianParticle/ dt, D, n
    
            ! Check whether file exists.
            inquire (file=file_path, iostat=check)
    
            ! Here we have some checks, this just makes sure that the file is around. 
            if (check /= 0) then
                write (stderr, '(3a)') 'Error: The file "', trim(file_path), '" does not exist.'
                return
            end if
    
            ! Open and read Namelist file.
            open (action='read', file=file_path, iostat=check, newunit=unit)
            read (nml=brownianParticle, iostat=check, unit=unit)
    
            ! This is to keep people aware of cases like : End of File runtime errors.
            if (check /= 0) then
                write (stderr, '(a)') 'Error: invalid Namelist format.'
            end if
    
            close (unit)
        end subroutine read_namelist
        
end program
