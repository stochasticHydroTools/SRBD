program EM 
    use Precision
    use BoxLibRNGs ! Random number generation
    implicit none

    integer, parameter          :: wp = r_sp, dim = 3
    real(wp), parameter         :: dt = 0.01_wp, a = 1.0_wp, visc = 3.0_wp, k = 5.0_wp, l0 = 1.0_wp, pi = 4.0_wp*ATAN(1.0_wp), KB = 1.38065E-23_wp, T = 300.0_wp
    real(wp)                    :: D_d, D_cm, mu_eff
    real(wp)                    :: mu_1 = 1.0_wp, mu_2 = 2.0_wp
    real(wp), dimension(dim)    :: r1, r2, rd, rcm
    integer                     :: n = 100
    
    r1 = (/0.0_wp, 0.0_wp, 0.0_wp/)
    r2 = (/1.5_wp, 1.5_wp, 1.5_wp/)

    D_d = KB * T * (mu_1 + mu_2)

    mu_eff = mu_1 + mu_2
    D_cm = KB * T * mu_1 * mu_2 / (mu_1 + mu_2)


    rd = r1 - r2
    rcm = 0.5 * (r1 + r2)

    ! Will update rd, rcm. 
    !call Euler_Maruyama(dt, n, a, visc, k, l0, r1, r2, dim)

    call implicitMidpoint(dt, n, mu_eff, k, D_cm, D_d, l0, rcm, rd)

    print *, rcm
    print *, rd

    contains 
     
        subroutine Euler_Maruyama(dt, nsteps, mu, k, D, l0, r1, r2)
            real(wp), intent(in)                        :: dt, mu, k, l0, D
            integer, intent(in)                         :: nsteps          
            real(wp), dimension(dim), intent(inout)     :: r1, r2
            
            ! Local variables           
            real(wp), dimension(dim)                    :: disp1, disp2, vel
            real(wp)                                    :: l12, sdev
            integer                                     :: i

            ! Brownian Motion with Deterministic Drift Realization. 
            do i = 1, nsteps
                call NormalRNGVec(numbers = disp1, n_numbers = dim) ! Mean zero and variance one
                call NormalRNGVec(numbers = disp2, n_numbers = dim) ! Mean zero and variance one

                l12 = norm2(r1-r2) 

                sdev = sqrt(2 * D * dt)


                ! If > 1 dimensions:
                if (dim /= 1) then 
                    vel = r1 - r2
                    vel = mu * k * (l12 - l0) * vel / l12

                    r1 = r1 - vel * dt + sdev*disp1 ! Apply one Euler-Maruyama Step to both r1, r2.
                    r2 = r2 + vel * dt + sdev*disp2 

                ! If 1 Dimension
                else 
                    r1 = r1 + mu * k * (r2 - r1 - l0) * dt + sdev*disp1
                    r2 = r2 + mu * k * (r1 - r2 - l0) * dt + sdev*disp2
                end if

            end do

        end subroutine
        
        
        subroutine explicitMidpoint(dt, nsteps, mu_eff, k, Dcm, D_d, l0, rcm, r_rel)
            real(wp), intent(in)                        :: dt, mu_eff, k, l0, Dcm, D_d
            integer, intent(in)                         :: nsteps
            real(wp), dimension(dim), intent(inout)     :: rcm, r_rel

            ! Local Variables
            real(wp), dimension(dim)                    :: disp1, disp2, rcm_pred, r_rel_pred
            real(wp)                                    :: l12, sdev_cm, sdev_d
            integer                                     :: i

            !mu_eff = mu1 * mu2 / (mu1 + mu2)

            ! Explicit Midpoint loop
            do i = 1, nsteps
                call NormalRNGVec(numbers = disp1, n_numbers = dim) ! Mean zero and variance one
                call NormalRNGVec(numbers = disp2, n_numbers = dim) ! Mean zero and variance one

                l12 = norm2(r_rel)
                sdev_cm = sqrt(2 * Dcm)
                sdev_d = sqrt(2 * D_d)

                ! Will want to apply one Explicit Midpoint step on BOTH rcm and rd

                ! COM STEP
                ! Predictor Step
                rcm_pred = rcm + sqrt(dt/2) * sdev_cm * disp1
                ! Corrector Step
                rcm = rcm + sqrt(dt/2) * sdev_cm * (disp1 + disp2)

                ! R DIFFERENCE STEP
                l12 = norm2(r_rel) ! Evaluate l12 when we are at x_n. Note that this is DEPENDANT on where rd is so l12 = l12(rd). Make sure to update appropriately
                !Predictor Step
                r_rel_pred = r_rel + dt * mu_eff * k * r_rel * (l0 - l12) / (l12 * 2) + sqrt(dt / 2) * sdev_d * disp1
                ! Corrector L has now changed, as l12 evaluated at the predictor stage is now different. 
                l12 =norm2(r_rel_pred) ! L evaluated at n + 1/2
                r_rel = r_rel + dt * mu_eff * k * r_rel_pred * (l0 - l12) / (l12) + sqrt(dt / 2) * sdev_d * (disp1 + disp2)
                

            end do
        end subroutine



        subroutine implicitMidpoint(dt, nsteps, mu_eff, k, Dcm, D_d, l0, rcm, rd)
            real(wp), intent(in)                        :: dt, mu_eff, k, l0, Dcm, D_d
            integer, intent(in)                         :: nsteps
            real(wp), dimension(dim), intent(inout)     :: rcm, rd

            ! Local Variables
            real(wp), dimension(dim)                    :: disp1, disp2, rcm_pred, r_d_pred
            real(wp)                                    :: l12, sdev_cm, sdev_d, L2_n, L_n
            integer                                     :: i

            !mu_eff = mu1 * mu2 / (mu1 + mu2)

            ! Implicit Midpoint loop
            do i = 1, nsteps
                call NormalRNGVec(numbers = disp1, n_numbers = dim) ! Mean zero and variance one
                call NormalRNGVec(numbers = disp2, n_numbers = dim) ! Mean zero and variance one

                sdev_cm = sqrt(2 * Dcm)
                sdev_d = sqrt(2 * D_d)

                ! Will want to apply one Implicit Midpoint step on BOTH rcm and rd

                ! COM STEP
                ! Predictor Step
                rcm_pred = rcm + sqrt(dt) * sdev_cm * disp1
                !   Corrector Step
                rcm = rcm + (sqrt(dt) / 2) * (sdev_cm + sdev_cm) * disp1

                ! R DIFFERENCE STEP
                l12 = norm2(rd) ! Evaluate l12 when we are at x_n. Note that this is DEPENDANT on where rd is so l12 = l12(rd). Make sure to update appropriately
                L_n = mu_eff * k * (l0 - l12) / l12

                !Predictor Step
                r_d_pred = rd + (dt / 2) * L_n * rd + sqrt(dt) * sdev_d * disp1
                r_d_pred = r_d_pred / (1 - dt * L_n / 2)

                ! Corrector L has now changed, as l12 evaluated at the predictor stage is now different. 

                rd = rd + (dt/2) * L_n * rd + (sqrt(dt) / 2) * (sdev_d + sdev_d) * disp1

                l12 = norm2(r_d_pred)  ! For evaluation of L_n+1
                L2_n = mu_eff * k * (l0 - l12) / l12

                rd = rd / (1 - dt * L2_n / 2)
            end do



        end subroutine
 
end program