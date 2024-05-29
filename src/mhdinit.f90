module mhdinit
    use parallel
    implicit none

    integer :: nx = 128, ny = 128, nz = 1, nvar = 8, nextern = 1

    real :: pi = 3.141592653589793

    integer :: ifield = 0, ipert = 0

    real :: Lx, Ly, dx, dy

    real :: T0 = 1.0, rho0 = 1.0  ! initial uniform temerature & rho 
    real :: B0 = 1.0, a0 = 0.05 !for current sheet

    real, allocatable, dimension(:) :: xgrid, ygrid

    !logical :: if_isothermal = .false.

    real :: adiabatic_index = 5./3. !, T_isothermal = 1.0
    logical :: if_resis = .false., if_AEB = .false.,if_corotating = .false.,&
        if_Hall = .false., if_resis_exp = .false., if_conserve_background = .false.,&
        if_visc = .false., if_visc_exp = .false., if_z_radial = .false.,&
        if_external_force = .false.
    real :: resistivity = 0.0, viscosity = 0.0, ion_inertial_length
    real :: corotating_angle = 0.0, cos_cor_ang, sin_cor_ang

    real :: time = 0.0

    !uu is the conservation variables:
    !rho, rho*(ux,uy,uz), bx,by,bz, e=(p/(k-1)+0.5*(rho*u^2+B^2))
    !note that in initiliazation, it is treated as primitive variables:
    !rho,ux,uy,uz,bx,by,bz,p
    real, allocatable, dimension(:,:,:,:) :: uu,uu_prim,flux,expand_term,current_density,&
        external_force
    complex, allocatable, dimension(:,:,:,:) :: uu_fourier,flux_fourier,expand_term_fourier,&
        current_density_fourier,external_force_fourier
    

    !arrays storing the wave number information
    real,allocatable,dimension(:) :: wave_number_x, wave_number_y, wave_number_z
    real,allocatable,dimension(:,:,:) :: k_square

    complex, allocatable, dimension(:,:,:,:) :: fnl, fnl_rk




    !user-defined----------------------
    real :: delta_b,db0,in_out_ratio = 1.0, initial_spectral_slope = 2.0, &
        delta_rho, delta_p
    real :: U_fast,U_slow,n_fast,n_slow,press0,shear_width,bx0,by0,bz0,theta
    real :: current_sheet_width = 0.075
    integer :: nmode = 16
    integer :: iMode = 128, mode_start = 64, mode_end = 512

    real :: jet_width=1.0, Vjet = 0.3
    real :: B0_SB,dB_SB,R1_SB,H_SB,Rm_SB,Rd_SB

    contains 

        subroutine grid_initialize
            ! must be called after parallel_initialize
            ! initialize grids and wave-numbers
            implicit none

            integer :: ix,iy,iz
            integer :: ixmin, ixmax, iymin, iymax, izmin, izmax 

            ixmin = xi_offset(myid_i+1) + 1
            ixmax = xi_offset(myid_i+1) + xi_size(myid_i+1)

            iymin = 1
            iymax = ny

            izmin = 1 
            izmax = 1

            allocate(xgrid(nx), ygrid(ny), &
                wave_number_x(nx), wave_number_y(ny),wave_number_z(nz),&
                k_square(ixmin:ixmax,iymin:iymax,izmin:izmax))

            dx = Lx / nx 
            dy = Ly / ny 

            do ix=1,nx 
                xgrid(ix) = 0 + (ix-1) * dx 

                if ( ix <= (nx/2+1) ) then
                    wave_number_x(ix) = 2*pi*(ix-1) / Lx
                else
                    wave_number_x(ix) = 2*pi*(ix-1-nx) / Lx
                endif
            enddo

            do iy=1,ny
                ygrid(iy) = 0 + (iy-1) * dy 

                if ( iy <= (ny/2+1) ) then
                    wave_number_y(iy) = 2*pi*(iy-1) / Ly
                else
                    wave_number_y(iy) = 2*pi*(iy-1-ny) / Ly
                endif
            enddo

            wave_number_z(1) = 0.0

            !define k_square
            do iz=izmin,izmax
                do iy=iymin,iymax 
                    do ix=ixmin,ixmax 
                        k_square(ix,iy,iz) = (wave_number_x(ix))**2 + &
                            (wave_number_y(iy))**2 
                    enddo
                enddo
            enddo
        end subroutine grid_initialize

        subroutine arrays_initialize
            implicit none 
            integer :: ixmin, ixmax, iymin, iymax, izmin, izmax 
            !--------------------------------------------------
            !arrays in real space
            !index range for real-space arrays
            ixmin = 1
            ixmax = nx 

            iymin = yi_offset(myid_i+1) + 1 
            iymax = yi_offset(myid_i+1) + yi_size(myid_i+1)

            izmin = 1 
            izmax = 1

            allocate(uu(ixmin:ixmax, iymin:iymax, izmin:izmax, 1:nvar), &
                uu_prim(ixmin:ixmax, iymin:iymax, izmin:izmax, 1:4), &
                flux(ixmin:ixmax, iymin:iymax, izmin:izmax, 1:18)) 

            if (if_AEB) then 
                allocate(expand_term(ixmin:ixmax, iymin:iymax, izmin:izmax,1))
            endif

            if (if_Hall) then 
                allocate(current_density(ixmin:ixmax, iymin:iymax, izmin:izmax,3))
            endif

            if (if_external_force) then 
                allocate(external_force(ixmin:ixmax, iymin:iymax, izmin:izmax,1:nextern))
            endif 
            !--------------------------------------------------


            !--------------------------------------------------
            !uu in fourier space
            ixmin = xi_offset(myid_i+1) + 1
            ixmax = xi_offset(myid_i+1) + xi_size(myid_i+1)

            iymin = 1
            iymax = ny

            izmin = 1 
            izmax = 1

            allocate(uu_fourier(ixmin:ixmax, iymin:iymax, izmin:izmax, 1:nvar),&
                flux_fourier(ixmin:ixmax, iymin:iymax, izmin:izmax, 1:18), &
                fnl(ixmin:ixmax, iymin:iymax, izmin:izmax, 1:nvar), &
                fnl_rk(ixmin:ixmax, iymin:iymax, izmin:izmax, 1:nvar))

            if (if_AEB) then 
                !allocate(uu_prim_fourier(ixmin:ixmax, iymin:iymax, izmin:izmax, 4:4))
                allocate(expand_term_fourier(ixmin:ixmax, iymin:iymax, izmin:izmax, 1))
            endif

            if (if_Hall) then 
                allocate(current_density_fourier(ixmin:ixmax, iymin:iymax, izmin:izmax, 3))
            endif

            if (if_external_force) then 
                allocate(external_force_fourier(ixmin:ixmax, iymin:iymax, izmin:izmax,1:nextern))
            endif 
            !--------------------------------------------------

        end subroutine arrays_initialize

        subroutine background_fields_initialize
            implicit none 
            real :: kx, ky, kz
            integer :: ix, iy, iz, ivar, ik, &
                ixmin,ixmax,iymin,iymax,izmin,izmax
            real :: ylow,ycent,yup,ninf = 1.0, a_cs 
            real :: ylow0, ylow1, yup0, yup1, ycent0, ycent1, y_cent
            real :: r_yz

            ixmin = 1
            ixmax = nx 

            iymin = yi_offset(myid_i+1) + 1 
            iymax = yi_offset(myid_i+1) + yi_size(myid_i+1)

            izmin = 1 
            izmax = 1

            uu(:,:,:,1) = rho0
            uu(:,:,:,8) = rho0*T0 !note this is pressure

            select case(ifield)

            Case(0)
                !add uniform background magnetic field
                uu(:,:,:,1) = rho0

                uu(:,:,:,5) = bx0 
                uu(:,:,:,6) = by0 
                uu(:,:,:,7) = bz0
                
                uu(:,:,:,8) = press0 
                
            Case(1)
                ! similar to Case(0) but use B0 and angle (x-y plane)
                uu(:,:,:,1) = rho0 

                uu(:,:,:,5) = b0 * cos(theta * PI / 180.)
                uu(:,:,:,6) = b0 * sin(theta * PI / 180.)
                uu(:,:,:,7) = bz0
                
                uu(:,:,:,8) = press0 
                
            case default
                continue
            end select

        end subroutine background_fields_initialize


        subroutine perturbation_initialize
            implicit none
            real :: kx, ky, kz
            integer :: ix, iy, iz, ivar, ik
            integer :: ixmin,ixmax,iymin,iymax,izmin,izmax

            real :: dbx, dby, dbz
            real,allocatable,dimension(:) :: phs, Br_sign
            integer,dimension(12) :: ir_arr 
            integer :: n_m,ir,ikx
            real :: pph,kk,yup,ylow,B_sign
            real :: amplitude_slope_index,ik_slope

            amplitude_slope_index = initial_spectral_slope/2.0

            ixmin = 1
            ixmax = nx 
            iymin = yi_offset(myid_i+1) + 1
            iymax = yi_offset(myid_i+1) + yi_size(myid_i+1)
            izmin = 1
            izmax = 1

            select case(ipert)

            case(1)
                !1D monochromatic, circularly-polarized, alfven wave
                kx =  iMode * 2 * pi / Lx
                do ix=ixmin,ixmax 
                    !uy
                    uu(ix,:,:,3) = uu(ix,:,:,3) + db0/sqrt(uu(ix,:,:,1)) &
                        * cos(kx * xgrid(ix)) 
                    !uz
                    uu(ix,:,:,4) = uu(ix,:,:,4) + db0/sqrt(uu(ix,:,:,1)) & 
                        * sin(kx * xgrid(ix))

                    !by
                    uu(ix,:,:,6) = uu(ix,:,:,6) - db0 * cos(kx * xgrid(ix))
                    !bz
                    uu(ix,:,:,7) = uu(ix,:,:,7) - db0 * sin(kx * xgrid(ix)) 
                enddo

            case(2)
                !2D gaussian perturbation in Pressure
                do iz = izmin,izmax
                    do iy = iymin,iymax
                        do ix = ixmin,ixmax 
                            uu(ix,iy,iz,8) = uu(ix,iy,iz,8) + delta_p * exp( -( (xgrid(ix)-0.5*Lx) &
                                /(0.01*Lx) )**2 - ((ygrid(iy)-0.5*Ly)/(0.01*Ly))**2)
                        enddo
                    enddo
                enddo
            

            case(3)
                !1D shear Alfven wave packet (uz & bz)
                do ix=ixmin,ixmax 
                    !uz
                    uu(ix,:,:,4) = uu(ix,:,:,4) + db0/sqrt(uu(ix,:,:,1)) & 
                        * exp( -( (xgrid(ix)-0.5*Lx)/(0.1*Lx) )**2)
                    !bz
                    uu(ix,:,:,7) = uu(ix,:,:,7) - db0 * exp( &
                        -( (xgrid(ix)-0.5*Lx)/(0.1*Lx) )**2)
                enddo

            case(31)
                !1D monochromatic shear Alfven wave (uz & bz)
                kx =  iMode * 2 * pi / Lx
                do ix=ixmin,ixmax 
                    !uz
                    uu(ix,:,:,4) = uu(ix,:,:,4) + db0/sqrt(uu(ix,:,:,1)) & 
                        * sin(kx * xgrid(ix))
                    !bz
                    uu(ix,:,:,7) = uu(ix,:,:,7) - db0 * sin(kx * xgrid(ix)) 
                enddo

            case(4)
                !1D gaussian pressure peak 
                do ix = ixmin,ixmax 
                    uu(ix,:,:,8) = uu(ix,:,:,8) + delta_p * exp( -( (xgrid(ix)-0.5*Lx) &
                        /(0.01*Lx) )**2)
                enddo

            case(5)
                !2D gaussian perturbation in Bz
                do iz = izmin,izmax
                    do iy = iymin,iymax
                        do ix = ixmin,ixmax 
                            uu(ix,iy,iz,7) = uu(ix,iy,iz,7) + db0 * exp( -( (xgrid(ix)-0.5*Lx) &
                                /(0.01*Lx) )**2 - ((ygrid(iy)-0.5*Ly)/(0.01*Ly))**2)
                        enddo
                    enddo
                enddo

            case(6)
                !2D gaussian perturbation in Bz
                do iz = izmin,izmax
                    do iy = iymin,iymax
                        do ix = ixmin,ixmax 
                            uu(ix,iy,iz,7) = uu(ix,iy,iz,7) + db0 * (xgrid(ix) - 0.5*Lx) * &
                            exp( -( (xgrid(ix)-0.5*Lx)/(0.01*Lx) )**2 - ((ygrid(iy)-0.5*Ly)/(0.01*Ly))**2)
                        enddo
                    enddo
                enddo


            case default
                continue 
            end select


        end subroutine perturbation_initialize

        subroutine initial_calc_conserve_variable
            !after initializing background fields and the perturbation
            !transform uu to be conserved quantities
            implicit none 

            uu_prim(:,:,:,1:3) = uu(:,:,:,2:4) !ux,uy,uz
            uu_prim(:,:,:,4) = uu(:,:,:,8)  !pressure

            !rho * u
            uu(:,:,:,2) = uu(:,:,:,1) * uu_prim(:,:,:,1)
            uu(:,:,:,3) = uu(:,:,:,1) * uu_prim(:,:,:,2)
            uu(:,:,:,4) = uu(:,:,:,1) * uu_prim(:,:,:,3)

            !energy density
            uu(:,:,:,8) = uu_prim(:,:,:,4)/(adiabatic_index-1) + &
                0.5*(uu(:,:,:,1)*( (uu_prim(:,:,:,1))**2 + &
                (uu_prim(:,:,:,2))**2 + (uu_prim(:,:,:,3))**2) + &
                (uu(:,:,:,5))**2 + (uu(:,:,:,6))**2 + &
                (uu(:,:,:,7))**2 )
        end subroutine 

end module mhdinit