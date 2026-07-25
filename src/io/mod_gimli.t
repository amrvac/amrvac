!> Module for reading in Legolas .ldat data as initial condition
module mod_gimli
    use, intrinsic :: iso_fortran_env
    use mod_global_parameters
    use mod_hd_phys, only: rho_hd=>rho_, mom_hd=>mom, p_hd=>p_
    use mod_mhd_phys, only: rho_mhd=>rho_, mom_mhd=>mom, p_mhd=>p_ 
    use mod_functions_bfield, only: mag_mhd=>mag

    implicit none
    public

    integer, parameter :: dp = real64
    integer :: idirmin0, idirmin

    complex(dp), parameter :: ic = (0.0d0, 1.0d0)

    integer :: ef_gridpts
    real(dp) :: k2, k3
    integer :: mhd_bool
    real(dp), allocatable :: ef_grid(:)
    complex(dp), allocatable :: rho(:)
    complex(dp), allocatable :: v1(:)
    complex(dp), allocatable :: v2(:)
    complex(dp), allocatable :: v3(:)
    complex(dp), allocatable :: p(:)
    complex(dp), allocatable :: B1(:)
    complex(dp), allocatable :: B2(:)
    complex(dp), allocatable :: B3(:)

    integer :: rho_idx, p_idx
    integer, allocatable :: mom_idx(:), mag_idx(:)

contains

    subroutine read_legolas_data(legolas_file, file_id)
        character(len=100), intent(in) :: legolas_file
        integer, intent(in) :: file_id

        open( &
        unit=file_id+mype, &
        file=legolas_file, &
        form='unformatted' &
        )

        ! need to read physics_type because mod_physics not yet loaded at this point
        read(file_id+mype) mhd_bool
        read(file_id+mype) ef_gridpts
        read(file_id+mype) k2, k3

        call allocate_arrays(ef_gridpts)
        read(file_id+mype) ef_grid
        read(file_id+mype) rho
        read(file_id+mype) v1
        read(file_id+mype) v2
        read(file_id+mype) v3
        read(file_id+mype) p
        if (mhd_bool == 1) then
            read(file_id+mype) B1
            read(file_id+mype) B2
            read(file_id+mype) B3
        end if

        if (mhd_bool == 1) then 
            read(file_id+mype) unit_length, unit_numberdensity, unit_temperature, &
            unit_density, unit_pressure, unit_velocity, unit_magneticfield, unit_time
        else 
            read(file_id+mype) unit_length, unit_numberdensity, unit_temperature, &
            unit_density, unit_pressure, unit_velocity, unit_time
        end if

        close(file_id+mype)
    end subroutine read_legolas_data

    subroutine allocate_arrays(gridpts)
        integer, intent(in) :: gridpts

        allocate(ef_grid(gridpts))
        allocate(rho(gridpts))
        allocate(v1, v2, p, mold=rho)
        if (ndim >= 3) allocate(v3, mold=rho)
        if (mhd_bool == 1) then
            allocate(B1, B2, mold=rho)
            if (ndim >= 3) allocate(B3, mold=rho)
        end if
    end subroutine allocate_arrays

    subroutine add_perturbation_to_w_array(ixI^L, ixO^L, w, w_index, x)
        integer, intent(in)     :: ixI^L, ixO^L
        real(dp), intent(inout) :: w(ixI^S, nw)
        integer, intent(in)     :: w_index
        real(dp), intent(in)    :: x(ixI^S, ndim)
        complex(dp) :: amplitude, exp_factor(ixI^S), quantity, values(ef_gridpts)
        integer  :: idx^D
        real(dp) :: k1 = 0.0d0

        call w_index_to_array(w_index, values)
        exp_factor = {exp(ic * k^D * x(ixI^S, ^D))|*}

        {do idx^D = ixImin^D, ixImax^D|\}
        call ef_amplitude(x(idx^D, 1), ef_grid, values, amplitude)
        quantity = amplitude * exp_factor(idx^D)
        w(idx^D, w_index) = w(idx^D, w_index) + quantity%re
        {end do|\}
    end subroutine add_perturbation_to_w_array

    subroutine init_gimli_indices()
        if (mhd_bool == 1) then
            rho_idx = rho_mhd
            p_idx   = p_mhd

            if (.not. allocated(mom_idx)) allocate(mom_idx(size(mom_mhd)))
            mom_idx = mom_mhd

            if (.not. allocated(mag_idx)) allocate(mag_idx(size(mag_mhd)))
            mag_idx = mag_mhd
        else
            rho_idx = rho_hd
            p_idx   = p_hd

            if (.not. allocated(mom_idx)) allocate(mom_idx(size(mom_hd)))
            mom_idx = mom_hd

        end if
    end subroutine init_gimli_indices

    subroutine w_index_to_array(w_index, array)
        integer, intent(in) :: w_index
        complex(dp), intent(inout) :: array(ef_gridpts)

        call init_gimli_indices()

        if (w_index == rho_idx) then
            array = rho
        else if (w_index == mom_idx(1)) then
            array = v1
        else if (w_index == mom_idx(2)) then
            array = v2
        end if

        if (ndim >= 3) then        
            if (w_index == mom_idx(3)) then
                array = v3
            end if
        end if
        if (w_index == p_idx) then
            array = p
        end if
        if (mhd_bool == 1) then
            if (w_index == mag_idx(1)) then
                array = B1
            else if (w_index == mag_idx(2)) then
                array = B2
            end if
            if (ndim >= 3) then
                if (w_index == mag_idx(3)) then
                    array = B3
                end if
            end if
        end if 
    end subroutine w_index_to_array

    subroutine ef_amplitude(x, grid, array, amplitude)
        use mod_geometry, only: coordinate, cylindrical
        real(dp), intent(in) :: x, grid(ef_gridpts)
        complex(dp), intent(in) :: array(ef_gridpts)
        complex(dp), intent(out) :: amplitude
        integer :: idl, idu
        real(dp) :: alpha

        if (x <= grid(1) .and. xprobmin1 <= x) then
            amplitude = array(1)
            
            if (coordinate==cylindrical .and. grid(1) > 0) then
                if (k2 .ne. 0.0d0) then
                    alpha = abs(k2)
                else
                    alpha = 2.0d0
                end if
                amplitude = array(1) * (x / grid(1))**alpha
            end if
        else if (x >= grid(size(grid)) .and. x <= xprobmax1) then
            amplitude = array(size(grid))
        else if (x < xprobmin1 .or. x > xprobmax1) then
            amplitude = 0.d0
        else
            idl = maxloc(grid, mask=(grid < x), dim=1)
            idu = minloc(grid, mask=(grid > x), dim=1)
            amplitude = array(idl) + (x - grid(idl)) * &
                (array(idu) - array(idl)) / (grid(idu) - grid(idl))
        end if
    end subroutine ef_amplitude

    ! subroutine custom analytics log file
    ! called in "usr_print_log => analytics_log" in mod_usr.t
    subroutine analytics_log
        use mod_global_parameters
        use mod_input_output, only: get_volume_average_func, printlog_default

        integer, parameter :: n_modes = 2
        integer, parameter :: my_unit = 123
        character(len=80)  :: fmt_string = '(9(es12.4))', filename
        logical, save      :: alive, visited=.false.
        double precision   :: volume, magn_avg
        double precision   :: Tmax, Tmin, vmax, B1max, B2max, B3max

        ! First output the standard log file:
        call printlog_default

        ! Now make the custom _c.log file:
        call get_minmax_temperature(Tmax,Tmin)
        call get_max_velocity(vmax)
        if (mhd_bool == 1) then
        call get_max_B(B1max, B2max, B3max)
            call get_volume_average_func(magnetic, magn_avg, volume)
        end if

        filename = trim(base_filename) // "_c.log"

        if (.not. visited) then
        ! Delete the log when not doing a restart run
            if (restart_from_file == undefined .or. reset_time) then
                open(unit=my_unit,file=trim(filename),form='formatted',status='replace')
                if (mhd_bool == 1) then
                    write(my_unit,'(a)') 'global_time Tmax Tmin vmax B1max B2max B3max mag_avg'
                else
                    write(my_unit,'(a)') 'global_time Tmax Tmin vmax'
                end if
            end if
            visited = .true.
        end if

        if (mype == 0) then
        write(filename,"(a)") filename
        inquire(file=filename,exist=alive)
        if(alive) then
            open(unit=my_unit,file=filename,form='formatted',status='old',access='append')
        else
            open(unit=my_unit,file=filename,form='formatted',status='new')
        endif

        ! if number of output doubles is increase, don't forget to change the fmt_string above
        if (mhd_bool == 1) then
            write(my_unit, fmt_string) global_time, Tmax, Tmin, vmax, B1max, B2max, B3max, magn_avg
        else
            write(my_unit, fmt_string) global_time, Tmax, Tmin, vmax
        end if
        close(my_unit)
        end if
    end subroutine analytics_log

    pure function magnetic(w_vec, w_size) result(magn_energy)
        use mod_global_parameters
        integer, intent(in)          :: w_size
        double precision, intent(in) :: w_vec(w_size)
        double precision             :: magn_energy

        magn_energy = 0.5d0 * sum(w_vec(mag_mhd(:))**2,dim=ndir+1)
    end function magnetic

    ! Calculate both min and max of temperature on grid in one go.
    subroutine get_minmax_temperature(Tmax, Tmin)
        use mod_global_parameters
        use mod_physics, only: phys_get_pthermal, phys_get_rho

        double precision, intent(out) :: Tmax, Tmin

        integer                       :: iigrid, igrid, iw
        double precision              :: Tmax_mype, Tmax_recv, Tmin_mype, Tmin_recv
        double precision              :: w(ixG^T,1:nw), wlocal(ixG^T,1:nw), xlocal(ixG^T,1:nw)
        double precision              :: Te(ixG^T), pth(ixG^T), rho(ixG^T)

        Tmax_mype = -bigdouble
        Tmin_mype = bigdouble

        !Loop over all the grids
        do iigrid = 1, igridstail
            igrid = igrids(iigrid)

            wlocal(ixG^T,1:nw) = ps(igrid)%w(ixG^T,1:nw)
            xlocal(ixG^T,1:ndim) = ps(igrid)%x(ixG^T,1:ndim)
            call phys_get_pthermal(wlocal,xlocal,ixG^LL,ixG^LL,pth)
            call phys_get_rho(wlocal,xlocal,ixG^LL,ixG^LL,rho)
            Te(ixG^T) = pth(ixG^T)/rho(ixG^T)

            ! Compare values on current grid to temporary max/min
            Tmax_mype = max(Tmax_mype,maxval(Te(ixM^T)))
            Tmin_mype = min(Tmin_mype,minval(Te(ixM^T)))
        end do

        ! Make the information available on all tasks
        call mpi_allreduce(Tmax_mype, Tmax_recv, 1, mpi_double_precision, &
            mpi_max, icomm, ierrmpi)
        call mpi_allreduce(Tmin_mype, Tmin_recv, 1, mpi_double_precision, &
            mpi_min, icomm, ierrmpi)

        Tmax = Tmax_recv
        Tmin = Tmin_recv

    end subroutine get_minmax_temperature

    subroutine get_max_velocity(vmax)
        use mod_global_parameters
        use mod_physics, only: phys_get_v

        double precision, intent(out) :: vmax

        integer                       :: iigrid, igrid, iw
        double precision              :: vmax_mype,vmax_recv
        double precision              :: v_vec(ixG^T,1:ndir), v(ixG^T)

        vmax_mype = -bigdouble

        !Loop over all the grids
        do iigrid = 1, igridstail
            igrid = igrids(iigrid)

            call phys_get_v(ps(igrid)%w, ps(igrid)%x, ixG^LL, ixG^LL, v_vec)

            v(ixM^T) = sqrt(sum(v_vec(ixM^T,:)**2,dim=ndir+1))

            vmax_mype =max(vmax_mype ,maxval(v(ixM^T)))

        end do

        ! Make the information available on all tasks
        call mpi_allreduce(vmax_mype, vmax_recv, 1, mpi_double_precision, &
            mpi_max, icomm, ierrmpi)

        vmax = vmax_recv

    end subroutine get_max_velocity

    subroutine get_max_B(B1max, B2max, B3max)
        use mod_global_parameters

        double precision, intent(out) :: B1max, B2max, B3max

        integer                       :: iigrid, igrid, iw
        double precision              :: B1max_mype, B1max_recv
        double precision              :: B2max_mype, B2max_recv
        double precision              :: B3max_mype, B3max_recv

        B1max_mype = -bigdouble
        B2max_mype = -bigdouble
        B3max_mype = -bigdouble

        !Loop over all the grids
        do iigrid = 1, igridstail
            igrid = igrids(iigrid)

            B1max_mype = max(B1max_mype, maxval(ps(igrid)%w(ixG^T, mag_mhd(1))))
            B2max_mype = max(B2max_mype, maxval(ps(igrid)%w(ixG^T, mag_mhd(2))))
            if (ndir >= 3) then
                B3max_mype = max(B3max_mype, maxval(ps(igrid)%w(ixG^T, mag_mhd(3))))
            end if
        end do

        ! Make the information available on all tasks
        call mpi_allreduce(B1max_mype, B1max_recv, 1, mpi_double_precision, &
            mpi_max, icomm, ierrmpi)
        call mpi_allreduce(B2max_mype, B2max_recv, 1, mpi_double_precision, &
            mpi_max, icomm, ierrmpi)
        call mpi_allreduce(B3max_mype, B3max_recv, 1, mpi_double_precision, &
            mpi_max, icomm, ierrmpi)

        B1max = B1max_recv
        B2max = B2max_recv
        B3max = B3max_recv

    end subroutine get_max_B

end module mod_gimli
!