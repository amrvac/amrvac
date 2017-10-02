!> Module containing all the particle routines
module mod_particles
  use mod_particle_base
  use mod_particle_advect
  use mod_particle_lorentz
  use mod_particle_gca
  use mod_particle_vay
  use mod_particle_lfimp
  use mod_particle_hc

  implicit none

contains

  !> Initialize particle data and parameters
  subroutine particles_init()
    use mod_global_parameters

    call particle_base_init()
    call init_particles_com()

    select case(physics_type_particles)
    case('advect')
      call advect_init()
    case('Lorentz')
      call Lorentz_init()
    case('gca')
      call gca_init()
    case('lfimp')
      call lfimp_init()
    case('Vay')
      call Vay_init()
    case('HC')
      call HC_init()
    case default
      if (mype == 0) then
        print *, "Unknown physics_type_particles", &
             trim(physics_type_particles)
        print *, "Options are: advect, gca, Lorentz, lfimp, Vay, HC"
        call mpistop("Unknown physics_type_particles")
      end if
    end select

  end subroutine particles_init

  !> Create initial particles
  subroutine particles_create()
    use mod_global_parameters
    use mod_timing

    ! Allocate grid variables
    tpartc_grid_0=MPI_WTIME()
    call init_gridvars()
    tpartc_grid = tpartc_grid + (MPI_WTIME()-tpartc_grid_0)

    select case(physics_type_particles)
    case('advect')
      call advect_create_particles()
    case('Lorentz', 'Vay', 'HC')
      ! The Vay mover can use the same routine
      call Lorentz_create_particles()
    case('gca')
      call gca_create_particles()
    case('lfimp')
      call lfimp_create_particles()
    case default
      if (mype == 0) then
        print *, "Unknown physics_type_particles", &
             trim(physics_type_particles)
        print *, "Options are: advect, gca, Lorentz, lfimp, Vay"
        call mpistop("Unknown physics_type_particles")
      end if
    end select

    ! Remove grid variables again
    call finish_gridvars()

  end subroutine particles_create

end module mod_particles
