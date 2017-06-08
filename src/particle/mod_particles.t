module mod_particles
  use mod_particle_base
  use mod_particle_advect
  use mod_particle_Lorentz
  use mod_particle_gca

  implicit none

contains

  !> Initialize particle data and parameters
  subroutine particles_init()
    call particle_base_init()
    call init_particles_com()

    select case(physics_type_particles)
    case('advect')
      call advect_init()
    case('Lorentz')
      call Lorentz_init()
    case('gca')
      call gca_init()
    case default
      call mpistop("unknown physics_type_particles (advect,gca,Lorentz)")
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
    case('Lorentz')
      call Lorentz_create_particles()
    case('gca')
      call gca_create_particles()
    case default
      call mpistop("unknown physics_type_particles (advect,gca,Lorentz)")
    end select

    ! Remove grid variables again
    call finish_gridvars()

  end subroutine particles_create

end module mod_particles
