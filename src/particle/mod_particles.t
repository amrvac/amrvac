!> Module containing all the particle routines
module mod_particles
  use mod_particle_base
  use mod_particle_advect
  use mod_particle_lorentz
  use mod_particle_gca
  use mod_particle_sample
  use mod_comm_lib, only: mpistop

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
    case('GCA')
      call gca_init()
    case('sample')
      call sample_init()
    case default
      if (mype == 0) then
        print *, "Unknown physics_type_particles", &
             trim(physics_type_particles)
        print *, "Options are: advect, Lorentz, GCA, sample"
        call mpistop("Unknown physics_type_particles")
      end if
    end select

  end subroutine particles_init

  !> Create initial particles
  subroutine particles_create()
    use mod_global_parameters

    select case(physics_type_particles)
    case('advect')
      call advect_create_particles()
    case('Lorentz')
      ! The Vay mover can use the same routine
      call Lorentz_create_particles()
    case('GCA')
      call gca_create_particles()
    case('sample')
      call sample_create_particles()
    case default
      if (mype == 0) then
        print *, "Unknown physics_type_particles", &
             trim(physics_type_particles)
        print *, "Options are: advect, Lorentz, GCA, sample"
        call mpistop("Unknown physics_type_particles")
      end if
    end select

  end subroutine particles_create

end module mod_particles
