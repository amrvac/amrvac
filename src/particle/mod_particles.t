module mod_particles
  use mod_particle_base
  use mod_particle_advect
  use mod_particle_Lorentz
  use mod_particle_gca

  implicit none

contains

  subroutine particles_init()
    print *, "calling particle_base_init"
    call particle_base_init()
    call init_particles_com()

    select case(physics_type_particles)
    case('advect')
      call particle_advect_init()
    case('Lorentz')
      call particle_Lorentz_init()
    case('gca')
      call particle_gca_init()
    case default
      call mpistop("unknown physics_type_particles (advect,gca,Lorentz)")
    end select

  end subroutine particles_init

end module mod_particles
