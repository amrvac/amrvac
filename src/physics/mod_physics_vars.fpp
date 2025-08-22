#:mute
#:include 'mod_physics_templates.fpp'
#:endmute

module mod_physics_vars
  use mod_global_parameters, only: name_len, ndim

  implicit none
  public
  
  double precision :: phys_gamma=5.d0/3.d0
  !$acc declare copyin(phys_gamma)

  !> String describing the physics type of the simulation
  character(len=name_len) :: physics_type = "${PHYS}$"
  !$acc declare copyin(physics_type)

  !> To use wider stencils in flux calculations. A value of 1 will extend it by
  !> one cell in both directions, in any dimension
  integer :: phys_wider_stencil = 0
  !$acc declare copyin(phys_wider_stencil)

  !> Whether the physics routines require diagonal ghost cells, for example for
  !> computing a curl.
  logical :: phys_req_diagonal = .true.
  !$acc declare copyin(phys_req_diagonal)

  !> Solve energy equation or not
  logical :: phys_energy=.false.
  !$acc declare copyin(phys_energy)
  
  !> Solve total energy equation or not
  logical :: phys_total_energy=.false.
  !$acc declare copyin(phys_total_energy)

  !> Solve internal energy instead of total energy
  logical :: phys_internal_e=.false.
  !$acc declare copyin(phys_internal_e)

  !> Solve partially ionized one-fluid plasma
  logical :: phys_partial_ionization=.false.
  !$acc declare copyin(phys_partial_ionization)

  !> if equilibrium pressure is splitted
  logical :: phys_equi_pe=.false.
  !$acc declare copyin(phys_equi_pe)

  @:phys_vars()

end module mod_physics_vars
