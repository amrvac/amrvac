module mod_usr
  use mod_global_parameters
  use mod_physics
  use mod_geometry

  use mod_hd

  implicit none
  private

  public :: usr_init

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_initialize
    use mod_variables

    usr_set_bc        => rm_set_bc
    usr_init_one_grid => rm_init_one_grid

    call set_coordinate_system("Cartesian_1D")

    call hd_activate()

  end subroutine usr_init

  subroutine rm_set_bc()
     {
     integer, parameter :: i_min^D = 2*^D - 1
     integer, parameter :: i_max^D = 2*^D
     \}
     bc_type%m = bc_cont
     bc_type%w = bc_cont
  end subroutine rm_set_bc

  ! Initialize one grid
  subroutine rm_init_one_grid(ixI^L,ixO^L,s)
    use mod_global_parameters
    use mod_physics
    use mod_hd
    use mod_hd_parameters

    integer, intent(in)        :: ixI^L, ixO^L
    type(state), intent(inout) :: s

    associate( w => s%w,          &
               x => s%mesh%x,     &
               m => s%metric%vars &
             )

    ! initialise prim
    s%is_prim = .True.
    w = 0.0d0
    
     ! iprob==1 rarefaction wave & shock
    if (iprob==1) then
        where (abs(x(ixO^S,1))<0.5d0)
           w(ixO^S,rho_)   = 1.5d0
           w(ixO^S,mom(1)) = 0.0d0
           w(ixO^S,e_)     = 1.5d0
        elsewhere
           w(ixO^S,rho_)   = 0.5d0
           w(ixO^S,mom(1)) = 0.0d0
           w(ixO^S,e_)     = 0.5d0
        end where
    ! iprob==2  shock & shock
    else if (iprob==2) then
        where (abs(x(ixO^S,1))<0.5d0)
           w(ixO^S,rho_)   = 1.0d0
           w(ixO^S,mom(1)) = 2.0d0
           w(ixO^S,e_)     = 1.5d0
        elsewhere
           w(ixO^S,rho_)   = 1.0d0
           w(ixO^S,mom(1)) = -2.0d0
           w(ixO^S,e_)     = 1.5d0
        end where
    ! iprob==3  rarefaction wave
    else if (iprob==3) then
        where (abs(x(ixO^S,1))<0.5d0)
           w(ixO^S,rho_)   = 1.0d0
           w(ixO^S,mom(1)) = -0.5d0
           w(ixO^S,e_)     = 1.0d0
        elsewhere
           w(ixO^S,rho_)   = 1.0d0
           w(ixO^S,mom(1)) = 0.5d0
           w(ixO^S,e_)     = 1.0d0
        end where
    else
        call mpistop("iprob not available!")
    end if
    
    ! Before, e_ was set to the pressure -> use e = p/(gamma-1)
    w(ixO^S, e_) = w(ixO^S, e_) / (hd_gamma - 1)

    call phys_to_conserved(ixI^L,ixO^L,s)
    end associate
  end subroutine rm_init_one_grid

end module mod_usr
