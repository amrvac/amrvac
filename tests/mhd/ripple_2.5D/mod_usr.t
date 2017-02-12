module mod_usr
  use mod_mhd
  implicit none

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 

    call set_coordinate_system("Cartesian_2.5D")

    call mhd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
  ! initialize one grid
    use mod_global_parameters
    use mod_physics

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixO^S,mom(1))=0.d0
    w(ixO^S,mom(2))=0.d0
    w(ixO^S,mom(3))=0.01d0*dexp(-(^D&x(ixO^S,^D)**2+)/1d-4)
    w(ixO^S,p_)    =1.d0+0.1*dexp(-(^D&x(ixO^S,^D)**2+)/1d-4)
    w(ixO^S,rho_)  =1.d0+0.1*dexp(-(^D&x(ixO^S,^D)**2+)/1d-4)
    w(ixO^S,mag(1))=1.d0
    w(ixO^S,mag(2))=0.d0
    w(ixO^S,mag(3))=0.d0
    
    call phys_to_conserved(ixI^L,ixO^L,w,x)
  
  end subroutine initonegrid_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_global_parameters
    use mod_physics

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision                   :: tmp(ixI^S)

    call phys_get_pthermal(w,x,ixI^L,ixO^L,tmp)
    ! output the temperature
    w(ixO^S,nw+1)=tmp(ixO^S)/w(ixO^S,rho_)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='Te'

  end subroutine specialvarnames_output

end module mod_usr
