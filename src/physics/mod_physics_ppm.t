module mod_physics_ppm

  implicit none
  public

  procedure(sub_ppm_flatcd), pointer      :: phys_ppm_flatcd => null()
  procedure(sub_ppm_flatsh), pointer      :: phys_ppm_flatsh => null()

  abstract interface
     subroutine sub_ppm_flatcd(ixI^L,ixO^L,ixL^L,ixR^L,w,d2w,drho,dp)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L,ixO^L,ixL^L,ixR^L
       double precision, intent(in)    :: w(ixI^S,nw),d2w(ixI^S,1:nwflux)
       double precision, intent(inout) :: drho(ixI^S),dp(ixI^S)
     end subroutine sub_ppm_flatcd

     subroutine sub_ppm_flatsh(ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L,idims,w,drho,dp,dv)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L
       integer, intent(in)             :: idims
       double precision, intent(in)    :: w(ixI^S,nw)
       double precision, intent(inout) :: drho(ixI^S),dp(ixI^S),dv(ixI^S)
     end subroutine sub_ppm_flatsh
  end interface

contains

  subroutine phys_ppm_check
    if (.not. associated(phys_ppm_flatcd)) &
         phys_ppm_flatcd => dummy_ppm_flatcd

    if (.not. associated(phys_ppm_flatsh)) &
         phys_ppm_flatsh => dummy_ppm_flatsh
  end subroutine phys_ppm_check

  subroutine dummy_ppm_flatcd(ixI^L,ixO^L,ixL^L,ixR^L,w,d2w,drho,dp)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L,ixL^L,ixR^L
    double precision, intent(in)    :: w(ixI^S,nw),d2w(ixI^S,1:nwflux)
    double precision, intent(inout) :: drho(ixI^S),dp(ixI^S)
    drho(ixO^S)=zero
    dp(ixO^S)=zero
  end subroutine dummy_ppm_flatcd

  subroutine dummy_ppm_flatsh(ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L,idims,w,drho,dp,dv)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(inout) :: drho(ixI^S),dp(ixI^S),dv(ixI^S)
    drho(ixO^S)=zero
    dp(ixO^S)=zero
    dv(ixO^S)=zero
  end subroutine dummy_ppm_flatsh

end module mod_physics_ppm
