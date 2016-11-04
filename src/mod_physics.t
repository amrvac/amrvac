!> This module defines the procedures of a physics module. It contains function
!> pointers for the various supported routines. An actual physics module has to
!> set these pointers to its implementation of these routines.
module mod_physics

  implicit none

  character(len=40) :: physics_type

  procedure(sub_init_params), pointer     :: phys_init_params => null()
  procedure(sub_check_params), pointer    :: phys_check_params => null()
  procedure(sub_read_params), pointer     :: phys_read_params => null()
  procedure(sub_convert), pointer         :: phys_to_conserved => null()
  procedure(sub_convert), pointer         :: phys_to_primitive => null()
  procedure(sub_convert), pointer         :: phys_convert_before_prolong => null()
  procedure(sub_convert), pointer         :: phys_convert_after_prolong => null()
  procedure(sub_convert), pointer         :: phys_convert_before_coarsen => null()
  procedure(sub_convert), pointer         :: phys_convert_after_coarsen => null()
  procedure(sub_average), pointer         :: phys_average => null()
  procedure(sub_get_eigenjump), pointer   :: phys_get_eigenjump => null()
  procedure(sub_rtimes), pointer          :: phys_rtimes => null()
  procedure(sub_get_v), pointer           :: phys_get_v => null()
  procedure(sub_get_cmax), pointer        :: phys_get_cmax => null()
  procedure(sub_get_flux), pointer        :: phys_get_flux => null()
  procedure(sub_get_dt), pointer          :: phys_get_dt => null()
  procedure(sub_add_source_geom), pointer :: phys_add_source_geom => null()
  procedure(sub_add_source), pointer      :: phys_add_source => null()
  procedure(sub_get_aux), pointer         :: phys_get_aux => null()
  procedure(sub_check_w), pointer         :: phys_check_w => null()
  procedure(sub_ppm_flatcd), pointer      :: phys_ppm_flatcd => null()
  procedure(sub_ppm_flatsh), pointer      :: phys_ppm_flatsh => null()
  procedure(sub_diffuse_hllcd), pointer   :: phys_diffuse_hllcd => null()
  procedure(sub_get_lCD), pointer         :: phys_get_lCD => null()
  procedure(sub_get_wCD), pointer         :: phys_get_wCD => null()

  abstract interface

     subroutine sub_init_params()
     end subroutine sub_init_params

     subroutine sub_check_params()
     end subroutine sub_check_params

     subroutine sub_read_params(file_unit, success)
       integer, intent(in) :: file_unit
       logical, intent(out) :: success
     end subroutine sub_read_params

     subroutine sub_convert(ixI^L, ixO^L, w, x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(inout) :: w(ixI^S, nw)
       double precision, intent(in)    :: x(ixI^S, 1:^ND)
     end subroutine sub_convert

     subroutine sub_average(wL, wR, x, ix^L, idim, wroe, workroe)
       use mod_global_parameters
       integer, intent(in)                          :: ix^L, idim
       double precision, dimension(ixG^T, nw)       :: wL, wR, wroe
       double precision, dimension(ixG^T, nworkroe) :: workroe
       double precision, dimension(ixG^T, 1:^ND)   :: x
     end subroutine sub_average

     subroutine sub_get_eigenjump(wL, wR, wC, x, ix^L, il, &
          idim, smalla, a, jump, workroe)
       use mod_global_parameters
       integer                                      :: ix^L, jx^L, ixC^L, il, idim
       double precision, dimension(ixG^T, nw)       :: wL, wR, wC
       double precision, dimension(ixG^T)           :: smalla, a, jump, v
       double precision, dimension(ixG^T, nworkroe) :: workroe
       double precision, dimension(ixG^T, 1:^ND)   :: x
     end subroutine sub_get_eigenjump

     subroutine sub_rtimes(q, w, ix^L, iw, il, idim, rq, workroe)
       use mod_global_parameters
       integer, intent(in)             :: ix^L, iw, il, idim
       double precision, intent(in)    :: w(ixG^T, nw), q(ixG^T)
       double precision, intent(inout) :: rq(ixG^T)
       double precision, intent(inout) :: workroe(ixG^T, nworkroe)
     end subroutine sub_rtimes

     subroutine sub_get_v(w, x, ixI^L, ixO^L, idim, v)
       use mod_global_parameters
       integer, intent(in)           :: ixI^L, ixO^L, idim
       double precision, intent(in)  :: w(ixI^S, nw), x(ixI^S, 1:^ND)
       double precision, intent(out) :: v(ixG^T)
     end subroutine sub_get_v

     subroutine sub_get_cmax(w, x, ixI^L, ixO^L, idim, cmax, cmin)
       use mod_global_parameters
       integer, intent(in)                       :: ixI^L, ixO^L, idim
       double precision, intent(in)              :: w(ixI^S, nw), x(ixI^S, 1:^ND)
       double precision, intent(inout)           :: cmax(ixG^T)
       double precision, optional, intent(inout) :: cmin(ixG^T)
     end subroutine sub_get_cmax

     subroutine sub_get_flux(w, x, ixI^L, ixO^L, iw, idim, f, transport)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, iw, idim
       double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:^ND)
       double precision, intent(inout) :: f(ixG^T)
       logical, intent(out)            :: transport
     end subroutine sub_get_flux

     subroutine sub_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: qdt, x(ixI^S, 1:^ND)
       double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
     end subroutine sub_add_source_geom

     subroutine sub_add_source(qdt, ixI^L, ixO^L, iw^LIM, qtC, wCT, qt, w, x, qsourcesplit)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
       double precision, intent(in)    :: qdt, qtC, qt, x(ixI^S, 1:^ND)
       double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
       logical, intent(in)             :: qsourcesplit
     end subroutine sub_add_source

     subroutine sub_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
       double precision, intent(inout) :: w(ixI^S, 1:nw), dtnew
     end subroutine sub_get_dt

     subroutine sub_get_aux(clipping,w,x,ixI^L,ixO^L,subname)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: x(ixI^S,1:ndim)
       double precision, intent(inout) :: w(ixI^S,nw)
       logical, intent(in)             :: clipping
       character(len=*)                :: subname
     end subroutine sub_get_aux

     subroutine sub_check_w(checkprimitive,ixI^L,ixO^L,w,flag)
       use mod_global_parameters
       logical             :: checkprimitive
       integer, intent(in) :: ixI^L, ixO^L
       double precision    :: w(ixI^S,nw)
       logical             :: flag(ixG^T)
     end subroutine sub_check_w

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

     subroutine sub_diffuse_hllcd(ixI^L,ixO^L,idims,wLC,wRC,fLC,fRC,patchf)
       use mod_global_parameters
       integer, intent(in)                                    :: ixI^L,ixO^L,idims
       double precision, dimension(ixG^T,1:nw), intent(in)    :: wRC,wLC
       double precision, dimension(ixG^T,1:nwflux),intent(in) :: fLC, fRC
       integer, dimension(ixG^T), intent(inout)               :: patchf
     end subroutine sub_diffuse_hllcd

     subroutine sub_get_lCD(wLC,wRC,fLC,fRC,cmin,cmax,idims,ixI^L,ixO^L, &
          whll,Fhll,lambdaCD,patchf)
       use mod_global_parameters
       integer, intent(in)                                      :: ixI^L,ixO^L,idims
       double precision, dimension(ixG^T,1:nw), intent(in)      :: wLC,wRC
       double precision, dimension(ixG^T,1:nwflux), intent(in)  :: fLC,fRC
       double precision, dimension(ixG^T), intent(in)           :: cmax,cmin
       integer, dimension(ixG^T), intent(inout)                 :: patchf
       double precision, dimension(ixG^T,1:nwflux), intent(out) :: Fhll,whll
       double precision, dimension(ixG^T), intent(out)          :: lambdaCD
     end subroutine sub_get_lCD

     subroutine sub_get_wCD(wLC,wRC,whll,vLC,vRC,fRC,fLC,Fhll,patchf,lambdaCD,cmin,cmax,&
          ixI^L,ixO^L,idims,f)
       use mod_global_parameters
       integer, intent(in)                                     :: ixI^L,ixO^L,idims
       double precision, dimension(ixG^T,1:nw), intent(in)     :: wRC,wLC
       double precision, dimension(ixG^T,1:nwflux), intent(in) :: whll, Fhll
       double precision, dimension(ixG^T), intent(in)          :: vRC, vLC,lambdaCD
       double precision, dimension(ixG^T), intent(in)          :: cmax,cmin
       double precision, dimension(ixG^T,1:nwflux), intent(in) :: fRC,fLC
       integer, dimension(ixG^T), intent(in)                   :: patchf
       double precision, dimension(ixG^T,1:nwflux),intent(out) :: f
     end subroutine sub_get_wCD

  end interface

contains

  subroutine phys_check_methods()
    ! Checks whether the required physics methods have been defined

    if (.not. associated(phys_init_params)) &
         phys_init_params => dummy_init_params

    if (.not. associated(phys_check_params)) &
         phys_check_params => dummy_check_params

    if (.not. associated(phys_read_params)) &
         call mpistop("Error: no read_params method has been specified")

    if (.not. associated(phys_to_conserved)) &
         phys_to_conserved => dummy_convert

    if (.not. associated(phys_to_primitive)) &
         phys_to_primitive => dummy_convert

    if (.not. associated(phys_convert_before_prolong)) &
         phys_convert_before_prolong => dummy_convert

    if (.not. associated(phys_convert_after_prolong)) &
         phys_convert_after_prolong => dummy_convert

    if (.not. associated(phys_convert_before_coarsen)) &
         phys_convert_before_coarsen => dummy_convert

    if (.not. associated(phys_convert_after_coarsen)) &
         phys_convert_after_coarsen => dummy_convert

    if (.not. associated(phys_average)) &
         call mpistop("Error: no average method has been specified")

    if (.not. associated(phys_get_eigenjump)) &
         call mpistop("Error: no eigenjump method has been specified")

    if (.not. associated(phys_rtimes)) &
         call mpistop("Error: no rtimes method has been specified")

    if (.not. associated(phys_get_v)) &
         call mpistop("Error: no get_v method has been specified")

    if (.not. associated(phys_get_cmax)) &
         call mpistop("Error: no get_cmax method has been specified")

    if (.not. associated(phys_get_flux)) &
         call mpistop("Error: no get_flux method has been specified")

    if (.not. associated(phys_get_dt)) &
         phys_get_dt => dummy_get_dt

    if (.not. associated(phys_add_source_geom)) &
         phys_add_source_geom => dummy_add_source_geom

    if (.not. associated(phys_add_source)) &
         phys_add_source => dummy_add_source

    if (.not. associated(phys_get_aux)) &
         phys_get_aux => dummy_get_aux

    if (.not. associated(phys_check_w)) &
         phys_check_w => dummy_check_w

    if (.not. associated(phys_ppm_flatcd)) &
         phys_ppm_flatcd => dummy_ppm_flatcd

    if (.not. associated(phys_ppm_flatsh)) &
         phys_ppm_flatsh => dummy_ppm_flatsh

    if (.not. associated(phys_diffuse_hllcd)) &
         phys_diffuse_hllcd => dummy_diffuse_hllcd

    if (.not. associated(phys_get_lCD)) &
         phys_get_lCD => dummy_get_lCD

    if (.not. associated(phys_get_wCD)) &
         phys_get_wCD => dummy_get_wCD

  end subroutine phys_check_methods

  subroutine dummy_init_params()
  end subroutine dummy_init_params

  subroutine dummy_check_params()
  end subroutine dummy_check_params

  subroutine dummy_convert(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)
  end subroutine dummy_convert

  subroutine dummy_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S, 1:^ND)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
  end subroutine dummy_add_source_geom

  subroutine dummy_add_source(qdt, ixI^L, ixO^L, iw^LIM, qtC, wCT, qt, w, x, qsourcesplit)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt, x(ixI^S, 1:^ND)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
    logical, intent(in)             :: qsourcesplit
  end subroutine dummy_add_source

  subroutine dummy_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
    double precision, intent(inout) :: w(ixI^S, 1:nw), dtnew

    dtnew = bigdouble
  end subroutine dummy_get_dt

  subroutine dummy_get_aux(clipping,w,x,ixI^L,ixO^L,subname)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,nw)
    logical, intent(in)             :: clipping
    character(len=*)                :: subname
  end subroutine dummy_get_aux

  subroutine dummy_check_w(checkprimitive,ixI^L,ixO^L,w,flag)
    use mod_global_parameters
    logical             :: checkprimitive
    integer, intent(in) :: ixI^L, ixO^L
    double precision    :: w(ixI^S,nw)
    logical             :: flag(ixG^T)
    flag(ixO^S)=.true.
  end subroutine dummy_check_w

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

  ! When method is hllcd or hllcd1 then: this subroutine is to impose enforce
  ! regions where we AVOID HLLC and use TVDLF instead: this is achieved by setting
  ! patchf to 4 in certain regions. An additional input parameter is nxdiffusehllc
  ! which sets the size of the fallback region. This nul version enforces TVDLF
  ! everywhere!
  subroutine dummy_diffuse_hllcd(ixI^L,ixO^L,idims,wLC,wRC,fLC,fRC,patchf)
    use mod_global_parameters
    integer, intent(in)                                    :: ixI^L,ixO^L,idims
    double precision, dimension(ixG^T,1:nw), intent(in)    :: wRC,wLC
    double precision, dimension(ixG^T,1:nwflux),intent(in) :: fLC, fRC
    integer, dimension(ixG^T), intent(inout)               :: patchf

    patchf(ixO^S) = 4 ! enforce TVDLF everywhere
  end subroutine dummy_diffuse_hllcd

  ! Calculate lambda at CD and set the patchf to know the orientation
  ! of the riemann fan and decide on the flux choice
  ! We also compute here the HLL flux and w value, for fallback strategy
  ! In this nul version, we simply compute nothing and ensure TVDLF fallback
  subroutine dummy_get_lCD(wLC,wRC,fLC,fRC,cmin,cmax,idims,ixI^L,ixO^L, &
       whll,Fhll,lambdaCD,patchf)
    use mod_global_parameters
    integer, intent(in)                                      :: ixI^L,ixO^L,idims
    double precision, dimension(ixG^T,1:nw), intent(in)      :: wLC,wRC
    double precision, dimension(ixG^T,1:nwflux), intent(in)  :: fLC,fRC
    double precision, dimension(ixG^T), intent(in)           :: cmax,cmin
    integer, dimension(ixG^T), intent(inout)                 :: patchf
    double precision, dimension(ixG^T,1:nwflux), intent(out) :: Fhll,whll
    double precision, dimension(ixG^T), intent(out)          :: lambdaCD

    ! Next must normally be computed
    Fhll(ixO^S,1:nwflux) = zero
    whll(ixO^S,1:nwflux) = zero
    lambdaCD(ixO^S)      = zero

    ! This actually ensures fallback to TVDLF
    patchf(ixO^S)=4
  end subroutine dummy_get_lCD

  ! compute the intermediate state U*. Only needed where patchf=-1/1. This nul
  ! version simply nullifies all values
  subroutine dummy_get_wCD(wLC,wRC,whll,vLC,vRC,fRC,fLC,Fhll,patchf,lambdaCD,cmin,cmax,&
       ixI^L,ixO^L,idims,f)
    use mod_global_parameters
    integer, intent(in)                                     :: ixI^L,ixO^L,idims
    double precision, dimension(ixG^T,1:nw), intent(in)     :: wRC,wLC
    double precision, dimension(ixG^T,1:nwflux), intent(in) :: whll, Fhll
    double precision, dimension(ixG^T), intent(in)          :: vRC, vLC,lambdaCD
    double precision, dimension(ixG^T), intent(in)          :: cmax,cmin
    double precision, dimension(ixG^T,1:nwflux), intent(in) :: fRC,fLC
    integer, dimension(ixG^T), intent(in)                   :: patchf
    double precision, dimension(ixG^T,1:nwflux),intent(out) :: f

    ! Next must normally be computed
    f(ixO^S,1:nwflux)  = zero
  end subroutine dummy_get_wCD

end module mod_physics
