!>   The module add viscous source terms and check time step
!>
!>   Viscous forces in the momentum equations:
!>   d m_i/dt +=  div (vc_mu * PI)
!>   !! Viscous work in the energy equation:
!>   !! de/dt +=  div (v . vc_mu * PI)
!>   where the PI stress tensor is
!>   PI_i,j = (dv_j/dx_i + dv_i/dx_j) - Kronecker delta_i,j*(2/3)*Sum_k dv_k/dx_k
!>   where vc_mu is the dynamic viscosity coefficient (g cm^-1 s^-1).
module mod_viscosity
  use mod_comm_lib, only: mpistop
  use mod_physics, only: phys_get_rho, phys_internal_e
  implicit none

  !> Viscosity coefficient
  double precision, public :: vc_mu = 1.d0

  !> Indices of the momentum density
  integer, allocatable, private, protected :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, private, protected              :: e_
  !> Indices of the velocity for the form of better vectorization
  integer, public, protected              :: v1_,v2_,v3_

  !> source split or not
  logical :: vc_split= .false.

  procedure(sub_add_source), pointer :: viscosity_add_source => null()
  ! Public methods
  public :: viscosity_add_source
  public :: viscosity_get_dt
  public :: viscosity_init

contains
  !> Read module parameters from a file
  subroutine vc_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /vc_list/ vc_mu, vc_split

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, vc_list, end=111)
111    close(unitpar)
    end do


  end subroutine vc_params_read

  !> Initialize the module
  subroutine viscosity_init(phys_wider_stencil)
    use mod_global_parameters
    use mod_geometry
    integer, intent(inout) :: phys_wider_stencil
    integer :: nwx,idir

    call vc_params_read(par_files)

    if(vc_split) any_source_split=.true.

    ! Determine flux variables
    nwx = 1                  ! rho (density)

    allocate(mom(ndir))
    do idir = 1, ndir
       nwx    = nwx + 1
       mom(idir) = nwx       ! momentum density
    end do
    v^C_=mom(^C);

    nwx = nwx + 1
    e_     = nwx          ! energy density

    select case (coordinate)
      case (Cartesian,Cartesian_stretched,Cartesian_expansion)
        viscosity_add_source => viscosity_add_source_Cartesian
      case (cylindrical)
        viscosity_add_source => viscosity_add_source_cylinder
      case (spherical)
        viscosity_add_source => viscosity_add_source_sphere
    end select

  end subroutine viscosity_init

  subroutine viscosity_add_source_Cartesian(qdt,ixI^L,ixO^L,wCT,wp,w,x,&
       energy,qsourcesplit,active)
  ! Add viscosity source in isotropic Newtonian fluids to w within ixO
  ! neglecting bulk viscosity
  ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods, only: usr_set_viscosity

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw), wp(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active

    double precision:: lambda(ixI^S,ndir,ndir),tmp(ixI^S),vlambda(ixI^S,ndir),divv23,nabla_v(ndim,ndir),vc_mus(ixI^S)
    integer:: ix^L,ix^D
    logical :: total_energy=.false.

    if(qsourcesplit .eqv. vc_split) then
      active = .true.
      if(energy.and..not.phys_internal_e) total_energy=.true.
      ! standard case, textbook viscosity
      ! Calculating viscosity sources
      ! involves second derivatives, two extra layers
      ix^L=ixO^L^LADD2;
      if({ ixImin^D>ixmin^D .or. ixImax^D<ixmax^D|.or.})&
        call mpistop("error for viscous source addition, 2 layers needed")
      ix^L=ixO^L^LADD1;
      if(associated(usr_set_viscosity)) then
        call usr_set_viscosity(ixI^L,ix^L,wp,x,vc_mus)
      else
        vc_mus=vc_mu
      end if
      {^IFTHREED
      {do ix^DB=ixmin^DB,ixmax^DB\}
        ! idim=1, idir=1
        lambda(ix^D,1,1)=(wp(ix1+1,ix2,ix3,v1_)-wp(ix1-1,ix2,ix3,v1_))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))
        ! idim=1, idir=2
        lambda(ix^D,1,2)=(wp(ix1+1,ix2,ix3,v2_)-wp(ix1-1,ix2,ix3,v2_))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))
        ! idim=1, idir=3
        lambda(ix^D,1,3)=(wp(ix1+1,ix2,ix3,v3_)-wp(ix1-1,ix2,ix3,v3_))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))
        ! idim=2, idir=1
        lambda(ix^D,2,1)=(wp(ix1,ix2+1,ix3,v1_)-wp(ix1,ix2-1,ix3,v1_))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))
        ! idim=2, idir=2
        lambda(ix^D,2,2)=(wp(ix1,ix2+1,ix3,v2_)-wp(ix1,ix2-1,ix3,v2_))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))
        ! idim=2, idir=3
        lambda(ix^D,2,3)=(wp(ix1,ix2+1,ix3,v3_)-wp(ix1,ix2-1,ix3,v3_))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))
        ! idim=3, idir=1
        lambda(ix^D,3,1)=(wp(ix1,ix2,ix3+1,v1_)-wp(ix1,ix2,ix3-1,v1_))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))
        ! idim=3, idir=2
        lambda(ix^D,3,2)=(wp(ix1,ix2,ix3+1,v2_)-wp(ix1,ix2,ix3-1,v2_))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))
        ! idim=3, idir=3
        lambda(ix^D,3,3)=(wp(ix1,ix2,ix3+1,v3_)-wp(ix1,ix2,ix3-1,v3_))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))
        if(phys_internal_e) then
          nabla_v(1,1)=lambda(ix^D,1,1)
          nabla_v(1,2)=lambda(ix^D,1,2)
          nabla_v(1,3)=lambda(ix^D,1,3)
          nabla_v(2,1)=lambda(ix^D,2,1)
          nabla_v(2,2)=lambda(ix^D,2,2)
          nabla_v(2,3)=lambda(ix^D,2,3)
          nabla_v(3,1)=lambda(ix^D,3,1)
          nabla_v(3,2)=lambda(ix^D,3,2)
          nabla_v(3,3)=lambda(ix^D,3,3)
        end if
        ! dv_i/d_j + dv_j/d_i
        lambda(ix^D,1,2)=(lambda(ix^D,1,2)+lambda(ix^D,2,1))*vc_mus(ix^D)
        lambda(ix^D,1,3)=(lambda(ix^D,1,3)+lambda(ix^D,3,1))*vc_mus(ix^D)
        lambda(ix^D,2,3)=(lambda(ix^D,2,3)+lambda(ix^D,3,2))*vc_mus(ix^D)
        lambda(ix^D,2,1)=lambda(ix^D,1,2)
        lambda(ix^D,3,1)=lambda(ix^D,1,3)
        lambda(ix^D,3,2)=lambda(ix^D,2,3)
        divv23=two*third*(lambda(ix^D,1,1)+lambda(ix^D,2,2)+lambda(ix^D,3,3))
        lambda(ix^D,1,1)=(two*lambda(ix^D,1,1)-divv23)*vc_mus(ix^D)
        lambda(ix^D,2,2)=(two*lambda(ix^D,2,2)-divv23)*vc_mus(ix^D)
        lambda(ix^D,3,3)=(two*lambda(ix^D,3,3)-divv23)*vc_mus(ix^D)
        if(total_energy) then
          vlambda(ix^D,1)=wp(ix^D,v1_)*lambda(ix^D,1,1)+wp(ix^D,v2_)*lambda(ix^D,2,1)+wp(ix^D,v3_)*lambda(ix^D,3,1)
          vlambda(ix^D,2)=wp(ix^D,v1_)*lambda(ix^D,1,2)+wp(ix^D,v2_)*lambda(ix^D,2,2)+wp(ix^D,v3_)*lambda(ix^D,3,2)
          vlambda(ix^D,3)=wp(ix^D,v1_)*lambda(ix^D,1,3)+wp(ix^D,v2_)*lambda(ix^D,2,3)+wp(ix^D,v3_)*lambda(ix^D,3,3)
        end if
        if(phys_internal_e) then
          w(ix^D,e_)=w(ix^D,e_)+&
             qdt*(lambda(ix^D,1,1)*nabla_v(1,1)+lambda(ix^D,2,1)*nabla_v(2,1)+lambda(ix^D,3,1)*nabla_v(3,1)&
                 +lambda(ix^D,1,2)*nabla_v(1,2)+lambda(ix^D,2,2)*nabla_v(2,2)+lambda(ix^D,3,2)*nabla_v(3,2)&
                 +lambda(ix^D,1,3)*nabla_v(1,3)+lambda(ix^D,2,3)*nabla_v(2,3)+lambda(ix^D,3,3)*nabla_v(3,3))
        end if
      {end do\}
      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        w(ix^D,mom(1))=((lambda(ix1+1,ix2,ix3,1,1)-lambda(ix1-1,ix2,ix3,1,1))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))+&
                        (lambda(ix1,ix2+1,ix3,2,1)-lambda(ix1,ix2-1,ix3,2,1))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))+&
                        (lambda(ix1,ix2,ix3+1,3,1)-lambda(ix1,ix2,ix3-1,3,1))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3)))*qdt+w(ix^D,mom(1))
        w(ix^D,mom(2))=((lambda(ix1+1,ix2,ix3,1,2)-lambda(ix1-1,ix2,ix3,1,2))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))+&
                        (lambda(ix1,ix2+1,ix3,2,2)-lambda(ix1,ix2-1,ix3,2,2))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))+&
                        (lambda(ix1,ix2,ix3+1,3,2)-lambda(ix1,ix2,ix3-1,3,2))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3)))*qdt+w(ix^D,mom(2))
        w(ix^D,mom(3))=((lambda(ix1+1,ix2,ix3,1,3)-lambda(ix1-1,ix2,ix3,1,3))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))+&
                        (lambda(ix1,ix2+1,ix3,2,3)-lambda(ix1,ix2-1,ix3,2,3))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))+&
                        (lambda(ix1,ix2,ix3+1,3,3)-lambda(ix1,ix2,ix3-1,3,3))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3)))*qdt+w(ix^D,mom(3))
      {end do\}
      }
      {^IFTWOD
      {do ix^DB=ixmin^DB,ixmax^DB\}
        ! idim=1, idir=1
        lambda(ix^D,1,1)=(wp(ix1+1,ix2,v1_)-wp(ix1-1,ix2,v1_))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))
        ! idim=1, idir=2
        lambda(ix^D,1,2)=(wp(ix1+1,ix2,v2_)-wp(ix1-1,ix2,v2_))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))
        ! idim=2, idir=1
        lambda(ix^D,2,1)=(wp(ix1,ix2+1,v1_)-wp(ix1,ix2-1,v1_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))
        ! idim=2, idir=2
        lambda(ix^D,2,2)=(wp(ix1,ix2+1,v2_)-wp(ix1,ix2-1,v2_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))
        if(phys_internal_e) then
          nabla_v(1,1)=lambda(ix^D,1,1)
          nabla_v(1,2)=lambda(ix^D,1,2)
          nabla_v(2,1)=lambda(ix^D,2,1)
          nabla_v(2,2)=lambda(ix^D,2,2)
          if(ndir==3) then
            nabla_v(1,3)=lambda(ix^D,1,3)
            nabla_v(2,3)=lambda(ix^D,2,3)
          end if
        end if
        ! dv_i/d_j + dv_j/d_i
        lambda(ix^D,1,2)=(lambda(ix^D,1,2)+lambda(ix^D,2,1))*vc_mus(ix^D)
        lambda(ix^D,2,1)=lambda(ix^D,1,2)
        divv23=two*third*(lambda(ix^D,1,1)+lambda(ix^D,2,2))
        lambda(ix^D,1,1)=(two*lambda(ix^D,1,1)-divv23)*vc_mus(ix^D)
        lambda(ix^D,2,2)=(two*lambda(ix^D,2,2)-divv23)*vc_mus(ix^D)
        if(ndir==3) then
          ! idim=1, idir=3
          lambda(ix1,ix2,1,3)=(wp(ix1+1,ix2,v3_)-wp(ix1-1,ix2,v3_))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))*vc_mus(ix^D)
          ! idim=2, idir=3
          lambda(ix1,ix2,2,3)=(wp(ix1,ix2+1,v3_)-wp(ix1,ix2-1,v3_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*vc_mus(ix^D)
          lambda(ix1,ix2,3,1)=lambda(ix1,ix2,1,3)
          lambda(ix1,ix2,3,2)=lambda(ix1,ix2,2,3)
          lambda(ix1,ix2,3,3)=-divv23*vc_mus(ix^D)
          if(total_energy) then
            vlambda(ix1,ix2,1)=wp(ix1,ix2,v1_)*lambda(ix1,ix2,1,1)+wp(ix1,ix2,v2_)*lambda(ix1,ix2,2,1)+wp(ix1,ix2,v3_)*lambda(ix1,ix2,3,1)
            vlambda(ix1,ix2,2)=wp(ix1,ix2,v1_)*lambda(ix1,ix2,1,2)+wp(ix1,ix2,v2_)*lambda(ix1,ix2,2,2)+wp(ix1,ix2,v3_)*lambda(ix1,ix2,3,2)
            vlambda(ix1,ix2,3)=wp(ix1,ix2,v1_)*lambda(ix1,ix2,1,3)+wp(ix1,ix2,v2_)*lambda(ix1,ix2,2,3)+wp(ix1,ix2,v3_)*lambda(ix1,ix2,3,3)
          end if
          if(phys_internal_e) then
            w(ix^D,e_)=w(ix^D,e_)+&
               qdt*(lambda(ix^D,1,1)*nabla_v(1,1)+lambda(ix^D,2,1)*nabla_v(2,1)&
                   +lambda(ix^D,1,2)*nabla_v(1,2)+lambda(ix^D,2,2)*nabla_v(2,2)&
                   +lambda(ix^D,1,3)*nabla_v(1,3)+lambda(ix^D,2,3)*nabla_v(2,3))
          end if
        else 
          if(total_energy) then
            vlambda(ix1,ix2,1)=wp(ix1,ix2,v1_)*lambda(ix1,ix2,1,1)+wp(ix1,ix2,v2_)*lambda(ix1,ix2,2,1)
            vlambda(ix1,ix2,2)=wp(ix1,ix2,v1_)*lambda(ix1,ix2,1,2)+wp(ix1,ix2,v2_)*lambda(ix1,ix2,2,2)
          end if
          if(phys_internal_e) then
            w(ix^D,e_)=w(ix^D,e_)+&
               qdt*(lambda(ix^D,1,1)*nabla_v(1,1)+lambda(ix^D,2,1)*nabla_v(2,1)&
                   +lambda(ix^D,1,2)*nabla_v(1,2)+lambda(ix^D,2,2)*nabla_v(2,2))
          end if
        end if
      {end do\}
      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        w(ix^D,mom(1))=((lambda(ix1+1,ix2,1,1)-lambda(ix1-1,ix2,1,1))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                        (lambda(ix1,ix2+1,2,1)-lambda(ix1,ix2-1,2,1))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2)))*qdt+w(ix^D,mom(1))
        w(ix^D,mom(2))=((lambda(ix1+1,ix2,1,2)-lambda(ix1-1,ix2,1,2))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                        (lambda(ix1,ix2+1,2,2)-lambda(ix1,ix2-1,2,2))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2)))*qdt+w(ix^D,mom(2))
        if(ndir==3) then
          w(ix1,ix2,mom(3))=((lambda(ix1+1,ix2,1,3)-lambda(ix1-1,ix2,1,3))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                             (lambda(ix1,ix2+1,2,3)-lambda(ix1,ix2-1,2,3))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2)))*qdt+w(ix1,ix2,mom(3))
        end if
      {end do\}
      }
      {^IFONED
      do ix1=ixmin1,ixmax1
        ! idim=1, idir=1
        lambda(ix1,1,1)=(wp(ix1+1,v1_)-wp(ix1-1,v1_))/(x(ix1+1,1)-x(ix1-1,1))
        ! dv_i/d_j + dv_j/d_i
        if(phys_internal_e) then
          nabla_v(1,1)=lambda(ix^D,1,1)
        end if
        divv23=two*third*lambda(ix1,1,1)
        lambda(ix1,1,1)=(two*lambda(ix1,1,1)-divv23)*vc_mus(ix^D)
        if(ndir==2) then
          ! idim=1, idir=2
          lambda(ix1,1,2)=(wp(ix1+1,v2_)-wp(ix1-1,v2_))/(x(ix1+1,1)-x(ix1-1,1))
          if(phys_internal_e)  nabla_v(1,2)=lambda(ix^D,1,2)
          lambda(ix1,1,2)=lambda(ix1,1,2)*vc_mus(ix^D)
          ! dv_i/d_j + dv_j/d_i
          lambda(ix1,2,1)=lambda(ix1,1,2)
          lambda(ix1,2,2)=-divv23*vc_mus(ix^D)
          if(total_energy) then
            vlambda(ix1,1)=wp(ix1,v1_)*lambda(ix1,1,1)+wp(ix1,v2_)*lambda(ix1,2,1)
            vlambda(ix1,2)=wp(ix1,v1_)*lambda(ix1,1,2)+wp(ix1,v2_)*lambda(ix1,2,2)
          end if
          if(phys_internal_e) then
            w(ix^D,e_)=w(ix^D,e_)+&
               qdt*(lambda(ix^D,1,1)*nabla_v(1,1)&
                   +lambda(ix^D,1,2)*nabla_v(1,2))
          end if
        else if(ndir==3) then
          ! idim=1, idir=2
          lambda(ix1,1,2)=(wp(ix1+1,v2_)-wp(ix1-1,v2_))/(x(ix1+1,1)-x(ix1-1,1))
          ! idim=1, idir=3
          lambda(ix1,1,3)=(wp(ix1+1,v3_)-wp(ix1-1,v3_))/(x(ix1+1,1)-x(ix1-1,1))
          if(phys_internal_e) then
            nabla_v(1,2)=lambda(ix^D,1,2)
            nabla_v(1,3)=lambda(ix^D,1,3)
          end if
          lambda(ix1,1,2)=lambda(ix1,1,2)*vc_mus(ix^D)
          lambda(ix1,1,3)=lambda(ix1,1,3)*vc_mus(ix^D)
          ! dv_i/d_j + dv_j/d_i
          lambda(ix1,2,1)=lambda(ix1,1,2)
          lambda(ix1,3,1)=lambda(ix1,1,3)
          lambda(ix1,2,2)=-divv23*vc_mus(ix^D)
          lambda(ix1,3,3)=-divv23*vc_mus(ix^D)
          if(total_energy) then
            vlambda(ix1,1)=wp(ix1,v1_)*lambda(ix1,1,1)+wp(ix1,v2_)*lambda(ix1,2,1)+wp(ix1,v3_)*lambda(ix1,3,1)
            vlambda(ix1,2)=wp(ix1,v1_)*lambda(ix1,1,2)+wp(ix1,v2_)*lambda(ix1,2,2)+wp(ix1,v3_)*lambda(ix1,3,2)
            vlambda(ix1,3)=wp(ix1,v1_)*lambda(ix1,1,3)+wp(ix1,v2_)*lambda(ix1,2,3)+wp(ix1,v3_)*lambda(ix1,3,3)
          end if
          if(phys_internal_e) then
            w(ix^D,e_)=w(ix^D,e_)+&
               qdt*(lambda(ix^D,1,1)*nabla_v(1,1)&
                   +lambda(ix^D,1,2)*nabla_v(1,2)&
                   +lambda(ix^D,1,3)*nabla_v(1,3))
          end if
        else 
          if(total_energy) then
            vlambda(ix1,1)=wp(ix1,v1_)*lambda(ix1,1,1)
          end if
          if(phys_internal_e) then
            w(ix^D,e_)=w(ix^D,e_)+&
               qdt*(lambda(ix^D,1,1)*nabla_v(1,1))
          end if
        end if
      end do
      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
      do ix1=ixOmin1,ixOmax1
        if(ndir==1) then
          w(ix1,mom(1))=(lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))*qdt+w(ix1,mom(1))
        else if(ndir==2) then
          w(ix1,mom(1))=(lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))*qdt+w(ix1,mom(1))
          w(ix1,mom(2))=(lambda(ix1+1,1,2)-lambda(ix1-1,1,2))/(x(ix1+1,1)-x(ix1-1,1))*qdt+w(ix1,mom(2))
        else
          w(ix1,mom(1))=(lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))*qdt+w(ix1,mom(1))
          w(ix1,mom(2))=(lambda(ix1+1,1,2)-lambda(ix1-1,1,2))/(x(ix1+1,1)-x(ix1-1,1))*qdt+w(ix1,mom(2))
          w(ix1,mom(3))=(lambda(ix1+1,1,3)-lambda(ix1-1,1,3))/(x(ix1+1,1)-x(ix1-1,1))*qdt+w(ix1,mom(3))
        end if
      end do
      }
      if(total_energy) then
        ! de/dt= +div(v.dot.[mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v *kr])
        ! thus e=e+d_i v_j tensor_ji
        call divvector(vlambda,ixI^L,ixO^L,tmp)
        w(ixO^S,e_)=w(ixO^S,e_)+tmp(ixO^S)*qdt
      end if
    end if

  end subroutine viscosity_add_source_Cartesian

  subroutine viscosity_add_source_sphere(qdt,ixI^L,ixO^L,wCT,wp,w,x,&
       energy,qsourcesplit,active)
  ! Add viscosity source in isotropic Newtonian fluids to w within ixO
  ! neglecting bulk viscosity
  ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods, only: usr_set_viscosity

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw), wp(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active

    double precision :: lambda(ixI^S,ndir,ndir),vlambda(ixI^S,ndir),invr(ixI^S),tctan(ixI^S),invrsin(ixI^S),nabla_v(ndir,ndir),vc_mus(ixI^S)
    double precision :: tsin,tcos,divv23
    integer:: ix^L,ix^D
    logical :: total_energy=.false.

    if(qsourcesplit .eqv. vc_split) then
      active = .true.
      if(energy.and..not.phys_internal_e) total_energy=.true.
      ! standard case, textbook viscosity
      ! Calculating viscosity sources
      ! involves second derivatives, two extra layers
      ix^L=ixO^L^LADD2;
      if({ ixImin^D>ixmin^D .or. ixImax^D<ixmax^D|.or.})&
        call mpistop("error for viscous source addition, 2 layers needed")
      ix^L=ixO^L^LADD1;

      ! construct lambda tensor: lambda_ij = gradv_ij + gradv_ji
      ! initialize
      if(associated(usr_set_viscosity)) then
        call usr_set_viscosity(ixI^L,ix^L,wp,x,vc_mus)
      else
        vc_mus=vc_mu
      end if
      {^IFTHREED
      {do ix^DB=ixmin^DB,ixmax^DB\}
        invr(ix^D)=1.d0/x(ix^D,1)
        tcos=dcos(x(ix^D,2))
        ! idim=1, idir=1
        lambda(ix^D,1,1)=(wp(ix1+1,ix2,ix3,v1_)-wp(ix1-1,ix2,ix3,v1_))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))
        ! idim=1, idir=2
        lambda(ix^D,1,2)=(wp(ix1+1,ix2,ix3,v2_)-wp(ix1-1,ix2,ix3,v2_))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))
        ! idim=1, idir=3
        lambda(ix^D,1,3)=(wp(ix1+1,ix2,ix3,v3_)-wp(ix1-1,ix2,ix3,v3_))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))
        ! idim=2, idir=1
        lambda(ix^D,2,1)=((wp(ix1,ix2+1,ix3,v1_)-wp(ix1,ix2-1,ix3,v1_))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))-wp(ix^D,v2_))*invr(ix^D)
        ! idim=2, idir=2
        lambda(ix^D,2,2)=((wp(ix1,ix2+1,ix3,v2_)-wp(ix1,ix2-1,ix3,v2_))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))+wp(ix^D,v1_))*invr(ix^D)
        ! idim=2, idir=3
        lambda(ix^D,2,3)=((wp(ix1,ix2+1,ix3,v3_)-wp(ix1,ix2-1,ix3,v3_))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))+wp(ix^D,v2_)*tcos)*invr(ix^D)
        tsin=dsin(x(ix^D,2))
        invrsin(ix^D)=invr(ix^D)/tsin
        tctan(ix^D)=tcos/tsin
        ! idim=3, idir=1
        lambda(ix^D,3,1)=((wp(ix1,ix2,ix3+1,v1_)-wp(ix1,ix2,ix3-1,v1_))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))-wp(ix^D,v3_)*tsin)*invrsin(ix^D)
        ! idim=3, idir=2
        lambda(ix^D,3,2)=((wp(ix1,ix2,ix3+1,v2_)-wp(ix1,ix2,ix3-1,v2_))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))-wp(ix^D,v3_)*tcos)*invrsin(ix^D)
        ! idim=3, idir=3
        lambda(ix^D,3,3)=((wp(ix1,ix2,ix3+1,v3_)-wp(ix1,ix2,ix3-1,v3_))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))+wp(ix^D,v1_)*tsin+&
            wp(ix^D,v2_)*tcos)*invrsin(ix^D)
        if(phys_internal_e) then
          nabla_v(1,1)=lambda(ix^D,1,1)
          nabla_v(1,2)=lambda(ix^D,1,2)
          nabla_v(1,3)=lambda(ix^D,1,3)
          nabla_v(2,1)=lambda(ix^D,2,1)
          nabla_v(2,2)=lambda(ix^D,2,2)
          nabla_v(2,3)=lambda(ix^D,2,3)
          nabla_v(3,1)=lambda(ix^D,3,1)
          nabla_v(3,2)=lambda(ix^D,3,2)
          nabla_v(3,3)=lambda(ix^D,3,3)
        end if
        ! dv_i/d_j + dv_j/d_i
        lambda(ix^D,1,2)=(lambda(ix^D,1,2)+lambda(ix^D,2,1))*vc_mus(ix^D)
        lambda(ix^D,1,3)=(lambda(ix^D,1,3)+lambda(ix^D,3,1))*vc_mus(ix^D)
        lambda(ix^D,2,3)=(lambda(ix^D,2,3)+lambda(ix^D,3,2))*vc_mus(ix^D)
        lambda(ix^D,2,1)=lambda(ix^D,1,2)
        lambda(ix^D,3,1)=lambda(ix^D,1,3)
        lambda(ix^D,3,2)=lambda(ix^D,2,3)
        divv23=two*third*(lambda(ix^D,1,1)+lambda(ix^D,2,2)+lambda(ix^D,3,3))
        lambda(ix^D,1,1)=(two*lambda(ix^D,1,1)-divv23)*vc_mus(ix^D)
        lambda(ix^D,2,2)=(two*lambda(ix^D,2,2)-divv23)*vc_mus(ix^D)
        lambda(ix^D,3,3)=(two*lambda(ix^D,3,3)-divv23)*vc_mus(ix^D)
        if(total_energy) then
          vlambda(ix^D,1)=wp(ix^D,v1_)*lambda(ix^D,1,1)+wp(ix^D,v2_)*lambda(ix^D,2,1)+wp(ix^D,v3_)*lambda(ix^D,3,1)
          vlambda(ix^D,2)=wp(ix^D,v1_)*lambda(ix^D,1,2)+wp(ix^D,v2_)*lambda(ix^D,2,2)+wp(ix^D,v3_)*lambda(ix^D,3,2)
          vlambda(ix^D,3)=wp(ix^D,v1_)*lambda(ix^D,1,3)+wp(ix^D,v2_)*lambda(ix^D,2,3)+wp(ix^D,v3_)*lambda(ix^D,3,3)
        end if
        if(phys_internal_e) then
          w(ix^D,e_)=w(ix^D,e_)+&
             qdt*(lambda(ix^D,1,1)*nabla_v(1,1)+lambda(ix^D,2,1)*nabla_v(2,1)+lambda(ix^D,3,1)*nabla_v(3,1)&
                 +lambda(ix^D,1,2)*nabla_v(1,2)+lambda(ix^D,2,2)*nabla_v(2,2)+lambda(ix^D,3,2)*nabla_v(3,2)&
                 +lambda(ix^D,1,3)*nabla_v(1,3)+lambda(ix^D,2,3)*nabla_v(2,3)+lambda(ix^D,3,3)*nabla_v(3,3))
        end if
      {end do\}
      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        w(ix^D,mom(1))=((lambda(ix1+1,ix2,ix3,1,1)-lambda(ix1-1,ix2,ix3,1,1))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))+&
                        (lambda(ix1,ix2+1,ix3,2,1)-lambda(ix1,ix2-1,ix3,2,1))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))*invr(ix^D)+&
                        (lambda(ix1,ix2,ix3+1,3,1)-lambda(ix1,ix2,ix3-1,3,1))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))*invrsin(ix^D)+&
                     two*lambda(ix^D,1,1)*invr(ix^D)+lambda(ix^D,2,1)*invr(ix^D)*tctan(ix^D))*qdt+w(ix^D,mom(1))
        w(ix^D,mom(2))=((lambda(ix1+1,ix2,ix3,1,2)-lambda(ix1-1,ix2,ix3,1,2))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))+&
                        (lambda(ix1,ix2+1,ix3,2,2)-lambda(ix1,ix2-1,ix3,2,2))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))*invr(ix^D)+&
                        (lambda(ix1,ix2,ix3+1,3,2)-lambda(ix1,ix2,ix3-1,3,2))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))*invrsin(ix^D)+&
                     two*lambda(ix^D,1,2)*invr(ix^D)+lambda(ix^D,2,2)*invr(ix^D)*tctan(ix^D)-lambda(ix^D,3,3)*tctan(ix^D)*invr(ix^D)**2)*qdt+&
                     w(ix^D,mom(2))
        w(ix^D,mom(3))=((lambda(ix1+1,ix2,ix3,1,3)-lambda(ix1-1,ix2,ix3,1,3))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))+&
                        (lambda(ix1,ix2+1,ix3,2,3)-lambda(ix1,ix2-1,ix3,2,3))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))*invr(ix^D)+&
                        (lambda(ix1,ix2,ix3+1,3,3)-lambda(ix1,ix2,ix3-1,3,3))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))*invrsin(ix^D)+&
                     two*lambda(ix^D,1,3)*invr(ix^D)+lambda(ix^D,2,3)*(invr(ix^D)+invr(ix^D)**2)*tctan(ix^D))*qdt+w(ix^D,mom(3))
      {end do\}
      }
      {^IFTWOD
      {do ix^DB=ixmin^DB,ixmax^DB\}
        invr(ix1,ix2)=1.d0/x(ix1,ix2,1)
        ! idim=1, idir=1
        lambda(ix1,ix2,1,1)=(wp(ix1+1,ix2,v1_)-wp(ix1-1,ix2,v1_))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))
        ! idim=1, idir=2
        lambda(ix1,ix2,1,2)=(wp(ix1+1,ix2,v2_)-wp(ix1-1,ix2,v2_))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))
        ! idim=2, idir=1
        lambda(ix1,ix2,2,1)=((wp(ix1,ix2+1,v1_)-wp(ix1,ix2-1,v1_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))-wp(ix1,ix2,v2_))*invr(ix1,ix2)
        ! idim=2, idir=2
        lambda(ix1,ix2,2,2)=((wp(ix1,ix2+1,v2_)-wp(ix1,ix2-1,v2_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))+wp(ix1,ix2,v1_))*invr(ix1,ix2)
        tcos=dcos(x(ix1,ix2,2))
        tsin=dsin(x(ix1,ix2,2))
        tctan(ix1,ix2)=tcos/tsin
        if(ndir==2) then
          if(phys_internal_e) then
            nabla_v(1,1)=lambda(ix^D,1,1)
            nabla_v(1,2)=lambda(ix^D,1,2)
            nabla_v(2,1)=lambda(ix^D,2,1)
            nabla_v(2,2)=lambda(ix^D,2,2)
          end if
          ! dv_i/d_j + dv_j/d_i
          lambda(ix1,ix2,1,2)=(lambda(ix1,ix2,1,2)+lambda(ix1,ix2,2,1))*vc_mus(ix^D)
          lambda(ix1,ix2,2,1)=lambda(ix1,ix2,1,2)
          divv23=two*third*(lambda(ix1,ix2,1,1)+lambda(ix1,ix2,2,2))
          lambda(ix1,ix2,1,1)=(two*lambda(ix1,ix2,1,1)-divv23)*vc_mus(ix^D)
          lambda(ix1,ix2,2,2)=(two*lambda(ix1,ix2,2,2)-divv23)*vc_mus(ix^D)
          if(total_energy) then
            vlambda(ix1,ix2,1)=wp(ix1,ix2,v1_)*lambda(ix1,ix2,1,1)+wp(ix1,ix2,v2_)*lambda(ix1,ix2,2,1)
            vlambda(ix1,ix2,2)=wp(ix1,ix2,v1_)*lambda(ix1,ix2,1,2)+wp(ix1,ix2,v2_)*lambda(ix1,ix2,2,2)
          end if
          if(phys_internal_e) then
            w(ix^D,e_)=w(ix^D,e_)+&
               qdt*(lambda(ix^D,1,1)*nabla_v(1,1)+lambda(ix^D,2,1)*nabla_v(2,1)&
                   +lambda(ix^D,1,2)*nabla_v(1,2)+lambda(ix^D,2,2)*nabla_v(2,2))
          end if
        else
          invrsin(ix1,ix2)=invr(ix1,ix2)/tsin
          ! idim=1, idir=3
          lambda(ix1,ix2,1,3)=(wp(ix1+1,ix2,v3_)-wp(ix1-1,ix2,v3_))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))
          ! idim=2, idir=3
          lambda(ix1,ix2,2,3)=((wp(ix1,ix2+1,v3_)-wp(ix1,ix2-1,v3_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))+wp(ix1,ix2,v2_)*tcos)*invr(ix1,ix2)
          ! idim=3, idir=1
          lambda(ix1,ix2,3,1)=-wp(ix1,ix2,v3_)*tsin*invrsin(ix1,ix2)
          ! idim=3, idir=2
          lambda(ix1,ix2,3,2)=-wp(ix1,ix2,v3_)*tcos*invrsin(ix1,ix2)
          ! idim=3, idir=3
          lambda(ix1,ix2,3,3)=(wp(ix1,ix2,v1_)*tsin+wp(ix1,ix2,v2_)*tcos)*invrsin(ix1,ix2)
          if(phys_internal_e) then
            nabla_v(1,1)=lambda(ix^D,1,1)
            nabla_v(1,2)=lambda(ix^D,1,2)
            nabla_v(1,3)=lambda(ix^D,1,3)
            nabla_v(2,1)=lambda(ix^D,2,1)
            nabla_v(2,2)=lambda(ix^D,2,2)
            nabla_v(2,3)=lambda(ix^D,2,3)
            nabla_v(3,1)=lambda(ix^D,3,1)
            nabla_v(3,2)=lambda(ix^D,3,2)
            nabla_v(3,3)=lambda(ix^D,3,3)
          end if
          ! dv_i/d_j + dv_j/d_i
          lambda(ix1,ix2,1,2)=(lambda(ix1,ix2,1,2)+lambda(ix1,ix2,2,1))*vc_mus(ix^D)
          lambda(ix1,ix2,1,3)=(lambda(ix1,ix2,1,3)+lambda(ix1,ix2,3,1))*vc_mus(ix^D)
          lambda(ix1,ix2,2,3)=(lambda(ix1,ix2,2,3)+lambda(ix1,ix2,3,2))*vc_mus(ix^D)
          lambda(ix1,ix2,2,1)=lambda(ix1,ix2,1,2)
          lambda(ix1,ix2,3,1)=lambda(ix1,ix2,1,3)
          lambda(ix1,ix2,3,2)=lambda(ix1,ix2,2,3)
          divv23=two*third*(lambda(ix^D,1,1)+lambda(ix^D,2,2)+lambda(ix^D,3,3))
          lambda(ix1,ix2,1,1)=(two*lambda(ix1,ix2,1,1)-divv23)*vc_mus(ix^D)
          lambda(ix1,ix2,2,2)=(two*lambda(ix1,ix2,2,2)-divv23)*vc_mus(ix^D)
          lambda(ix1,ix2,3,3)=(two*lambda(ix1,ix2,3,3)-divv23)*vc_mus(ix^D)
          if(total_energy) then
            vlambda(ix1,ix2,1)=wp(ix1,ix2,v1_)*lambda(ix1,ix2,1,1)+wp(ix1,ix2,v2_)*lambda(ix1,ix2,2,1)+wp(ix1,ix2,v3_)*lambda(ix1,ix2,3,1)
            vlambda(ix1,ix2,2)=wp(ix1,ix2,v1_)*lambda(ix1,ix2,1,2)+wp(ix1,ix2,v2_)*lambda(ix1,ix2,2,2)+wp(ix1,ix2,v3_)*lambda(ix1,ix2,3,2)
            vlambda(ix1,ix2,3)=wp(ix1,ix2,v1_)*lambda(ix1,ix2,1,3)+wp(ix1,ix2,v2_)*lambda(ix1,ix2,2,3)+wp(ix1,ix2,v3_)*lambda(ix1,ix2,3,3)
          end if
          if(phys_internal_e) then
            w(ix^D,e_)=w(ix^D,e_)+&
               qdt*(lambda(ix^D,1,1)*nabla_v(1,1)+lambda(ix^D,2,1)*nabla_v(2,1)+lambda(ix^D,3,1)*nabla_v(3,1)&
                   +lambda(ix^D,1,2)*nabla_v(1,2)+lambda(ix^D,2,2)*nabla_v(2,2)+lambda(ix^D,3,2)*nabla_v(3,2)&
                   +lambda(ix^D,1,3)*nabla_v(1,3)+lambda(ix^D,2,3)*nabla_v(2,3)+lambda(ix^D,3,3)*nabla_v(3,3))
          end if
        end if
      {end do\}
      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if(ndir==2) then
          w(ix1,ix2,mom(1))=((lambda(ix1+1,ix2,1,1)-lambda(ix1-1,ix2,1,1))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                             (lambda(ix1,ix2+1,2,1)-lambda(ix1,ix2-1,2,1))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix1,ix2)+&
                          two*lambda(ix1,ix2,1,1)*invr(ix1,ix2)+lambda(ix1,ix2,2,1)*invr(ix1,ix2)*tctan(ix1,ix2))*qdt+w(ix1,ix2,mom(1))
          w(ix1,ix2,mom(2))=((lambda(ix1+1,ix2,1,2)-lambda(ix1-1,ix2,1,2))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                             (lambda(ix1,ix2+1,2,2)-lambda(ix1,ix2-1,2,2))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix1,ix2)+&
                          two*lambda(ix1,ix2,1,2)*invr(ix1,ix2)+lambda(ix1,ix2,2,2)*invr(ix1,ix2)*tctan(ix1,ix2))*qdt+w(ix1,ix2,mom(2))
        else
          w(ix1,ix2,mom(1))=((lambda(ix1+1,ix2,1,1)-lambda(ix1-1,ix2,1,1))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                             (lambda(ix1,ix2+1,2,1)-lambda(ix1,ix2-1,2,1))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix1,ix2)+&
                          two*lambda(ix1,ix2,1,1)*invr(ix1,ix2)+lambda(ix1,ix2,2,1)*invr(ix1,ix2)*tctan(ix1,ix2))*qdt+w(ix1,ix2,mom(1))
          w(ix1,ix2,mom(2))=((lambda(ix1+1,ix2,1,2)-lambda(ix1-1,ix2,1,2))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                             (lambda(ix1,ix2+1,2,2)-lambda(ix1,ix2-1,2,2))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix1,ix2)+&
                          two*lambda(ix1,ix2,1,2)*invr(ix1,ix2)+lambda(ix1,ix2,2,2)*invr(ix1,ix2)*tctan(ix1,ix2)-&
                              lambda(ix1,ix2,3,3)*tctan(ix1,ix2)*invr(ix1,ix2)**2)*qdt+w(ix1,ix2,mom(2))
          w(ix1,ix2,mom(3))=((lambda(ix1+1,ix2,1,3)-lambda(ix1-1,ix2,1,3))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                             (lambda(ix1,ix2+1,2,3)-lambda(ix1,ix2-1,2,3))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix1,ix2)+&
                          two*lambda(ix1,ix2,1,3)*invr(ix1,ix2)+lambda(ix1,ix2,2,3)*(invr(ix1,ix2)+invr(ix1,ix2)**2)*tctan(ix1,ix2))*qdt+&
                          w(ix1,ix2,mom(3))
        end if
      {end do\}
      }
      {^IFONED
      do ix1=ixmin1,ixmax1
        invr(ix1)=1.d0/x(ix1,1)
        ! idim=1, idir=1
        lambda(ix1,1,1)=(wp(ix1+1,v1_)-wp(ix1-1,v1_))/(x(ix1+1,1)-x(ix1-1,1))
        if(ndir==1) then
          if(phys_internal_e) then
            nabla_v(1,1)=lambda(ix^D,1,1)
          end if
          divv23=two*third*lambda(ix1,1,1)
          lambda(ix1,1,1)=(two*lambda(ix1,1,1)-divv23)*vc_mus(ix^D)
          if(total_energy) then
            vlambda(ix1,1)=wp(ix1,v1_)*lambda(ix1,1,1)
          end if
          if(phys_internal_e) then
            w(ix^D,e_)=w(ix^D,e_)+&
               qdt*lambda(ix^D,1,1)*nabla_v(1,1)
          end if
        else if(ndir==2) then
          ! idim=1, idir=2
          lambda(ix1,1,2)=(wp(ix1+1,v2_)-wp(ix1-1,v2_))/(x(ix1+1,1)-x(ix1-1,1))
          lambda(ix1,2,1)=0.d0
          lambda(ix1,2,2)=-divv23
          if(phys_internal_e) then
            nabla_v(1,1)=lambda(ix^D,1,1)
            nabla_v(1,2)=lambda(ix^D,1,2)
            nabla_v(2,1)=lambda(ix^D,2,1)
            nabla_v(2,2)=lambda(ix^D,2,2)
          end if
          ! dv_i/d_j + dv_j/d_i
          lambda(ix1,1,2)=lambda(ix1,1,2)*vc_mus(ix^D)
          lambda(ix1,2,1)=lambda(ix1,1,2)
          divv23=two*third*(lambda(ix1,1,1)+lambda(ix1,2,2))
          lambda(ix1,1,1)=(two*lambda(ix1,1,1)-divv23)*vc_mus(ix^D)
          lambda(ix1,2,2)=(two*lambda(ix1,2,2)-divv23)*vc_mus(ix^D)
          if(total_energy) then
            vlambda(ix1,1)=wp(ix1,v1_)*lambda(ix1,1,1)+wp(ix1,v2_)*lambda(ix1,2,1)
            vlambda(ix1,2)=wp(ix1,v1_)*lambda(ix1,1,2)+wp(ix1,v2_)*lambda(ix1,2,2)
          end if
          if(phys_internal_e) then
            w(ix^D,e_)=w(ix^D,e_)+&
               qdt*(lambda(ix^D,1,1)*nabla_v(1,1)+lambda(ix^D,2,1)*nabla_v(2,1)&
                   +lambda(ix^D,1,2)*nabla_v(1,2)+lambda(ix^D,2,2)*nabla_v(2,2))
          end if
        else
          tcos=0.d0
          tsin=1.d0
          tctan(ix1)=tcos/tsin
          invrsin(ix1)=invr(ix1)/tsin
          lambda(ix1,1,2)=(wp(ix1+1,v2_)-wp(ix1-1,v2_))/(x(ix1+1,1)-x(ix1-1,1))
          lambda(ix1,1,3)=(wp(ix1+1,v3_)-wp(ix1-1,v3_))/(x(ix1+1,1)-x(ix1-1,1))
          lambda(ix1,2,1)=0.d0
          lambda(ix1,2,2)=-divv23
          lambda(ix1,2,3)=wp(ix1,v2_)*tcos*invr(ix1)
          lambda(ix1,3,1)=-wp(ix1,v3_)*tsin*invrsin(ix1)
          lambda(ix1,3,2)=-wp(ix1,v3_)*tcos*invrsin(ix1)
          lambda(ix1,3,3)=(wp(ix1,v1_)*tsin+wp(ix1,v2_)*tcos)*invrsin(ix1)
          if(phys_internal_e) then
            nabla_v(1,1)=lambda(ix^D,1,1)
            nabla_v(1,2)=lambda(ix^D,1,2)
            nabla_v(1,3)=lambda(ix^D,1,3)
            nabla_v(2,1)=lambda(ix^D,2,1)
            nabla_v(2,2)=lambda(ix^D,2,2)
            nabla_v(2,3)=lambda(ix^D,2,3)
            nabla_v(3,1)=lambda(ix^D,3,1)
            nabla_v(3,2)=lambda(ix^D,3,2)
            nabla_v(3,3)=lambda(ix^D,3,3)
          end if
          ! dv_i/d_j + dv_j/d_i
          lambda(ix1,1,2)=(lambda(ix1,1,2)+lambda(ix1,2,1))*vc_mus(ix^D)
          lambda(ix1,1,3)=(lambda(ix1,1,3)+lambda(ix1,3,1))*vc_mus(ix^D)
          lambda(ix1,2,3)=(lambda(ix1,2,3)+lambda(ix1,3,2))*vc_mus(ix^D)
          lambda(ix1,3,1)=lambda(ix1,1,3)
          lambda(ix1,2,3)=lambda(ix1,3,2)
          divv23=two*third*(lambda(ix1,1,1)+lambda(ix1,2,2)+lambda(ix1,3,3))
          lambda(ix1,1,1)=(two*lambda(ix1,1,1)-divv23)*vc_mus(ix^D)
          lambda(ix1,2,2)=(two*lambda(ix1,2,2)-divv23)*vc_mus(ix^D)
          lambda(ix1,3,3)=(two*lambda(ix1,3,3)-divv23)*vc_mus(ix^D)
          if(total_energy) then
            vlambda(ix1,1)=wp(ix1,v1_)*lambda(ix1,1,1)+wp(ix1,v2_)*lambda(ix1,2,1)+wp(ix1,v3_)*lambda(ix1,3,1)
            vlambda(ix1,2)=wp(ix1,v1_)*lambda(ix1,1,2)+wp(ix1,v2_)*lambda(ix1,2,2)+wp(ix1,v3_)*lambda(ix1,3,2)
            vlambda(ix1,3)=wp(ix1,v1_)*lambda(ix1,1,3)+wp(ix1,v2_)*lambda(ix1,2,3)+wp(ix1,v3_)*lambda(ix1,3,3)
          end if
          if(phys_internal_e) then
            w(ix^D,e_)=w(ix^D,e_)+&
               qdt*(lambda(ix^D,1,1)*nabla_v(1,1)+lambda(ix^D,2,1)*nabla_v(2,1)+lambda(ix^D,3,1)*nabla_v(3,1)&
                   +lambda(ix^D,1,2)*nabla_v(1,2)+lambda(ix^D,2,2)*nabla_v(2,2)+lambda(ix^D,3,2)*nabla_v(3,2)&
                   +lambda(ix^D,1,3)*nabla_v(1,3)+lambda(ix^D,2,3)*nabla_v(2,3)+lambda(ix^D,3,3)*nabla_v(3,3))
          end if
        end if
      end do
      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
      do ix1=ixOmin1,ixOmax1
        if(ndir==1) then
          w(ix1,mom(1))=((lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))+&
                      two*lambda(ix1,1,1)*invr(ix1))*qdt+w(ix1,mom(1))
        else if(ndir==2) then
          w(ix1,mom(1))=((lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))+&
                      two*lambda(ix1,1,1)*invr(ix1)+lambda(ix1,2,1)*invr(ix1)*tctan(ix1))*qdt+w(ix1,mom(1))
          w(ix1,mom(2))=((lambda(ix1+1,1,2)-lambda(ix1-1,1,2))/(x(ix1+1,1)-x(ix1-1,1))+&
                      two*lambda(ix1,1,2)*invr(ix1)+lambda(ix1,2,2)*invr(ix1)*tctan(ix1))*qdt+w(ix1,mom(2))
        else
          w(ix1,mom(1))=((lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))+&
                      two*lambda(ix1,1,1)*invr(ix1)+lambda(ix1,2,1)*invr(ix1)*tctan(ix1))*qdt+w(ix1,mom(1))
          w(ix1,mom(2))=((lambda(ix1+1,1,2)-lambda(ix1-1,1,2))/(x(ix1+1,1)-x(ix1-1,1))+&
                      two*lambda(ix1,1,2)*invr(ix1)+lambda(ix1,2,2)*invr(ix1)*tctan(ix1)-lambda(ix1,3,3)*tctan(ix1)*invr(ix1)**2)*qdt+&
                      w(ix1,mom(2))
          w(ix1,mom(3))=((lambda(ix1+1,1,3)-lambda(ix1-1,1,3))/(x(ix1+1,1)-x(ix1-1,1))+&
                      two*lambda(ix1,1,3)*invr(ix1)+lambda(ix1,2,3)*(invr(ix1)+invr(ix1)**2)*tctan(ix1))*qdt+w(ix1,mom(3))
        end if
      end do
      }
      if(total_energy) then
        ! de/dt= +div(v.dot.[mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v *kr])
        ! thus e=e+d_i v_j tensor_ji
        call divvector(vlambda,ixI^L,ixO^L,invr)
        w(ixO^S,e_)=w(ixO^S,e_)+invr(ixO^S)*qdt
      end if

    end if

  end subroutine viscosity_add_source_sphere

  subroutine viscosity_add_source_cylinder(qdt,ixI^L,ixO^L,wCT,wp,w,x,&
       energy,qsourcesplit,active)
  ! Add viscosity source in isotropic Newtonian fluids to w within ixO
  ! neglecting bulk viscosity
  ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods, only: usr_set_viscosity

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw), wp(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active

    double precision:: lambda(ixI^S,ndir,ndir),vlambda(ixI^S,ndir),invr(ixI^S),nabla_v(ndir,ndir),vc_mus(ixI^S)
    double precision :: divv23
    integer:: ix^L,ix^D
    logical :: total_energy=.false.

    if(qsourcesplit .eqv. vc_split) then
      active = .true.
      if(energy.and..not.phys_internal_e) total_energy=.true.
      ! standard case, textbook viscosity
      ! Calculating viscosity sources
      ! involves second derivatives, two extra layers
      ix^L=ixO^L^LADD2;
      if({ ixImin^D>ixmin^D .or. ixImax^D<ixmax^D|.or.})&
        call mpistop("error for viscous source addition, 2 layers needed")
      ix^L=ixO^L^LADD1;

      ! construct lambda tensor: lambda_ij = gradv_ij + gradv_ji
      ! initialize
      if(associated(usr_set_viscosity)) then
        call usr_set_viscosity(ixI^L,ix^L,wp,x,vc_mus)
      else
        vc_mus=vc_mu
      end if
      {^IFTHREED
      {do ix^DB=ixmin^DB,ixmax^DB\}
        invr(ix^D)=1.d0/x(ix^D,1)
        ! idim=1, idir=1
        lambda(ix^D,1,1)=(wp(ix1+1,ix2,ix3,v1_)-wp(ix1-1,ix2,ix3,v1_))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))
        ! idim=1, idir=2
        lambda(ix^D,1,2)=(wp(ix1+1,ix2,ix3,v2_)-wp(ix1-1,ix2,ix3,v2_))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))
        ! idim=1, idir=3
        lambda(ix^D,1,3)=(wp(ix1+1,ix2,ix3,v3_)-wp(ix1-1,ix2,ix3,v3_))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))
        ! idim=2, idir=1
        lambda(ix^D,2,1)=((wp(ix1,ix2+1,ix3,v1_)-wp(ix1,ix2-1,ix3,v1_))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))-wp(ix^D,v2_))*invr(ix^D)
        ! idim=2, idir=2
        lambda(ix^D,2,2)=((wp(ix1,ix2+1,ix3,v2_)-wp(ix1,ix2-1,ix3,v2_))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))+wp(ix^D,v1_))*invr(ix^D)
        ! idim=2, idir=3
        lambda(ix^D,2,3)=(wp(ix1,ix2+1,ix3,v3_)-wp(ix1,ix2-1,ix3,v3_))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))*invr(ix^D)
        ! idim=3, idir=1
        lambda(ix^D,3,1)=(wp(ix1,ix2,ix3+1,v1_)-wp(ix1,ix2,ix3-1,v1_))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))*invr(ix^D)
        ! idim=3, idir=2
        lambda(ix^D,3,2)=(wp(ix1,ix2,ix3+1,v2_)-wp(ix1,ix2,ix3-1,v2_))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))*invr(ix^D)
        ! idim=3, idir=3
        lambda(ix^D,3,3)=(wp(ix1,ix2,ix3+1,v3_)-wp(ix1,ix2,ix3-1,v3_))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))*invr(ix^D)
        if(phys_internal_e) then
          nabla_v(1,1)=lambda(ix^D,1,1)
          nabla_v(1,2)=lambda(ix^D,1,2)
          nabla_v(1,3)=lambda(ix^D,1,3)
          nabla_v(2,1)=lambda(ix^D,2,1)
          nabla_v(2,2)=lambda(ix^D,2,2)
          nabla_v(2,3)=lambda(ix^D,2,3)
          nabla_v(3,1)=lambda(ix^D,3,1)
          nabla_v(3,2)=lambda(ix^D,3,2)
          nabla_v(3,3)=lambda(ix^D,3,3)
        end if
        ! dv_i/d_j + dv_j/d_i
        lambda(ix^D,1,2)=(lambda(ix^D,1,2)+lambda(ix^D,2,1))*vc_mus(ix^D)
        lambda(ix^D,1,3)=(lambda(ix^D,1,3)+lambda(ix^D,3,1))*vc_mus(ix^D)
        lambda(ix^D,2,3)=(lambda(ix^D,2,3)+lambda(ix^D,3,2))*vc_mus(ix^D)
        lambda(ix^D,2,1)=lambda(ix^D,1,2)
        lambda(ix^D,3,1)=lambda(ix^D,1,3)
        lambda(ix^D,3,2)=lambda(ix^D,2,3)
        divv23=two*third*(lambda(ix^D,1,1)+lambda(ix^D,2,2)+lambda(ix^D,3,3))
        lambda(ix^D,1,1)=(two*lambda(ix^D,1,1)-divv23)*vc_mus(ix^D)
        lambda(ix^D,2,2)=(two*lambda(ix^D,2,2)-divv23)*vc_mus(ix^D)
        lambda(ix^D,3,3)=(two*lambda(ix^D,3,3)-divv23)*vc_mus(ix^D)
        if(total_energy) then
          vlambda(ix^D,1)=wp(ix^D,v1_)*lambda(ix^D,1,1)+wp(ix^D,v2_)*lambda(ix^D,2,1)+wp(ix^D,v3_)*lambda(ix^D,3,1)
          vlambda(ix^D,2)=wp(ix^D,v1_)*lambda(ix^D,1,2)+wp(ix^D,v2_)*lambda(ix^D,2,2)+wp(ix^D,v3_)*lambda(ix^D,3,2)
          vlambda(ix^D,3)=wp(ix^D,v1_)*lambda(ix^D,1,3)+wp(ix^D,v2_)*lambda(ix^D,2,3)+wp(ix^D,v3_)*lambda(ix^D,3,3)
        end if
        if(phys_internal_e) then
          w(ix^D,e_)=w(ix^D,e_)+&
             qdt*(lambda(ix^D,1,1)*nabla_v(1,1)+lambda(ix^D,2,1)*nabla_v(2,1)+lambda(ix^D,3,1)*nabla_v(3,1)&
                 +lambda(ix^D,1,2)*nabla_v(1,2)+lambda(ix^D,2,2)*nabla_v(2,2)+lambda(ix^D,3,2)*nabla_v(3,2)&
                 +lambda(ix^D,1,3)*nabla_v(1,3)+lambda(ix^D,2,3)*nabla_v(2,3)+lambda(ix^D,3,3)*nabla_v(3,3))
        end if
      {end do\}
      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        w(ix^D,mom(1))=((lambda(ix1+1,ix2,ix3,1,1)-lambda(ix1-1,ix2,ix3,1,1))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))+&
                        (lambda(ix1,ix2+1,ix3,2,1)-lambda(ix1,ix2-1,ix3,2,1))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))*invr(ix^D)+&
                        (lambda(ix1,ix2,ix3+1,3,1)-lambda(ix1,ix2,ix3-1,3,1))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))+&
                        (lambda(ix^D,1,1)-lambda(ix^D,2,2))*invr(ix^D))*qdt+w(ix^D,mom(1))
        w(ix^D,mom(2))=((lambda(ix1+1,ix2,ix3,1,2)-lambda(ix1-1,ix2,ix3,1,2))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))+&
                        (lambda(ix1,ix2+1,ix3,2,2)-lambda(ix1,ix2-1,ix3,2,2))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))*invr(ix^D)+&
                        (lambda(ix1,ix2,ix3+1,3,2)-lambda(ix1,ix2,ix3-1,3,2))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))+&
                        two*lambda(ix^D,1,2)*invr(ix^D))*qdt+w(ix^D,mom(2))
        w(ix^D,mom(3))=((lambda(ix1+1,ix2,ix3,1,3)-lambda(ix1-1,ix2,ix3,1,3))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))+&
                        (lambda(ix1,ix2+1,ix3,2,3)-lambda(ix1,ix2-1,ix3,2,3))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))*invr(ix^D)+&
                        (lambda(ix1,ix2,ix3+1,3,3)-lambda(ix1,ix2,ix3-1,3,3))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))+&
                         lambda(ix^D,1,3)*invr(ix^D))*qdt+w(ix^D,mom(3))
      {end do\}
      }
      {^IFTWOD
      {do ix^DB=ixmin^DB,ixmax^DB\}
        invr(ix^D)=1.d0/x(ix^D,1)
        ! idim=1, idir=1
        lambda(ix^D,1,1)=(wp(ix1+1,ix2,v1_)-wp(ix1-1,ix2,v1_))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))
        ! idim=1, idir=2
        lambda(ix^D,1,2)=(wp(ix1+1,ix2,v2_)-wp(ix1-1,ix2,v2_))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))
        if(phi_==2) then
          ! idim=2, idir=1
          lambda(ix^D,2,1)=((wp(ix1,ix2+1,v1_)-wp(ix1,ix2-1,v1_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))-wp(ix^D,v2_))*invr(ix^D)
          ! idim=2, idir=2
          lambda(ix^D,2,2)=((wp(ix1,ix2+1,v2_)-wp(ix1,ix2-1,v2_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))+wp(ix^D,v1_))*invr(ix^D)
          if(ndir==3) then
            lambda(ix^D,1,3)=(wp(ix1+1,ix2,v3_)-wp(ix1-1,ix2,v3_))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))
            ! idim=2, idir=3
            lambda(ix^D,2,3)=(wp(ix1,ix2+1,v3_)-wp(ix1,ix2-1,v3_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix^D)
            ! idim=3, idir=1
            lambda(ix^D,3,1)=0.d0
            ! idim=3, idir=2
            lambda(ix^D,3,2)=0.d0
            ! idim=3, idir=3
            lambda(ix^D,3,3)=0.d0
          end if
        else
          ! idim=2, idir=1
          lambda(ix^D,2,1)=(wp(ix1,ix2+1,v1_)-wp(ix1,ix2-1,v1_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix^D)
          ! idim=2, idir=2
          lambda(ix^D,2,2)=(wp(ix1,ix2+1,v2_)-wp(ix1,ix2-1,v2_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix^D)
          if(ndir==3) then
            lambda(ix^D,1,3)=(wp(ix1+1,ix2,v3_)-wp(ix1-1,ix2,v3_))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))
            ! idim=2, idir=3
            lambda(ix^D,2,3)=(wp(ix1,ix2+1,v3_)-wp(ix1,ix2-1,v3_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix^D)
            ! idim=3, idir=1
            lambda(ix^D,3,1)=-wp(ix^D,v3_)*invr(ix^D)
            ! idim=3, idir=2
            lambda(ix^D,3,2)=0.d0
            ! idim=3, idir=3
            lambda(ix^D,3,3)=wp(ix^D,v1_)*invr(ix^D)
          end if
        end if
        if(ndir==3) then
          if(phys_internal_e) then
            nabla_v(1,1)=lambda(ix^D,1,1)
            nabla_v(1,2)=lambda(ix^D,1,2)
            nabla_v(1,3)=lambda(ix^D,1,3)
            nabla_v(2,1)=lambda(ix^D,2,1)
            nabla_v(2,2)=lambda(ix^D,2,2)
            nabla_v(2,3)=lambda(ix^D,2,3)
            nabla_v(3,1)=lambda(ix^D,3,1)
            nabla_v(3,2)=lambda(ix^D,3,2)
            nabla_v(3,3)=lambda(ix^D,3,3)
          end if
          ! dv_i/d_j + dv_j/d_i
          lambda(ix^D,1,2)=(lambda(ix^D,1,2)+lambda(ix^D,2,1))*vc_mus(ix^D)
          lambda(ix^D,1,3)=(lambda(ix^D,1,3)+lambda(ix^D,3,1))*vc_mus(ix^D)
          lambda(ix^D,2,3)=(lambda(ix^D,2,3)+lambda(ix^D,3,2))*vc_mus(ix^D)
          lambda(ix^D,2,1)=lambda(ix^D,1,2)
          lambda(ix^D,3,1)=lambda(ix^D,1,3)
          lambda(ix^D,3,2)=lambda(ix^D,2,3)
          divv23=two*third*(lambda(ix^D,1,1)+lambda(ix^D,2,2)+lambda(ix^D,3,3))
          lambda(ix^D,1,1)=(two*lambda(ix^D,1,1)-divv23)*vc_mus(ix^D)
          lambda(ix^D,2,2)=(two*lambda(ix^D,2,2)-divv23)*vc_mus(ix^D)
          lambda(ix^D,3,3)=(two*lambda(ix^D,3,3)-divv23)*vc_mus(ix^D)
          if(total_energy) then
            vlambda(ix^D,1)=wp(ix^D,v1_)*lambda(ix^D,1,1)+wp(ix^D,v2_)*lambda(ix^D,2,1)+wp(ix^D,v3_)*lambda(ix^D,3,1)
            vlambda(ix^D,2)=wp(ix^D,v1_)*lambda(ix^D,1,2)+wp(ix^D,v2_)*lambda(ix^D,2,2)+wp(ix^D,v3_)*lambda(ix^D,3,2)
            vlambda(ix^D,3)=wp(ix^D,v1_)*lambda(ix^D,1,3)+wp(ix^D,v2_)*lambda(ix^D,2,3)+wp(ix^D,v3_)*lambda(ix^D,3,3)
          end if
          if(phys_internal_e) then
            w(ix^D,e_)=w(ix^D,e_)+&
               qdt*(lambda(ix^D,1,1)*nabla_v(1,1)+lambda(ix^D,2,1)*nabla_v(2,1)+lambda(ix^D,3,1)*nabla_v(3,1)&
                   +lambda(ix^D,1,2)*nabla_v(1,2)+lambda(ix^D,2,2)*nabla_v(2,2)+lambda(ix^D,3,2)*nabla_v(3,2)&
                   +lambda(ix^D,1,3)*nabla_v(1,3)+lambda(ix^D,2,3)*nabla_v(2,3)+lambda(ix^D,3,3)*nabla_v(3,3))
          end if
        else 
          if(phys_internal_e) then
            nabla_v(1,1)=lambda(ix^D,1,1)
            nabla_v(1,2)=lambda(ix^D,1,2)
            nabla_v(2,1)=lambda(ix^D,2,1)
            nabla_v(2,2)=lambda(ix^D,2,2)
          end if
          ! dv_i/d_j + dv_j/d_i
          lambda(ix^D,1,2)=(lambda(ix^D,1,2)+lambda(ix^D,2,1))*vc_mus(ix^D)
          lambda(ix^D,2,1)=lambda(ix^D,1,2)
          divv23=two*third*(lambda(ix^D,1,1)+lambda(ix^D,2,2))
          lambda(ix^D,1,1)=(two*lambda(ix^D,1,1)-divv23)*vc_mus(ix^D)
          lambda(ix^D,2,2)=(two*lambda(ix^D,2,2)-divv23)*vc_mus(ix^D)
          if(total_energy) then
            vlambda(ix^D,1)=wp(ix^D,v1_)*lambda(ix^D,1,1)+wp(ix^D,v2_)*lambda(ix^D,2,1)
            vlambda(ix^D,2)=wp(ix^D,v1_)*lambda(ix^D,1,2)+wp(ix^D,v2_)*lambda(ix^D,2,2)
          end if
          if(phys_internal_e) then
            w(ix^D,e_)=w(ix^D,e_)+&
               qdt*(lambda(ix^D,1,1)*nabla_v(1,1)+lambda(ix^D,2,1)*nabla_v(2,1)&
                   +lambda(ix^D,1,2)*nabla_v(1,2)+lambda(ix^D,2,2)*nabla_v(2,2))
          end if
        end if
      {end do\}
      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if(ndir==2) then
          if(phi_==2) then
            w(ix^D,mom(1))=((lambda(ix1+1,ix2,1,1)-lambda(ix1-1,ix2,1,1))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                            (lambda(ix1,ix2+1,2,1)-lambda(ix1,ix2-1,2,1))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix^D)+&
                            (lambda(ix^D,1,1)-lambda(ix^D,2,2))*invr(ix^D))*qdt+w(ix^D,mom(1))
            w(ix^D,mom(2))=((lambda(ix1+1,ix2,1,2)-lambda(ix1-1,ix2,1,2))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                            (lambda(ix1,ix2+1,2,2)-lambda(ix1,ix2-1,2,2))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix^D)+&
                         two*lambda(ix^D,1,2)*invr(ix^D))*qdt+w(ix^D,mom(2))
          else
            w(ix^D,mom(1))=((lambda(ix1+1,ix2,1,1)-lambda(ix1-1,ix2,1,1))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                            (lambda(ix1,ix2+1,2,1)-lambda(ix1,ix2-1,2,1))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))+&
                             lambda(ix^D,1,1)*invr(ix^D))*qdt+w(ix^D,mom(1))
            w(ix^D,mom(2))=((lambda(ix1+1,ix2,1,2)-lambda(ix1-1,ix2,1,2))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                            (lambda(ix1,ix2+1,2,2)-lambda(ix1,ix2-1,2,2))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))+&
                             lambda(ix^D,1,2)*invr(ix^D))*qdt+w(ix^D,mom(2))
          end if
        else
          if(phi_==2) then
            w(ix^D,mom(1))=((lambda(ix1+1,ix2,1,1)-lambda(ix1-1,ix2,1,1))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                            (lambda(ix1,ix2+1,2,1)-lambda(ix1,ix2-1,2,1))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix^D)+&
                            (lambda(ix^D,1,1)-lambda(ix^D,2,2))*invr(ix^D))*qdt+w(ix^D,mom(1))
            w(ix^D,mom(2))=((lambda(ix1+1,ix2,1,2)-lambda(ix1-1,ix2,1,2))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                            (lambda(ix1,ix2+1,2,2)-lambda(ix1,ix2-1,2,2))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix^D)+&
                            two*lambda(ix^D,1,2)*invr(ix^D))*qdt+w(ix^D,mom(2))
            w(ix^D,mom(3))=((lambda(ix1+1,ix2,1,3)-lambda(ix1-1,ix2,1,3))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                            (lambda(ix1,ix2+1,2,3)-lambda(ix1,ix2-1,2,3))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix^D)+&
                             lambda(ix^D,1,3)*invr(ix^D))*qdt+w(ix^D,mom(3))
          else
            w(ix^D,mom(1))=((lambda(ix1+1,ix2,1,1)-lambda(ix1-1,ix2,1,1))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                            (lambda(ix1,ix2+1,2,1)-lambda(ix1,ix2-1,2,1))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))+&
                            (lambda(ix^D,1,1)-lambda(ix^D,2,2))*invr(ix^D))*qdt+w(ix^D,mom(1))
            w(ix^D,mom(2))=((lambda(ix1+1,ix2,1,3)-lambda(ix1-1,ix2,1,3))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                            (lambda(ix1,ix2+1,2,3)-lambda(ix1,ix2-1,2,3))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))+&
                            two*lambda(ix^D,1,3)*invr(ix^D))*qdt+w(ix^D,mom(2))
            w(ix^D,mom(3))=((lambda(ix1+1,ix2,1,2)-lambda(ix1-1,ix2,1,2))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                            (lambda(ix1,ix2+1,2,2)-lambda(ix1,ix2-1,2,2))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))+&
                             lambda(ix^D,1,2)*invr(ix^D))*qdt+w(ix^D,mom(3))
          end if
        end if
      {end do\}
      }
      {^IFONED
      {do ix^DB=ixmin^DB,ixmax^DB\}
        invr(ix^D)=1.d0/x(ix^D,1)
        if(ndir==1) then
          ! idim=1, idir=1
          lambda(ix^D,1,1)=(wp(ix1+1,v1_)-wp(ix1-1,v1_))/(x(ix1+1,1)-x(ix1-1,1))
          if(phys_internal_e) then
            nabla_v(1,1)=lambda(ix^D,1,1)
          end if
          divv23=two*third*lambda(ix^D,1,1)
          lambda(ix^D,1,1)=(two*lambda(ix^D,1,1)-divv23)*vc_mus(ix^D)
          if(total_energy) then
            vlambda(ix^D,1)=wp(ix^D,v1_)*lambda(ix^D,1,1)
          end if
          if(phys_internal_e) then
            w(ix^D,e_)=w(ix^D,e_)+&
               qdt*lambda(ix^D,1,1)*nabla_v(1,1)
          end if
        else if(ndir==2) then
          ! idim=1, idir=1
          lambda(ix^D,1,1)=(wp(ix1+1,v1_)-wp(ix1-1,v1_))/(x(ix1+1,1)-x(ix1-1,1))
          ! idim=1, idir=2
          lambda(ix^D,1,2)=(wp(ix1+1,v2_)-wp(ix1-1,v2_))/(x(ix1+1,1)-x(ix1-1,1))
          if(phi_==2) then
            ! idim=2, idir=1
            lambda(ix^D,2,1)=-wp(ix^D,v2_)*invr(ix^D)
            ! idim=2, idir=2
            lambda(ix^D,2,2)=+wp(ix^D,v1_)*invr(ix^D)
          else
            ! idim=2, idir=1
            lambda(ix^D,2,1)=0.d0
            ! idim=2, idir=2
            lambda(ix^D,2,2)=0.d0
          end if
          if(phys_internal_e) then
            nabla_v(1,1)=lambda(ix^D,1,1)
            nabla_v(1,2)=lambda(ix^D,1,2)
            nabla_v(2,1)=lambda(ix^D,2,1)
            nabla_v(2,2)=lambda(ix^D,2,2)
          end if
          ! dv_i/d_j + dv_j/d_i
          lambda(ix^D,1,2)=(lambda(ix^D,1,2)+lambda(ix^D,2,1))*vc_mus(ix^D)
          lambda(ix^D,2,1)=lambda(ix^D,1,2)
          divv23=two*third*(lambda(ix^D,1,1)+lambda(ix^D,2,2))
          lambda(ix^D,1,1)=(two*lambda(ix^D,1,1)-divv23)*vc_mus(ix^D)
          lambda(ix^D,2,2)=(two*lambda(ix^D,2,2)-divv23)*vc_mus(ix^D)
          if(total_energy) then
            vlambda(ix^D,1)=wp(ix^D,v1_)*lambda(ix^D,1,1)+wp(ix^D,v2_)*lambda(ix^D,2,1)
            vlambda(ix^D,2)=wp(ix^D,v1_)*lambda(ix^D,1,2)+wp(ix^D,v2_)*lambda(ix^D,2,2)
          end if
          if(phys_internal_e) then
            w(ix^D,e_)=w(ix^D,e_)+&
               qdt*(lambda(ix^D,1,1)*nabla_v(1,1)+lambda(ix^D,2,1)*nabla_v(2,1)&
                   +lambda(ix^D,1,2)*nabla_v(1,2)+lambda(ix^D,2,2)*nabla_v(2,2))
          end if
        else
          ! idim=1, idir=1
          lambda(ix^D,1,1)=(wp(ix1+1,v1_)-wp(ix1-1,v1_))/(x(ix1+1,1)-x(ix1-1,1))
          ! idim=1, idir=2
          lambda(ix^D,1,2)=(wp(ix1+1,v2_)-wp(ix1-1,v2_))/(x(ix1+1,1)-x(ix1-1,1))
          ! idim=1, idir=3
          lambda(ix^D,1,3)=(wp(ix1+1,v3_)-wp(ix1-1,v3_))/(x(ix1+1,1)-x(ix1-1,1))
          if(phi_==2) then
            ! idim=2, idir=1
            lambda(ix^D,2,1)=-wp(ix^D,v2_)*invr(ix^D)
            ! idim=2, idir=2
            lambda(ix^D,2,2)=+wp(ix^D,v1_)*invr(ix^D)
            ! idim=2, idir=3
            lambda(ix^D,2,3)=0.d0
            ! idim=3, idir=1
            lambda(ix^D,3,1)=0.d0
            ! idim=3, idir=2
            lambda(ix^D,3,2)=0.d0
            ! idim=3, idir=3
            lambda(ix^D,3,3)=0.d0
          else
            ! idim=2, idir=1
            lambda(ix^D,2,1)=0.d0
            ! idim=2, idir=2
            lambda(ix^D,2,2)=0.d0
            ! idim=2, idir=3
            lambda(ix^D,2,3)=0.d0
            ! idim=3, idir=1
            lambda(ix^D,3,1)=-wp(ix^D,v3_)*invr(ix^D)
            ! idim=3, idir=2
            lambda(ix^D,3,2)=0.d0
            ! idim=3, idir=3
            lambda(ix^D,3,3)=wp(ix^D,v1_)*invr(ix^D)
          end if
          if(phys_internal_e) then
            nabla_v(1,1)=lambda(ix^D,1,1)
            nabla_v(1,2)=lambda(ix^D,1,2)
            nabla_v(1,3)=lambda(ix^D,1,3)
            nabla_v(2,1)=lambda(ix^D,2,1)
            nabla_v(2,2)=lambda(ix^D,2,2)
            nabla_v(2,3)=lambda(ix^D,2,3)
            nabla_v(3,1)=lambda(ix^D,3,1)
            nabla_v(3,2)=lambda(ix^D,3,2)
            nabla_v(3,3)=lambda(ix^D,3,3)
          end if
          ! dv_i/d_j + dv_j/d_i
          lambda(ix^D,1,2)=(lambda(ix^D,1,2)+lambda(ix^D,2,1))*vc_mus(ix^D)
          lambda(ix^D,1,3)=(lambda(ix^D,1,3)+lambda(ix^D,3,1))*vc_mus(ix^D)
          lambda(ix^D,2,3)=(lambda(ix^D,2,3)+lambda(ix^D,3,2))*vc_mus(ix^D)
          lambda(ix^D,2,1)=lambda(ix^D,1,2)
          lambda(ix^D,3,1)=lambda(ix^D,1,3)
          lambda(ix^D,3,2)=lambda(ix^D,2,3)
          divv23=two*third*(lambda(ix^D,1,1)+lambda(ix^D,2,2)+lambda(ix^D,3,3))
          lambda(ix^D,1,1)=(two*lambda(ix^D,1,1)-divv23)*vc_mus(ix^D)
          lambda(ix^D,2,2)=(two*lambda(ix^D,2,2)-divv23)*vc_mus(ix^D)
          lambda(ix^D,3,3)=(two*lambda(ix^D,3,3)-divv23)*vc_mus(ix^D)
          if(total_energy) then
            vlambda(ix^D,1)=wp(ix^D,v1_)*lambda(ix^D,1,1)+wp(ix^D,v2_)*lambda(ix^D,2,1)+wp(ix^D,v3_)*lambda(ix^D,3,1)
            vlambda(ix^D,2)=wp(ix^D,v1_)*lambda(ix^D,1,2)+wp(ix^D,v2_)*lambda(ix^D,2,2)+wp(ix^D,v3_)*lambda(ix^D,3,2)
            vlambda(ix^D,3)=wp(ix^D,v1_)*lambda(ix^D,1,3)+wp(ix^D,v2_)*lambda(ix^D,2,3)+wp(ix^D,v3_)*lambda(ix^D,3,3)
          end if
          if(phys_internal_e) then
            w(ix^D,e_)=w(ix^D,e_)+&
               qdt*(lambda(ix^D,1,1)*nabla_v(1,1)+lambda(ix^D,2,1)*nabla_v(2,1)+lambda(ix^D,3,1)*nabla_v(3,1)&
                   +lambda(ix^D,1,2)*nabla_v(1,2)+lambda(ix^D,2,2)*nabla_v(2,2)+lambda(ix^D,3,2)*nabla_v(3,2)&
                   +lambda(ix^D,1,3)*nabla_v(1,3)+lambda(ix^D,2,3)*nabla_v(2,3)+lambda(ix^D,3,3)*nabla_v(3,3))
          end if
        end if
      {end do\}
      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if(ndir==1) then
          w(ix^D,mom(1))=(lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))*qdt+w(ix^D,mom(1))
        else if(ndir==2) then
          if(phi_==2) then
            w(ix^D,mom(1))=((lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))+&
                            (lambda(ix^D,1,1)-lambda(ix^D,2,2))*invr(ix^D))*qdt+w(ix^D,mom(1))
            w(ix^D,mom(2))=((lambda(ix1+1,1,2)-lambda(ix1-1,1,2))/(x(ix1+1,1)-x(ix1-1,1))+&
                            two*lambda(ix^D,1,2)*invr(ix^D))*qdt+w(ix^D,mom(2))
          else
            w(ix^D,mom(1))=((lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))+&
                            (lambda(ix^D,1,1)-lambda(ix^D,2,2))*invr(ix^D))*qdt+w(ix^D,mom(1))
          end if
        else
          if(phi_==2) then
            w(ix^D,mom(1))=((lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))+&
                            (lambda(ix^D,1,1)-lambda(ix^D,2,2))*invr(ix^D))*qdt+w(ix^D,mom(1))
            w(ix^D,mom(2))=((lambda(ix1+1,1,2)-lambda(ix1-1,1,2))/(x(ix1+1,1)-x(ix1-1,1))+&
                            two*lambda(ix^D,1,2)*invr(ix^D))*qdt+w(ix^D,mom(2))
            w(ix^D,mom(3))=((lambda(ix1+1,1,3)-lambda(ix1-1,1,3))/(x(ix1+1,1)-x(ix1-1,1))+&
                             lambda(ix^D,1,3)*invr(ix^D))*qdt+w(ix^D,mom(3))
          else
            w(ix^D,mom(1))=((lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))+&
                            (lambda(ix^D,1,1)-lambda(ix^D,2,2))*invr(ix^D))*qdt+w(ix^D,mom(1))
            w(ix^D,mom(2))=((lambda(ix1+1,1,3)-lambda(ix1-1,1,3))/(x(ix1+1,1)-x(ix1-1,1))+&
                            two*lambda(ix^D,1,3)*invr(ix^D))*qdt+w(ix^D,mom(2))
            w(ix^D,mom(3))=((lambda(ix1+1,1,2)-lambda(ix1-1,1,2))/(x(ix1+1,1)-x(ix1-1,1))+&
                             lambda(ix^D,1,2)*invr(ix^D))*qdt+w(ix^D,mom(3))
          end if
        end if
      {end do\}
      }
      if(total_energy) then
        ! de/dt= +div(v.dot.[mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v *kr])
        ! thus e=e+d_i v_j tensor_ji
        call divvector(vlambda,ixI^L,ixO^L,invr)
        w(ixO^S,e_)=w(ixO^S,e_)+invr(ixO^S)*qdt
      end if

    end if

  end subroutine viscosity_add_source_cylinder

  subroutine viscosity_get_dt(wprim,ixI^L,ixO^L,dtnew,dx^D,x)
    ! Check diffusion time limit for dt < dtdiffpar * dx**2 / (mu/rho)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods, only: usr_set_viscosity

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in) :: wprim(ixI^S,1:nw)
    double precision, intent(inout) :: dtnew

    double precision :: tmp(ixI^S),rho(ixI^S),vc_mus(ixI^S)
    double precision:: dtdiff_visc, dxinv2(1:ndim), max_mu
    integer:: idim

    ! Calculate the kinematic viscosity tmp=mu/rho
    ! here vc_mu must be non-zero!!!
    ! allow for handling of split of densities by calling get_rho
    call phys_get_rho(wprim,x,ixI^L,ixO^L,rho)
    if(associated(usr_set_viscosity)) then
      call usr_set_viscosity(ixI^L,ixO^L,wprim,x,vc_mus)
      tmp(ixO^S)=vc_mus(ixO^S)/rho(ixO^S)
    else
      tmp(ixO^S)=vc_mu/rho(ixO^S)
    end if

    if(slab_uniform)then
      ^D&dxinv2(^D)=one/dx^D**2;
      do idim=1,ndim
         dtdiff_visc=dtdiffpar/maxval(tmp(ixO^S)*dxinv2(idim))
         ! limit the time step
         dtnew=min(dtnew,dtdiff_visc)
      enddo
    else
      do idim=1,ndim
         max_mu=maxval(tmp(ixO^S)/block%ds(ixO^S,idim)**2)
         dtdiff_visc=dtdiffpar/max_mu
         ! limit the time step
         dtnew=min(dtnew,dtdiff_visc)
      enddo
    endif

  end subroutine viscosity_get_dt

  subroutine sub_add_source(qdt,ixI^L,ixO^L,wCT,wp,w,x,&
       energy,qsourcesplit,active)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw), wp(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active
  end subroutine sub_add_source

end module mod_viscosity
