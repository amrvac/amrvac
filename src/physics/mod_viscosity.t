!>   The module add viscous source terms and check time step
!>
!>   Viscous forces in the momentum equations:
!>   d m_i/dt +=  - div (vc_mu * PI)
!>   !! Viscous work in the energy equation:
!>   !! de/dt    += - div (v . vc_mu * PI)
!>   where the PI stress tensor is
!>   PI_i,j = - (dv_j/dx_i + dv_i/dx_j) + (2/3)*Sum_k dv_k/dx_k
!>   where vc_mu is the dynamic viscosity coefficient (g cm^-1 s^-1).
module mod_viscosity
  implicit none

  !> Viscosity coefficient
  double precision, public :: vc_mu = 1.d0

  !> fourth order
  logical :: vc_4th_order = .false.

  !> source split or not
  logical :: vc_split= .false.

  !> whether to compute the viscous terms as
  !> fluxes (ie in the div on the LHS), or not (by default)
  logical :: viscInDiv= .false.

  !> Index of the density (in the w array)
  integer, private, parameter              :: rho_ = 1

  !> Indices of the momentum density
  integer, allocatable, private, protected :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, private, protected              :: e_

  ! Public methods
  public :: visc_get_flux_prim

contains
  !> Read this module"s parameters from a file
  subroutine vc_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /vc_list/ vc_mu, vc_4th_order, vc_split, viscInDiv

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, vc_list, end=111)
111    close(unitpar)
    end do

  end subroutine vc_params_read

  !> Initialize the module
  subroutine viscosity_init(phys_wider_stencil,phys_req_diagonal)
    use mod_global_parameters
    integer, intent(inout) :: phys_wider_stencil
    logical, intent(inout) :: phys_req_diagonal
    integer :: nwx,idir

    call vc_params_read(par_files)

    ! Determine flux variables
    nwx = 1                  ! rho (density)

    allocate(mom(ndir))
    do idir = 1, ndir
       nwx    = nwx + 1
       mom(idir) = nwx       ! momentum density
    end do

    nwx = nwx + 1
    e_     = nwx          ! energy density

    if (viscInDiv) then
      ! to compute the derivatives from left and right upwinded values
      phys_wider_stencil = 1
      phys_req_diagonal = .true.  ! viscInDiv
    end if

  end subroutine viscosity_init

  subroutine viscosity_add_source(qdt,ixI^L,ixO^L,wCT,w,x,&
       energy,qsourcesplit,active)
  ! Add viscosity source in isotropic Newtonian fluids to w within ixO
  ! neglecting bulk viscosity
  ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
    use mod_global_parameters
    use mod_geometry
    use mod_physics, only: phys_solve_eaux

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active

    integer:: ix^L,idim,idir,jdir,iw
    double precision:: lambda(ixI^S,ndir,ndir),tmp(ixI^S),tmp2(ixI^S),v(ixI^S,ndir),vlambda(ixI^S,ndir)

    if (viscInDiv) return

    if(qsourcesplit .eqv. vc_split) then
      active = .true.
      ! standard case, textbook viscosity
      ! Calculating viscosity sources
      if(.not.vc_4th_order) then
        ! involves second derivatives, two extra layers
        ix^L=ixO^L^LADD2;
        if({ ixImin^D>ixmin^D .or. ixImax^D<ixmax^D|.or.})&
          call mpistop("error for viscous source addition, 2 layers needed")
        ix^L=ixO^L^LADD1;
      else
        ! involves second derivatives, four extra layers
        ix^L=ixO^L^LADD4;
        if({ ixImin^D>ixmin^D .or. ixImax^D<ixmax^D|.or.})&
          call mpistop("error for viscous source addition"//&
          "requested fourth order gradients: 4 layers needed")
        ix^L=ixO^L^LADD2;
      end if

      ! get velocity
      do idir=1,ndir
        v(ixI^S,idir)=wCT(ixI^S,mom(idir))/wCT(ixI^S,rho_)
      end do

      ! construct lambda tensor: lambda_ij = gradv_ij + gradv_ji
      ! initialize
      lambda=zero

      !next construct
      do idim=1,ndim; do idir=1,ndir
      ! Calculate velocity gradient tensor within ixL: gradv= grad v,
      ! thus gradv_ij=d_j v_i
        tmp(ixI^S)=v(ixI^S,idir)
        ! Correction for Christoffel terms in non-cartesian
        if (coordinate==cylindrical .and. idim==r_  .and. idir==phi_  ) tmp(ixI^S) = tmp(ixI^S)/x(ixI^S,1)
        if (coordinate==spherical) then
          if     (idim==r_  .and. (idir==2 .or. idir==phi_)) then
            tmp(ixI^S) = tmp(ixI^S)/x(ixI^S,1)
{^NOONED
          elseif (idim==2  .and. idir==phi_) then
            tmp(ixI^S)=tmp(ixI^S)/dsin(x(ixI^S,2))
}
          endif
        endif
        call gradient(tmp,ixI^L,ix^L,idim,tmp2)
        ! Correction for Christoffel terms in non-cartesian
        if (coordinate==cylindrical .and. idim==r_  .and. idir==phi_  ) tmp2(ix^S)=tmp2(ix^S)*x(ix^S,1)
        if (coordinate==cylindrical .and. idim==phi_ .and. idir==phi_ ) tmp2(ix^S)=tmp2(ix^S)+v(ix^S,r_)/x(ix^S,1)
        if (coordinate==spherical) then
          if (idim==r_  .and. (idir==2 .or. idir==phi_)) then
            tmp2(ix^S) = tmp2(ix^S)*x(ix^S,1)
{^NOONED
          elseif (idim==2  .and. idir==phi_ ) then
            tmp2(ix^S)=tmp2(ix^S)*dsin(x(ix^S,2))
          elseif (idim==2   .and. idir==2   ) then
            tmp2(ix^S)=tmp2(ix^S)+v(ix^S,r_)/x(ix^S,1)
          elseif (idim==phi_.and. idir==phi_) then
            tmp2(ix^S)=tmp2(ix^S)+v(ix^S,r_)/x(ix^S,1)+v(ix^S,2)/(x(ix^S,1)*dtan(x(ix^S,2)))
}
          endif
        endif
        lambda(ix^S,idim,idir)= lambda(ix^S,idim,idir)+ tmp2(ix^S)
        lambda(ix^S,idir,idim)= lambda(ix^S,idir,idim)+ tmp2(ix^S)
      enddo; enddo;

      ! Multiply lambda with viscosity coefficient and dt
      lambda(ix^S,1:ndir,1:ndir)=lambda(ix^S,1:ndir,1:ndir)*vc_mu*qdt

      !calculate div v term through trace action separately
      ! rq : it is safe to use the trace rather than compute the divergence
      !      since we always retrieve the divergence (even with the
      !      Christoffel terms)
      tmp=0.d0
      do idir=1,ndir
         tmp(ix^S)=tmp(ix^S)+lambda(ix^S,idir,idir)
      end do
      tmp(ix^S)=tmp(ix^S)/3.d0

      !substract trace from diagonal elements
      do idir=1,ndir
         lambda(ix^S,idir,idir)=lambda(ix^S,idir,idir)-tmp(ix^S)
      enddo

      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
      ! hence m_j=m_j+d_i tensor_ji
      do idir=1,ndir
        do idim=1,ndim
              tmp(ix^S)=lambda(ix^S,idir,idim)
              ! Correction for divergence of a tensor
              if (coordinate==cylindrical .and. idim==r_ .and. (idir==r_ .or. idir==z_)) tmp(ix^S) = tmp(ix^S)*x(ix^S,1)
              if (coordinate==cylindrical .and. idim==r_ .and. idir==phi_              ) tmp(ix^S) = tmp(ix^S)*x(ix^S,1)**two
              if (coordinate==spherical) then
                if (idim==r_ .and. idir==r_                 ) tmp(ix^S) = tmp(ix^S)*x(ix^S,1)**two
                if (idim==r_ .and. (idir==2 .or. idir==phi_)) tmp(ix^S) = tmp(ix^S)*x(ix^S,1)**3.d0
{^NOONED
                if (idim==2  .and. (idir==r_ .or. idir==2))   tmp(ix^S) = tmp(ix^S)*dsin(x(ix^S,2))
                if (idim==2  .and. idir==phi_               ) tmp(ix^S) = tmp(ix^S)*dsin(x(ix^S,2))**two
}
              endif
              call gradient(tmp,ixI^L,ixO^L,idim,tmp2)
              ! Correction for divergence of a tensor
              if (coordinate==cylindrical .and. idim==r_ .and. (idir==r_ .or. idir==z_)) tmp2(ixO^S) = tmp2(ixO^S)/x(ixO^S,1)
              if (coordinate==cylindrical .and. idim==r_ .and. idir==phi_              ) tmp2(ixO^S) = tmp2(ixO^S)/(x(ixO^S,1)**two)
              if (coordinate==spherical) then
                if (idim==r_ .and. idir==r_                 ) tmp2(ixO^S) = tmp2(ixO^S)/(x(ixO^S,1)**two)
                if (idim==r_ .and. (idir==2 .or. idir==phi_)) tmp2(ixO^S) = tmp2(ixO^S)/(x(ixO^S,1)**3.d0)
{^NOONED
                if (idim==2  .and. (idir==r_ .or. idir==2))   tmp2(ixO^S) = tmp2(ixO^S)/(dsin(x(ixO^S,2)))
                if (idim==2  .and. idir==phi_               ) tmp2(ixO^S) = tmp2(ixO^S)/(dsin(x(ixO^S,2))**two)
}
              endif
              w(ixO^S,mom(idir))=w(ixO^S,mom(idir))+tmp2(ixO^S)
        enddo
        ! Correction for geometrical terms in the div of a tensor
        if (coordinate==cylindrical .and. idir==r_  ) w(ixO^S,mom(idir))=w(ixO^S,mom(idir))-lambda(ixO^S,phi_,phi_)/x(ixO^S,1)
        if (coordinate==spherical   .and. idir==r_  ) w(ixO^S,mom(idir))=w(ixO^S,mom(idir))-(lambda(ixO^S,2,2)+lambda(ixO^S,phi_,phi_))/x(ixO^S,1)
{^NOONED
        if (coordinate==spherical   .and. idir==2   ) w(ixO^S,mom(idir))=w(ixO^S,mom(idir))-lambda(ixO^S,phi_,phi_)/(x(ixO^S,1)/dtan(x(ixO^S,2)))
}
      end do

      if(energy) then
        ! de/dt= +div(v.dot.[mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v *kr])
        ! thus e=e+d_i v_j tensor_ji
        vlambda=0.d0
        do idim=1,ndim
          do idir=1,ndir
             vlambda(ixI^S,idim)=vlambda(ixI^S,idim)+v(ixI^S,idir)*lambda(ixI^S,idir,idim)
          end do
        enddo
        call divvector(vlambda,ixI^L,ixO^L,tmp2)
        w(ixO^S,e_)=w(ixO^S,e_)+tmp2(ixO^S)
        if(phys_solve_eaux) w(ixO^S,iw_eaux)=w(ixO^S,iw_eaux)+tmp2(ixO^S)
      end if
    end if

  end subroutine viscosity_add_source

  subroutine viscosity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    ! Check diffusion time limit for dt < dtdiffpar * dx**2 / (mu/rho)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision, intent(inout) :: dtnew

    double precision :: tmp(ixI^S)
    double precision:: dtdiff_visc, dxinv2(1:ndim)
    integer:: idim

    ! Calculate the kinematic viscosity tmp=mu/rho

    tmp(ixO^S)=vc_mu/w(ixO^S,rho_)

    ^D&dxinv2(^D)=one/dx^D**2;
    do idim=1,ndim
       dtdiff_visc=dtdiffpar/maxval(tmp(ixO^S)*dxinv2(idim))
       ! limit the time step
       dtnew=min(dtnew,dtdiff_visc)
    enddo

  end subroutine viscosity_get_dt

  ! viscInDiv
  ! Get the viscous stress tensor terms in the idim direction
  ! Beware : a priori, won't work for ndir /= ndim
  ! Rq : we work with primitive w variables here
  ! Rq : ixO^L is already extended by 1 unit in the direction we work on

  subroutine visc_get_flux_prim(w, x, ixI^L, ixO^L, idim, f, energy)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:^ND)
    double precision, intent(inout) :: f(ixI^S, nwflux)
    logical, intent(in) :: energy
    integer                         :: idir, i
    double precision :: v(ixI^S,1:ndir)

    double precision                :: divergence(ixI^S)

    double precision:: lambda(ixI^S,ndir) !, tmp(ixI^S) !gradV(ixI^S,ndir,ndir)

    if (.not. viscInDiv) return

    do i=1,ndir
     v(ixI^S,i)=w(ixI^S,i+1)
    enddo
    call divvector(v,ixI^L,ixO^L,divergence)

    call get_crossgrad(ixI^L,ixO^L,x,w,idim,lambda)
    lambda(ixO^S,idim) = lambda(ixO^S,idim) - (2.d0/3.d0) * divergence(ixO^S)

    ! Compute the idim-th row of the viscous stress tensor
    do idir = 1, ndir
      f(ixO^S, mom(idir)) = f(ixO^S, mom(idir)) - vc_mu * lambda(ixO^S,idir)
      if (energy) f(ixO^S, e_) = f(ixO^S, e_) - vc_mu * lambda(ixO^S,idir) * v(ixI^S,idir)
    enddo

  end subroutine visc_get_flux_prim

  ! Compute the cross term ( d_i v_j + d_j v_i in Cartesian BUT NOT IN
  ! CYLINDRICAL AND SPHERICAL )
  subroutine get_crossgrad(ixI^L,ixO^L,x,w,idim,cross)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(out)   :: cross(ixI^S,ndir)
    integer :: idir
    double precision :: tmp(ixI^S), v(ixI^S)

    if (ndir/=ndim) call mpistop("This formula are probably wrong for ndim/=ndir")
    ! Beware also, we work w/ the angle as the 3rd component in cylindrical
    ! and the colatitude as the 2nd one in spherical
    cross(ixI^S,:)=zero
    tmp(ixI^S)=zero
    select case(coordinate)
    case (Cartesian,Cartesian_stretched)
      call cart_cross_grad(ixI^L,ixO^L,x,w,idim,cross)
    case (cylindrical)
      if (idim==1) then
        ! for rr and rz
        call cart_cross_grad(ixI^L,ixO^L,x,w,idim,cross)
        ! then we overwrite rth w/ the correct expression
        {^NOONED
        v(ixI^S)=w(ixI^S,mom(1)) ! v_r
        call gradient(v,ixI^L,ixO^L,2,tmp) ! d_th (rq : already contains 1/r)
        cross(ixI^S,2)=tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(2))/x(ixI^S,1)  ! v_th / r
        call gradient(v,ixI^L,ixO^L,1,tmp) ! d_r
        cross(ixI^S,2)=cross(ixI^S,2)+tmp(ixI^S)*x(ixI^S,1)
        }
      elseif (idim==2) then
        ! thr (idem as above)
        v(ixI^S)=w(ixI^S,mom(1)) ! v_r
        {^NOONED
        call gradient(v,ixI^L,ixO^L,2,tmp) ! d_th
        cross(ixI^S,1)=tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(2))/x(ixI^S,1)  ! v_th / r
        call gradient(v,ixI^L,ixO^L,1,tmp) ! d_r
        cross(ixI^S,1)=cross(ixI^S,1)+tmp(ixI^S)*x(ixI^S,1)
        ! thth
        v(ixI^S)=w(ixI^S,mom(2)) ! v_th
        call gradient(v,ixI^L,ixO^L,2,tmp) ! d_th
        cross(ixI^S,2)=two*tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(1)) ! v_r
        cross(ixI^S,2)=cross(ixI^S,2)+two*v(ixI^S)/x(ixI^S,1) ! + 2 vr/r
        !thz
        v(ixI^S)=w(ixI^S,mom(3)) ! v_z
        call gradient(v,ixI^L,ixO^L,2,tmp) ! d_th
        }
        cross(ixI^S,3)=tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(2))  ! v_th
        {^IFTHREED
        call gradient(v,ixI^L,ixO^L,3,tmp) ! d_z
        }
        cross(ixI^S,3)=cross(ixI^S,3)+tmp(ixI^S)
        {^IFTHREED
      elseif (idim==3) then
        ! for zz and rz
        call cart_cross_grad(ixI^L,ixO^L,x,w,idim,cross)
        ! then we overwrite zth w/ the correct expression
        !thz
        v(ixI^S)=w(ixI^S,mom(3)) ! v_z
        call gradient(v,ixI^L,ixO^L,2,tmp) ! d_th
        cross(ixI^S,2)=tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(2))  ! v_th
        call gradient(v,ixI^L,ixO^L,3,tmp) ! d_z
        cross(ixI^S,2)=cross(ixI^S,2)+tmp(ixI^S)
        }
      endif
    case (spherical)
      if (idim==1) then
        ! rr (normal, simply 2 * dr vr)
        v(ixI^S)=w(ixI^S,mom(1)) ! v_r
        call gradient(v,ixI^L,ixO^L,1,tmp) ! d_r
        cross(ixI^S,1)=two*tmp(ixI^S)
        !rth
        v(ixI^S)=w(ixI^S,mom(1)) ! v_r
        {^NOONED
        call gradient(v,ixI^L,ixO^L,2,tmp) ! d_th (rq : already contains 1/r)
        cross(ixI^S,2)=tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(2))/x(ixI^S,1)  ! v_th / r
        call gradient(v,ixI^L,ixO^L,1,tmp) ! d_r
        cross(ixI^S,2)=cross(ixI^S,2)+tmp(ixI^S)*x(ixI^S,1)
        }
        {^IFTHREED
        !rph
        v(ixI^S)=w(ixI^S,mom(1)) ! v_r
        call gradient(v,ixI^L,ixO^L,3,tmp) ! d_phi (rq : contains 1/rsin(th))
        cross(ixI^S,3)=tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(3))/x(ixI^S,1) ! v_phi / r
        call gradient(v,ixI^L,ixO^L,1,tmp) ! d_r
        cross(ixI^S,3)=cross(ixI^S,3)+tmp(ixI^S)*x(ixI^S,1)
        }
      elseif (idim==2) then
        ! thr
        v(ixI^S)=w(ixI^S,mom(1)) ! v_r
        {^NOONED
        call gradient(v,ixI^L,ixO^L,2,tmp) ! d_th (rq : already contains 1/r)
        cross(ixI^S,1)=tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(2))/x(ixI^S,1)  ! v_th / r
        call gradient(v,ixI^L,ixO^L,1,tmp) ! d_r
        cross(ixI^S,1)=cross(ixI^S,1)+tmp(ixI^S)*x(ixI^S,1)
        ! thth
        v(ixI^S)=w(ixI^S,mom(2)) ! v_th
        call gradient(v,ixI^L,ixO^L,2,tmp) ! d_th
        cross(ixI^S,2)=two*tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(1)) ! v_r
        cross(ixI^S,2)=cross(ixI^S,2)+two*v(ixI^S)/x(ixI^S,1) ! + 2 vr/r
        }
        {^IFTHREED
        !thph
        v(ixI^S)=w(ixI^S,mom(2)) ! v_th
        call gradient(v,ixI^L,ixO^L,3,tmp) ! d_phi (rq : contains 1/rsin(th))
        cross(ixI^S,3)=tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(3))/dsin(x(ixI^S,2)) ! v_ph / sin(th)
        call gradient(v,ixI^L,ixO^L,2,tmp) ! d_th
        cross(ixI^S,3)=cross(ixI^S,3)+tmp(ixI^S)*dsin(x(ixI^S,2))
        }
        {^IFTHREED
      elseif (idim==3) then
        !phr
        v(ixI^S)=w(ixI^S,mom(1)) ! v_r
        call gradient(v,ixI^L,ixO^L,3,tmp) ! d_phi
        cross(ixI^S,1)=tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(3))/x(ixI^S,1) ! v_phi / r
        call gradient(v,ixI^L,ixO^L,1,tmp) ! d_r
        cross(ixI^S,1)=cross(ixI^S,1)+tmp(ixI^S)*x(ixI^S,1)
        !phth
        v(ixI^S)=w(ixI^S,mom(2)) ! v_th
        call gradient(v,ixI^L,ixO^L,3,tmp) ! d_phi
        cross(ixI^S,2)=tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(3))/dsin(x(ixI^S,2)) ! v_ph / sin(th)
        call gradient(v,ixI^L,ixO^L,2,tmp) ! d_th
        cross(ixI^S,2)=cross(ixI^S,2)+tmp(ixI^S)*dsin(x(ixI^S,2))
        !phph
        v(ixI^S)=w(ixI^S,mom(3)) ! v_ph
        call gradient(v,ixI^L,ixO^L,3,tmp) ! d_phi
        cross(ixI^S,3)=two*tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(1)) ! v_r
        cross(ixI^S,3)=cross(ixI^S,3)+two*v(ixI^S)/x(ixI^S,1) ! + 2 vr/r
        v(ixI^S)=w(ixI^S,mom(2)) ! v_th
        cross(ixI^S,3)=cross(ixI^S,3)+two*v(ixI^S)/(x(ixI^S,1)*dtan(x(ixI^S,2))) ! + 2 vth/(rtan(th))
        }
      endif
    case default
      call mpistop("Unknown geometry specified")
    end select

  end subroutine get_crossgrad

  !> yields d_i v_j + d_j v_i for a given i, OK in Cartesian and for some
  !> tensor terms in cylindrical (rr & rz) and in spherical (rr)
  subroutine cart_cross_grad(ixI^L,ixO^L,x,w,idim,cross)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:^ND)
    double precision, intent(out)   :: cross(ixI^S,ndir)
    integer :: idir
    double precision :: tmp(ixI^S), v(ixI^S)

    v(ixI^S)=w(ixI^S,mom(idim))
    do idir=1,ndir
      call gradient(v,ixI^L,ixO^L,idir,tmp)
      cross(ixO^S,idir)=tmp(ixO^S)
    enddo
    do idir=1,ndir
      v(ixI^S)=w(ixI^S,mom(idir))
      call gradient(v,ixI^L,ixO^L,idim,tmp)
      cross(ixO^S,idir)=cross(ixO^S,idir)+tmp(ixO^S)
    enddo

  end subroutine cart_cross_grad

  subroutine visc_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)
    use mod_global_parameters
    use mod_geometry
    ! w and wCT conservative variables here
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S, 1:ndim)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
    ! to change and to set as a parameter in the parfile once the possibility to
    ! solve the equations in an angular momentum conserving form has been
    ! implemented (change tvdlf.t eg)
    double precision :: v(ixI^S,1:ndir), vv(ixI^S), divergence(ixI^S)
    double precision :: tmp(ixI^S),tmp1(ixI^S)
    integer          :: i

    if (.not. viscInDiv) return

    select case (coordinate)
    case (cylindrical)
      ! get the velocity components
      do i=1,ndir
       v(ixI^S,i)=wCT(ixI^S,mom(i))/wCT(ixI^S,rho_)
      enddo
      ! thth tensor term - - -
        ! 1st the cross grad term
{^NOONED
      vv(ixI^S)=v(ixI^S,2) ! v_th
      call gradient(vv,ixI^L,ixO^L,2,tmp1) ! d_th
      tmp(ixO^S)=two*(tmp1(ixO^S)+v(ixI^S,1)/x(ixO^S,1)) ! 2 ( d_th v_th / r + vr/r )
        ! 2nd the divergence
      call divvector(v,ixI^L,ixO^L,divergence)
      tmp(ixO^S) = tmp(ixO^S) - (2.d0/3.d0) * divergence(ixO^S)
      ! s[mr]=-thth/radius
      w(ixO^S,mom(1))=w(ixO^S,mom(1))-qdt*vc_mu*tmp(ixO^S)/x(ixO^S,1)
      if (.not. angmomfix) then
        ! rth tensor term - - -
        vv(ixI^S)=v(ixI^S,1) ! v_r
        call gradient(vv,ixI^L,ixO^L,2,tmp1) ! d_th
        tmp(ixO^S)=tmp1(ixO^S)
        vv(ixI^S)=v(ixI^S,2)/x(ixI^S,1)  ! v_th / r
        call gradient(vv,ixI^L,ixO^L,1,tmp1) ! d_r
        tmp(ixO^S)=tmp(ixO^S)+tmp1(ixO^S)*x(ixO^S,1)
        ! s[mphi]=+rth/radius
        w(ixO^S,mom(2))=w(ixO^S,mom(2))+qdt*vc_mu*tmp(ixO^S)/x(ixO^S,1)
      endif
}
    case (spherical)
      ! get the velocity components
      do i=1,ndir
       v(ixI^S,i)=wCT(ixI^S,mom(i))/wCT(ixI^S,rho_)
      enddo
      ! thth tensor term - - -
      ! 1st the cross grad term
      vv(ixI^S)=v(ixI^S,2) ! v_th
{^NOONED
      call gradient(vv,ixI^L,ixO^L,2,tmp1) ! d_th
      tmp(ixO^S)=two*(tmp1(ixO^S)+v(ixO^S,1)/x(ixO^S,1)) ! 2 ( 1/r * d_th v_th + vr/r )
      ! 2nd the divergence
      call divvector(v,ixI^L,ixO^L,divergence)
      tmp(ixO^S) = tmp(ixO^S) - (2.d0/3.d0) * divergence(ixO^S)
      ! s[mr]=-thth/radius
      w(ixO^S,mom(1))=w(ixO^S,mom(1))-qdt*vc_mu*tmp(ixO^S)/x(ixO^S,1)
}
      ! phiphi tensor term - - -
      ! 1st the cross grad term
      vv(ixI^S)=v(ixI^S,3) ! v_ph
{^IFTHREED
      call gradient(vv,ixI^L,ixO^L,3,tmp1) ! d_phi
      tmp(ixO^S)=two*(tmp1(ixO^S)+v(ixO^S,1)/x(ixO^S,1)+v(ixO^S,2)/(x(ixO^S,1)*dtan(x(ixO^S,2)))) ! 2 ( 1/rsinth * d_ph v_ph + vr/r + vth/rtanth )
}
      ! 2nd the divergence
      tmp(ixO^S) = tmp(ixO^S) - (2.d0/3.d0) * divergence(ixO^S)
      ! s[mr]=-phiphi/radius
      w(ixO^S,mom(1))=w(ixO^S,mom(1))-qdt*vc_mu*tmp(ixO^S)/x(ixO^S,1)
      ! s[mth]=-cotanth*phiphi/radius
{^NOONED
      w(ixO^S,mom(2))=w(ixO^S,mom(2))-qdt*vc_mu*tmp(ixO^S)/(x(ixO^S,1)*dtan(x(ixO^S,2)))
}
      if (.not. angmomfix) then
        ! rth tensor term - - -
        vv(ixI^S)=v(ixI^S,1) ! v_r
        call gradient(vv,ixI^L,ixO^L,2,tmp) ! d_th (rq : already contains 1/r)
        vv(ixI^S)=v(ixI^S,2)/x(ixI^S,1)  ! v_th / r
        call gradient(vv,ixI^L,ixO^L,1,tmp1) ! d_r
        tmp(ixO^S)=tmp(ixO^S)+tmp1(ixO^S)*x(ixO^S,1)
        ! s[mth]=+rth/radius
        w(ixO^S,mom(2))=w(ixO^S,mom(2))+qdt*vc_mu*tmp(ixO^S)/x(ixO^S,1)
        ! rphi tensor term - - -
        vv(ixI^S)=v(ixI^S,1) ! v_r
{^IFTHREED
        call gradient(vv,ixI^L,ixO^L,3,tmp) ! d_phi (rq : contains 1/rsin(th))
}
        vv(ixI^S)=v(ixI^S,3)/x(ixI^S,1) ! v_phi / r
        call gradient(vv,ixI^L,ixO^L,1,tmp1) ! d_r
        tmp(ixO^S)=tmp(ixO^S)+tmp1(ixO^S)*x(ixO^S,1)
        ! s[mphi]=+rphi/radius
        w(ixO^S,mom(3))=w(ixO^S,mom(3))+qdt*vc_mu*tmp(ixO^S)/x(ixO^S,1)
        ! phith tensor term - - -
        vv(ixI^S)=v(ixI^S,2) ! v_th
{^IFTHREED
        call gradient(vv,ixI^L,ixO^L,3,tmp) ! d_phi
}
{^NOONED
        vv(ixI^S)=v(ixI^S,3)/dsin(x(ixI^S,2)) ! v_ph / sin(th)
        call gradient(vv,ixI^L,ixO^L,2,tmp1) ! d_th
        tmp(ixO^S)=tmp(ixO^S)+tmp1(ixO^S)*dsin(x(ixO^S,2))
        ! s[mphi]=+cotanth*phith/radius
        w(ixO^S,mom(3))=w(ixO^S,mom(3))+qdt*vc_mu*tmp(ixO^S)/(x(ixO^S,1)*dtan(x(ixO^S,2)))
}
      endif
    end select

  end subroutine visc_add_source_geom

end module mod_viscosity
