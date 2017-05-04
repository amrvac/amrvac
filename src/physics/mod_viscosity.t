!>   THE MODULE ADD VISCOUS SOURCE TERMS AND CHECK TIME STEP
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

  !> The adiabatic index
  double precision, private, protected :: vc_gamma

  !> fourth order
  logical :: vc_4th_order = .false.

  !> source split or not
  logical :: vc_split= .false.

  !> Index of the density (in the w array)
  integer, private, parameter              :: rho_ = 1

  !> Indices of the momentum density
  integer, allocatable, private, protected :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, private, protected              :: e_

contains
  !> Read this module"s parameters from a file
  subroutine vc_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /vc_list/ vc_mu, vc_4th_order, vc_split

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, vc_list, end=111)
111    close(unitpar)
    end do

  end subroutine vc_params_read

  !> Initialize the module
  subroutine viscosity_init(phys_gamma)
    use mod_global_parameters
    integer :: nwx,idir
    double precision, intent(in) :: phys_gamma

    vc_gamma=phys_gamma

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

  end subroutine viscosity_init

  subroutine viscosity_add_source(qdt,ixI^L,ixO^L,wCT,w,x,energy,qsourcesplit)
  ! Add viscosity source to w within ixO 
    use mod_global_parameters
    
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    
    integer:: ix^L,idim,idir,jdir,iw
    double precision:: lambda(ixI^S,ndir,ndir),tmp(ixI^S),tmp2(ixI^S),v(ixI^S,ndir)

    if(qsourcesplit .eqv. vc_split) then
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
          call mpistop("error for viscous source addition, &
          requested fourth order gradients: 4 layers needed")
        ix^L=ixO^L^LADD2;
      end if

      ! get velocity
      do idir=1,ndir
        v(ixI^S,idir)=wCT(ixI^S,mom(idir))/wCT(ixI^S,rho_)
      end do

      ! construct lambda tensor: lambda_ij = gradv_ij + gradv_ji
      ! initialize
      lambda(ix^S,1:ndir,1:ndir)=zero
      
      !next construct 
      do idim=1,ndim; do idir=1,ndir
      ! Calculate velocity gradient tensor within ixL: gradv= grad v, 
      ! thus gradv_ij=d_j v_i
        tmp(ixI^S)=v(ixI^S,idir)
        select case(typegrad)
        case("central")
          call gradient(tmp,ixI^L,ix^L,idim,tmp2)
        case("limited")
          call gradientS(tmp,ixI^L,ix^L,idim,tmp2)
        end select
        lambda(ix^S,idim,idir)= lambda(ix^S,idim,idir)+ tmp2(ix^S)
        lambda(ix^S,idir,idim)= lambda(ix^S,idir,idim)+ tmp2(ix^S)
      enddo; enddo;
      
      ! Multiply lambda with viscosity coefficient and dt
      lambda=lambda*vc_mu*qdt
      
      !calculate div v term through trace action separately
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
              select case(typegrad)
              case("central")
                call gradient(tmp,ixI^L,ixO^L,idim,tmp2)
              case("limited")
                call gradientS(tmp,ixI^L,ixO^L,idim,tmp2)
              end select
              w(ixO^S,mom(idir))=w(ixO^S,mom(idir))+tmp2(ixO^S)
        enddo
      end do

      if(energy) then
        ! de/dt= +div(v.dot.[mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v *kr])
        ! thus e=e+d_i v_j tensor_ji
        do idim=1,ndim
          tmp=0.d0
          do idir=1,ndir
             tmp(ix^S)=tmp(ix^S)+v(ix^S,idir)*lambda(ix^S,idir,idim)
          end do
          call gradient(tmp,ixI^L,ixO^L,idim,tmp2)
          if (solve_pthermal) then
            w(ixO^S,e_)=w(ixO^S,e_)+tmp2(ixO^S)*(vc_gamma-1.d0)
          else
            w(ixO^S,e_)=w(ixO^S,e_)+tmp2(ixO^S)
          endif
        enddo
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

end module mod_viscosity
