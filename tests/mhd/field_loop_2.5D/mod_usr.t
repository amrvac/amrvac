! Test of an advecting field loop, copied from: "An unsplit Godunov method for
! ideal MHD via constrained transport", Gardiner et al. 2005
! (http://dx.doi.org/10.1016/j.jcp.2004.11.016).
module mod_usr
  use mod_mhd

  implicit none

  ! Amplitude of vector potential
  double precision, parameter :: A0 = 1.0d-3

  ! Radius of field loop
  double precision, parameter :: R0 = 0.3d0

  ! Initialize as numerical gradient from a vector potential
  logical, parameter :: init_from_vectorpot = .true.

  integer :: i_divb_2, i_divb_4, i_B_err

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods
    use mod_physics

    usr_init_one_grid => initonegrid_usr
    usr_modify_output => set_output_vars
    usr_refine_grid   => my_refine
    usr_init_vector_potential=>initvecpot_usr

    call set_coordinate_system('Cartesian_2.5D')
    call mhd_activate()

    i_divb_2 = var_set_extravar("divb-2", "divb-2")
    i_divb_4 = var_set_extravar("divb-4", "divb-4")
    i_B_err = var_set_extravar("B_err", "B_err")
  end subroutine usr_init

  ! initialize one grid
  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: v0(ndim)
    double precision                :: bfield(ixI^S, ndir)

    select case (iprob)
    case (1)
       v0(:) = [2.0d0, 1.0d0]
    case (2)
       v0(:) = [0.0d0, 0.0d0]
    case default
       call mpistop("Invalid iprob")
    end select

    w(ixO^S,rho_)   = 1.0d0 ! Density
    w(ixO^S,mom(1)) = v0(1) ! Vx
    w(ixO^S,mom(2)) = v0(2) ! Vy
    w(ixO^S,e_)     = 1.0d0 ! Pressure
    if(stagger_grid) then
      call b_from_vector_potential(ixGs^LL,ixI^L,ixO^L,block%ws,x)
      call mhd_face_to_center(ixO^L,block)
    else 
      bfield=0.d0
      call bfield_solution(ixI^L, ixO^L, x, v0, 0.0d0, bfield)
      w(ixO^S, mag(:)) = bfield(ixO^S, :)
    end if
    call mhd_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine initonegrid_usr

  subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    integer, intent(in)                :: ixI^L, ixC^L,idir
    double precision, intent(in)       :: xC(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)

    double precision :: rx(ixC^S)

    A(ixC^S) = 0.d0
    if (idir==3) then
      rx(ixC^S)=sqrt(xC(ixC^S,1)**2+xC(ixC^S,2)**2)
      where(rx(ixC^S)<=R0)
        A(ixC^S) = A0*(R0-rx(ixC^S))
      end where
    end if

  end subroutine initvecpot_usr

  subroutine bfield_solution(ixI^L, ixO^L, x, v, t, bfield)
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: v(ndim), t
    double precision, intent(inout) :: bfield(ixI^S,ndir)
    double precision                :: x1(ixO^S), x2(ixO^S)
    double precision                :: Az(ixI^S)

    if (init_from_vectorpot) then
       call vectorpot_z(ixI^L, ixI^L, x, v, t, Az)
       call gradient(Az, ixI^L, ixO^L, 2, bfield(ixI^S, 1))
       call gradient(Az, ixI^L, ixO^L, 1, bfield(ixI^S, 2))
       bfield(ixO^S, 2) = -bfield(ixO^S, 2)
    else
       ! Determine coordinates modulo domain size
       x1 = x(ixO^S,1) - v(1) * t - xprobmin1
       x1 = xprobmin1 + modulo(x1, xprobmax1-xprobmin1)
       x2 = x(ixO^S,2) - v(2) * t - xprobmin2
       x2 = xprobmin2 + modulo(x2, xprobmax2-xprobmin2)
      where (x1**2 + x2**2 < R0**2)
        bfield(ixO^S, 1) = -A0 * x2/sqrt(x1**2 + x2**2)
        bfield(ixO^S, 2) = A0 * x1/sqrt(x1**2 + x2**2)
      elsewhere
        bfield(ixO^S, 1) = 0.0d0
        bfield(ixO^S, 2) = 0.0d0
      end where
    end if
  end subroutine bfield_solution

  subroutine vectorpot_z(ixI^L, ixO^L, x, v, t, Az)
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: v(ndim), t
    double precision, intent(inout) :: Az(ixO^S)
    double precision                :: x1(ixO^S), x2(ixO^S)

    ! Determine coordinates modulo domain size
    x1 = x(ixO^S,1) - v(1) * t - xprobmin1
    x1 = xprobmin1 + modulo(x1, xprobmax1-xprobmin1)
    x2 = x(ixO^S,2) - v(2) * t - xprobmin2
    x2 = xprobmin2 + modulo(x2, xprobmax2-xprobmin2)

    where (x1**2 + x2**2 < R0**2)
       Az(ixO^S) = A0 * (R0 - sqrt(x1**2 + x2**2))
    elsewhere
       Az(ixO^S) = 0.0d0
    end where
  end subroutine vectorpot_z

  subroutine set_output_vars(ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,nw)
    double precision                :: divb(ixI^S)
    double precision                :: v0(ndim)
    double precision                :: bfield(ixI^S, ndir)
    double precision                :: ws(ixGs^T, nws)

    integer :: idim,hxO^L

    select case (iprob)
    case (1)
       v0(:) = [2.0d0, 1.0d0]
    case (2)
       v0(:) = [0.0d0, 0.0d0]
    case default
       call mpistop("Invalid iprob")
    end select

    ! output divB second order
    call get_divb(w,ixI^L,ixO^L,divb,1)
    w(ixO^S,i_divb_2)=divb(ixO^S)

    ! output divB fourth order
    call get_divb(w,ixI^L,ixO^L,divb,2)
    w(ixO^S,i_divb_4)=divb(ixO^S)

    if(stagger_grid) then
      call b_from_vector_potential(ixGs^LL,ixI^L,ixO^L,ws,x)
      do idim=1,ndim
        hxO^L=ixO^L-kr(idim,^D);
        bfield(ixO^S,idim)=half*(ws(ixO^S,idim)+ws(hxO^S,idim))
      end do
    else
      call bfield_solution(ixI^L, ixO^L, x, v0, qt, bfield)
    end if

    w(ixO^S,i_B_err) = sqrt( &
         (w(ixO^S, mag(1)) - bfield(ixO^S, 1))**2 + &
         (w(ixO^S, mag(2)) - bfield(ixO^S, 2))**2)
  end subroutine set_output_vars

  ! Refine left half of the domain, to test divB methods with refinement
  subroutine my_refine(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters
    integer, intent(in)          :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout)       :: refine, coarsen

    refine = 0
    coarsen = 0
    if (any(x(ixO^S, 1) < 0.0d0)) then
       refine = 1
       coarsen = -1
    else
       refine = -1
       coarsen = -1
    end if
  end subroutine my_refine

end module mod_usr
