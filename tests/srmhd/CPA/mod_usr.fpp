module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  ! Circularly polarized Alfven wave (CPA-3D), special relativistic MHD.
  ! Setup of Del Zanna et al. 2024, Fluids 9, 16, Sect. 4 (exact solution from
  ! Del Zanna, Zanotti, Bucciantini & Londrillo 2007, A&A 473, 11), with the wave
  ! propagating along the main diagonal of the periodic cube as in Mignone et al.
  ! The wave is an exact nonlinear solution: after one period T = 2 pi/(k v_A)
  ! the solution must coincide with the initial condition.
  double precision, parameter :: cpa_rho = 1.0d0   ! rest-mass density
  double precision, parameter :: cpa_p   = 1.0d0   ! thermal pressure
  double precision, parameter :: cpa_B0  = 1.0d0   ! background field along X
  double precision, parameter :: cpa_eta = 1.0d0   ! wave amplitude (B_perp/B0)

  ! Orthonormal basis of the wave frame: e_X along the main diagonal
  ! (rotation angles tan(alpha)=1, sin(beta)=1/sqrt(3), cf. Mignone et al.)
  double precision, parameter :: eX(3) = [ 0.57735026918962584d0, 0.57735026918962584d0, 0.57735026918962584d0]
  double precision, parameter :: eY(3) = [-0.70710678118654746d0, 0.70710678118654746d0, 0.0d0]
  double precision, parameter :: eZ(3) = [-0.40824829046386307d0,-0.40824829046386307d0, 0.81649658092772615d0]

  ! Error norms are computed for these primitive variables
  integer, parameter :: nerr = 8
  character(len=4), parameter :: errnames(nerr) = &
       ['rho ', 'p   ', 'vx  ', 'vy  ', 'vz  ', 'Bx  ', 'By  ', 'Bz  ']

contains

  subroutine usr_init()

    usr_init_one_grid => initonegrid_usr
    usr_print_log     => cpa_print_errors   ! active with typefilelog='special'

    call phys_activate()

  end subroutine usr_init

  !> Alfven speed of the exact relativistic CPA wave (Del Zanna et al. 2007)
  double precision function cpa_vA()
    double precision :: rhoh, denom
    rhoh  = cpa_rho + srmhd_gamma/(srmhd_gamma-1.0d0)*cpa_p
    denom = rhoh + cpa_B0**2*(1.0d0+cpa_eta**2)
    cpa_vA = sqrt( cpa_B0**2/denom &
         / (0.5d0*(1.0d0 + sqrt(1.0d0 - (2.0d0*cpa_eta*cpa_B0**2/denom)**2))) )
  end function cpa_vA

  !> Exact solution at position xx and time t:
  !> prim = (rho, p, vx, vy, vz, Bx, By, Bz), v = 3-velocity.
  !> The wave travels in +X with phase velocity v_A.
  subroutine cpa_exact(xx, t, vA, prim)
    double precision, intent(in)  :: xx(3), t, vA
    double precision, intent(out) :: prim(nerr)
    double precision :: k(3), phase, BY, BZ
    integer :: idir

    ! k_x = k_y = k_z = 2 pi / L  (|k| = 2 pi sqrt(3) for L = 1)
    k(1) = 2.0d0*dpi/(xprobmax1-xprobmin1)
    k(2) = 2.0d0*dpi/(xprobmax2-xprobmin2)
    k(3) = 2.0d0*dpi/(xprobmax3-xprobmin3)

    phase = sum(k*xx) - sqrt(sum(k**2))*vA*t

    ! wave frame: B_X = B0, B_Y = eta B0 cos, B_Z = eta B0 sin,
    !             v_X = 0,  v_{Y,Z} = -v_A B_{Y,Z}/B0
    BY = cpa_eta*cpa_B0*cos(phase)
    BZ = cpa_eta*cpa_B0*sin(phase)

    prim(1) = cpa_rho
    prim(2) = cpa_p
    do idir = 1, 3
       prim(2+idir) = -vA/cpa_B0*(BY*eY(idir) + BZ*eZ(idir))
       prim(5+idir) = cpa_B0*eX(idir) + BY*eY(idir) + BZ*eZ(idir)
    end do
  end subroutine cpa_exact

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
         ixGmin3:ixGmax3,1:nw)

    double precision :: vA, prim(nerr), lfac
    integer :: i1, i2, i3, idir

    vA = cpa_vA()

    do i3 = ixGmin3, ixGmax3
       do i2 = ixGmin2, ixGmax2
          do i1 = ixGmin1, ixGmax1
             call cpa_exact(x(i1,i2,i3,1:3), 0.0d0, vA, prim)
             ! Lorentz factor (v^2 = v_A^2 eta^2 is uniform)
             lfac = 1.0d0/sqrt(1.0d0 - sum(prim(3:5)**2))
             w(i1,i2,i3,rho_) = prim(1)
             w(i1,i2,i3,p_)   = prim(2)
             do idir = 1, 3
                w(i1,i2,i3,mom(idir)) = lfac*prim(2+idir)   ! four-velocity
                w(i1,i2,i3,mag(idir)) = prim(5+idir)
             end do
          end do
       end do
    end do

    w(:,:,:,psi_) = 0

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

  !> L1 and Linf errors of the primitive variables with respect to the exact
  !> solution at the current time. Called at every log output (and at the end
  !> of the run) when typefilelog = 'special'. The host copy of the state is
  !> up to date when this is called. Appends one line per call to
  !> <base_filename>errors.csv
  subroutine cpa_print_errors()
    double precision, allocatable :: wp(:,:,:,:)
    double precision :: vA, prim(nerr), num(nerr), lfac, err
    double precision :: l1loc(nerr), linfloc(nerr), l1(nerr), linf(nerr)
    double precision :: ncells_loc, ncells
    integer :: iigrid, igrid, i1, i2, i3, iv, idir, ntot
    logical :: exists
    integer, parameter :: myunit = 91
    character(len=1024) :: filename

    vA = cpa_vA()
    l1loc = 0.0d0
    linfloc = 0.0d0
    ncells_loc = 0.0d0

    allocate(wp(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw))

    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       wp = ps(igrid)%w
       call phys_to_primitive(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
            ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,wp,ps(igrid)%x)

       do i3 = ixMlo3, ixMhi3
          do i2 = ixMlo2, ixMhi2
             do i1 = ixMlo1, ixMhi1
                call cpa_exact(ps(igrid)%x(i1,i2,i3,1:3), global_time, vA, prim)
                lfac = sqrt(1.0d0 + sum(wp(i1,i2,i3,mom(1:3))**2))
                num(1) = wp(i1,i2,i3,rho_)
                num(2) = wp(i1,i2,i3,p_)
                do idir = 1, 3
                   num(2+idir) = wp(i1,i2,i3,mom(idir))/lfac   ! 3-velocity
                   num(5+idir) = wp(i1,i2,i3,mag(idir))
                end do
                do iv = 1, nerr
                   err = abs(num(iv) - prim(iv))
                   l1loc(iv)   = l1loc(iv) + err
                   linfloc(iv) = max(linfloc(iv), err)
                end do
             end do
          end do
       end do
       ncells_loc = ncells_loc + dble((ixMhi1-ixMlo1+1)*(ixMhi2-ixMlo2+1)&
            *(ixMhi3-ixMlo3+1))
    end do

    deallocate(wp)

    call MPI_ALLREDUCE(l1loc, l1, nerr, MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierrmpi)
    call MPI_ALLREDUCE(linfloc, linf, nerr, MPI_DOUBLE_PRECISION, MPI_MAX, icomm, ierrmpi)
    call MPI_ALLREDUCE(ncells_loc, ncells, 1, MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierrmpi)
    l1 = l1/ncells   ! uniform grid: L1 = mean absolute error

    if (mype == 0) then
       ntot = nint(ncells**(1.0d0/3.0d0))
       filename = trim(base_filename)//'errors.csv'
       inquire(file=trim(filename), exist=exists)
       open(myunit, file=trim(filename), position='append', status='unknown')
       if (.not. exists) then
          write(myunit,'(a)',advance='no') 'N,it,time'
          do iv = 1, nerr
             write(myunit,'(a)',advance='no') ',L1_'//trim(errnames(iv))
          end do
          do iv = 1, nerr
             write(myunit,'(a)',advance='no') ',Linf_'//trim(errnames(iv))
          end do
          write(myunit,*)
       end if
       write(myunit,'(i0,",",i0,",",es24.16,16(",",es24.16))') ntot, it, &
            global_time, l1, linf
       close(myunit)
       write(*,'(a,es12.4,a,8es11.3)') ' CPA errors t=', global_time, &
            '  L1(rho,p,v,B):', l1
    end if

  end subroutine cpa_print_errors

end module mod_usr
