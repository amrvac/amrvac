!> multi-D SIR module with diffusion (physics routines)
!> SIR stands for Susceptible - Infected - Recovered
!> equations implemented are from:
!> International Journal of PDEs, Volume 2014, article ID 186437, Lofti et al.
!> http://dx.doi.org/10.1155/2014/186437
!>
module mod_sir_phys

  implicit none
  private

  integer, protected, public :: su_ = 1
  integer, protected, public :: in_ = 2
  integer, protected, public :: re_ = 3

  !> Whether particles module is added
  logical, public, protected              :: sir_particles = .false.

  !> Diffusion coefficients for s, i, r

  !> Diffusion coefficient for susceptible population (s)
  double precision, public, protected :: D1 = 0.01d0
  !> Diffusion coefficient for infected population (i)
  double precision, public, protected :: D2 = 0.01d0
  !> Diffusion coefficient for recovered population (r)
  double precision, public, protected :: D3 = 0.01d0

  !> population recruitment rate Lambda
  double precision, public, protected :: sir_Lambda  = 0.0d0
  !> death rate due to disease d
  double precision, public, protected :: sir_d = 0.0d0
  !> natural death rate mu
  double precision, public, protected :: sir_mu = 0.0d0
  !> recovery rate of the infected r
  double precision, public, protected :: sir_r = 0.1d0
  !> infection coefficient beta
  double precision, public, protected :: sir_beta = 0.1d0
  !> alpha parameters
  double precision, public, protected :: sir_alfa1 = 0.1d0
  double precision, public, protected :: sir_alfa2 = 0.1d0
  double precision, public, protected :: sir_alfa3 = 0.1d0

  !> Whether to handle the source term in split fashion
  logical :: sir_source_split = .false.

  ! Public methods
  public :: sir_phys_init

contains

  subroutine sir_params_read(files)
    use mod_global_parameters, only: unitpar, smalldouble
    character(len=*), intent(in) :: files(:)
    integer                      :: n
    character(len=20)            :: equation_name

    namelist /sir_list/ D1, D2, D3, sir_Lambda, sir_d, sir_mu, sir_r, sir_beta, &
                         sir_alfa1, sir_alfa2, sir_alfa3, sir_particles

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status='old')
       read(unitpar, sir_list, end=111)
111    close(unitpar)
    end do

    if(sir_Lambda<0.0d0) call mpistop("Lambda must be >=0")
    if(sir_d<0.0d0)      call mpistop("d must be >=0")
    if(sir_mu<0.0d0)     call mpistop("mu must be >=0")
    if(sir_r<0.0d0)      call mpistop("r must be >=0")
    if(sir_beta<0.0d0)   call mpistop("beta must be >=0")

    if(sir_alfa1<0.0d0) call mpistop("alfa1 must be >=0")
    if(sir_alfa2<0.0d0) call mpistop("alfa2 must be >=0")
    if(sir_alfa3<0.0d0) call mpistop("alfa3 must be >=0")

    if(D1<0.0d0) call mpistop("D1 must be >=0")
    if(D2<0.0d0) call mpistop("D2 must be >=0")
    if(D3<0.0d0) call mpistop("D3 must be >=0")

    if (dabs(D1+D2+D3)<smalldouble) then
      print *,"vanishing diffusion coefficients:",D1, D2, D3
      print *,"setting all diffusion coefficients to zero"
      D1=0.0d0
      D2=0.0d0
      D3=0.0d0
    endif
   
  end subroutine sir_params_read

  !> Write this module's parameters to a snapshot
  subroutine sir_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 0
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er
    integer                             :: idim

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)
  end subroutine sir_write_info

  subroutine sir_phys_init()
    use mod_global_parameters
    use mod_physics
    use mod_multigrid_coupling
    use mod_particles, only: particles_init

    call sir_params_read(par_files)

    physics_type = "sir"
    phys_energy  = .false.
    phys_req_diagonal = .false.
    use_particles = sir_particles

    ! Treat all variables as a density
    su_ = var_set_fluxvar("s", "s")
    in_ = var_set_fluxvar("i", "i")
    re_ = var_set_fluxvar("r", "r")

    ! set number of variables which need update ghostcells
    nwgc=nwflux

    ! Disable flux conservation near AMR boundaries, since we have no fluxes
    fix_conserve_global = .false.

    phys_get_cmax     => sir_get_cmax
    phys_get_cbounds  => sir_get_cbounds
    phys_get_flux     => sir_get_flux
    phys_to_conserved => sir_to_conserved
    phys_to_primitive => sir_to_primitive
    phys_add_source   => sir_add_source
    phys_get_dt       => sir_get_dt
    phys_write_info   => sir_write_info
    phys_check_params => sir_check_params
    phys_implicit_update   => sir_implicit_update
    phys_evaluate_implicit => sir_evaluate_implicit

    ! Initialize particles module
    if (sir_particles) then
       call particles_init()
       phys_req_diagonal = .true.
    end if

  end subroutine sir_phys_init

  subroutine sir_check_params
    use mod_multigrid_coupling
    use mod_global_parameters
    integer :: n

    if (any(flux_method /= fs_source)) then
       call mpistop("mod_sir requires flux_scheme = source")
    endif

    if (dabs(D1+D2+D3)<smalldouble .and. use_imex_scheme) then
       call mpistop("mod_sir for zero diffusion: explicit only : no IMEX possible")
    endif

    if (use_imex_scheme) then
       if((dabs(D1)<smalldouble.or.dabs(D2)<smalldouble).or.dabs(D3)<smalldouble) then
           call mpistop("imex requires all positive diffusion coefficients D1,D2,D3")
       endif
    endif

    if (use_imex_scheme) then
       use_multigrid=.true.
       ! Set boundary conditions for the multigrid solver
       do n = 1, 2*ndim
          select case (typeboundary(su_, n))
          case ('symm')
             ! d/dx s = 0
             mg%bc(n, mg_iphi)%bc_type = mg_bc_neumann
             mg%bc(n, mg_iphi)%bc_value = 0.0_dp
          case ('asymm')
             ! s = 0
             mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
             mg%bc(n, mg_iphi)%bc_value = 0.0_dp
          case ('cont')
             ! d/dx s = 0
             mg%bc(n, mg_iphi)%bc_type = mg_bc_neumann
             mg%bc(n, mg_iphi)%bc_value = 0.0_dp
          case ('periodic')
             ! Nothing to do here
          case default
             print *, "divb_multigrid warning: unknown b.c.: ", &
                  trim(typeboundary(su_, n))
             mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
             mg%bc(n, mg_iphi)%bc_value = 0.0_dp
          end select
       end do
    end if

  end subroutine sir_check_params

  subroutine sir_to_conserved(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)

    ! Do nothing (primitive and conservative are equal for sir module)
  end subroutine sir_to_conserved

  subroutine sir_to_primitive(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)

    ! Do nothing (primitive and conservative are equal for sir module)
  end subroutine sir_to_primitive

  subroutine sir_get_cmax(w, x, ixI^L, ixO^L, idim, cmax)
    use mod_global_parameters
    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    double precision, intent(inout)           :: cmax(ixI^S)

    cmax(ixO^S) = 0.0d0
  end subroutine sir_get_cmax

  subroutine sir_get_cbounds(wLC, wRC, wLp, wRp, x, ixI^L, ixO^L, idim, cmax, cmin)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wLC(ixI^S, nw), wRC(ixI^S,nw)
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)
    double precision, intent(inout) :: cmax(ixI^S)
    double precision, intent(inout), optional :: cmin(ixI^S)

    if (present(cmin)) then
       cmin(ixO^S) = 0.0d0
       cmax(ixO^S) = 0.0d0
    else
       cmax(ixO^S) = 0.0d0
    end if

  end subroutine sir_get_cbounds

  subroutine sir_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(inout) :: dtnew
    double precision                :: maxrate, DD,D12
    double precision                :: maxi_siterm,maxs_siterm,maxrate_si,maxrate_sir
    logical :: debug

    debug=.false.
    if (dabs(D1+D2+D3)<smalldouble) then
       dtnew = bigdouble
    else
       ! use dtdiffpar to set the courant parameter for diffusion term
       ! if explicit must be less than one
       ! dt < dx^2 / (2 * ndim * diffusion_coeff)
       D12=max(D1,D2)
       DD=max(D12,D3)
       dtnew = dtdiffpar * minval([ dx^D ])**2 / (2 * ndim * DD)
    endif
    if(mype==0.and.debug)then
        print *,'dtnew from diffusion=',dtnew
    endif

    ! Estimate time step for reactions
    ! in case of all coeffs zero: pure diffusion
    ! use courantpar to set coefficient
    if (.not.(sir_mu+sir_beta+sir_mu+sir_d+sir_r <= smalldouble)) then
       maxi_siterm=maxval(w(ixO^S,in_)/  &
           (1.0d0+sir_alfa1*w(ixO^S,su_)+sir_alfa2*w(ixO^S,in_)+sir_alfa3*w(ixO^S,su_)*w(ixO^S,in_)))
       maxs_siterm=maxval(w(ixO^S,su_)/ &
           (1.0d0+sir_alfa1*w(ixO^S,su_)+sir_alfa2*w(ixO^S,in_)+sir_alfa3*w(ixO^S,su_)*w(ixO^S,in_)))
       maxrate_si=max(sir_mu+sir_beta*maxi_siterm, &
                sir_beta*maxs_siterm-sir_mu-sir_d-sir_r)
       maxrate_sir=max(sir_mu,maxrate_si)
       dtnew = min(dtnew, courantpar/maxrate_sir)
       if(mype==0.and.debug)then
          print *,'dtnew from reactions=',courantpar/maxrate_sir
       endif
    endif

  end subroutine sir_get_dt

  ! There is no transport flux
  subroutine sir_get_flux(wC, w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wC(ixI^S, 1:nw)
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)
    double precision, intent(out)   :: f(ixI^S, nwflux)

    f(ixO^S, :) = 0.0d0
  end subroutine sir_get_flux

  ! w[iw]= w[iw]+qdt*S[wCT, qtC, x] where S is the source based on wCT within ixO
  subroutine sir_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,active)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision                :: lpl_s(ixO^S), lpl_i(ixO^S), lpl_r(ixO^S)
    double precision                :: siterm(ixO^S)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active

    if (qsourcesplit .eqv. sir_source_split) then
       if (.not.use_imex_scheme) then
          call sir_laplacian(ixI^L, ixO^L, wCT(ixI^S, su_), lpl_s)
          call sir_laplacian(ixI^L, ixO^L, wCT(ixI^S, in_), lpl_i)
          call sir_laplacian(ixI^L, ixO^L, wCT(ixI^S, re_), lpl_r)
       else
          ! all IMEX variants: only reactions
          lpl_s = 0.0d0
          lpl_i = 0.0d0
          lpl_r = 0.0d0
       end if

       siterm(ixO^S)= wCT(ixO^S, su_) * wCT(ixO^S, in_) / &
               (1.0d0+sir_alfa1*wCT(ixO^S, su_)+sir_alfa2*wCT(ixO^S, in_)+sir_alfa3*wCT(ixO^S, su_)*wCT(ixO^S, in_))

       w(ixO^S, su_) = w(ixO^S, su_) + qdt * (D1 * lpl_s  &
                  + sir_Lambda - sir_mu*wCT(ixO^S, su_) &
                  -sir_beta*siterm(ixO^S))
       w(ixO^S, in_) = w(ixO^S, in_) + qdt * (D2 * lpl_i  &
                  - (sir_mu+sir_d+sir_r)*wCT(ixO^S, in_) &
                  +sir_beta*siterm(ixO^S))
       w(ixO^S, re_) = w(ixO^S, re_) + qdt * (D3 * lpl_r  &
                  +sir_r*wCT(ixO^S, in_)-sir_mu*wCT(ixO^S, re_)) 

       active = .true.
    end if

  end subroutine sir_add_source

  !> Compute the Laplacian using a standard second order scheme. For now this
  !> method only works in slab geometries.
  subroutine sir_laplacian(ixI^L,ixO^L,var,lpl)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: var(ixI^S)
    double precision, intent(out) :: lpl(ixO^S)
    integer                       :: idir, jxO^L, hxO^L
    double precision              :: h_inv2

    if (slab) then
       lpl(ixO^S) = 0.0d0
       do idir = 1, ndim
          hxO^L=ixO^L-kr(idir,^D);
          jxO^L=ixO^L+kr(idir,^D);
          h_inv2 = 1/dxlevel(idir)**2
          lpl(ixO^S) = lpl(ixO^S) + h_inv2 * &
               (var(jxO^S) - 2 * var(ixO^S) + var(hxO^S))
       end do
    else
       call mpistop("sir_laplacian not implemented in this geometry")
    end if
  end subroutine sir_laplacian

  subroutine put_laplacians_onegrid(ixI^L,ixO^L,w)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)

    double precision                :: lpl_s(ixO^S), lpl_i(ixO^S), lpl_r(ixO^S)

    call sir_laplacian(ixI^L, ixO^L, w(ixI^S, su_), lpl_s)
    call sir_laplacian(ixI^L, ixO^L, w(ixI^S, in_), lpl_i)
    call sir_laplacian(ixI^L, ixO^L, w(ixI^S, re_), lpl_r)

    w(ixO^S,su_)=D1*lpl_s
    w(ixO^S,in_)=D2*lpl_i
    w(ixO^S,re_)=D3*lpl_r

  end subroutine put_laplacians_onegrid

  !> inplace update of psa==>F_im(psa)
  subroutine sir_evaluate_implicit(qtC,psa)
    use mod_global_parameters
    type(state), target :: psa(max_blocks)
    double precision, intent(in) :: qtC

    integer :: iigrid, igrid, level
    integer :: ixO^L

    ixO^L=ixG^LL^LSUB1;
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
       call put_laplacians_onegrid(ixG^LL,ixO^L,psa(igrid)%w)
    end do
    !$OMP END PARALLEL DO

  end subroutine sir_evaluate_implicit

  !> Implicit solve of psa=psb+dtfactor*dt*F_im(psa)
  subroutine sir_implicit_update(dtfactor,qdt,qtC,psa,psb)
    use mod_global_parameters
    use mod_forest
    use mod_multigrid_coupling

    type(state), target :: psa(max_blocks)
    type(state), target :: psb(max_blocks)
    double precision, intent(in) :: qdt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: dtfactor

    integer                      :: n
    double precision             :: res, max_residual, lambda

    integer                        :: iw_to,iw_from
    integer                        :: iigrid, igrid, id
    integer                        :: nc, lvl
    type(tree_node), pointer       :: pnode
    real(dp)                       :: fac

    ! Avoid setting a very restrictive limit to the residual when the time step
    ! is small (as the operator is ~ 1/(D * qdt))
    if (qdt < dtmin) then
        if(mype==0)then
            print *,'skipping implicit solve: dt too small!'
            print *,'Currently at time=',global_time,' time step=',qdt,' dtmin=',dtmin
        endif
        return
    endif
    max_residual = 1d-7/qdt

    mg%operator_type = mg_helmholtz
    call mg_set_methods(mg)

    if (.not. mg%is_allocated) call mpistop("multigrid tree not allocated yet")

    ! First handle the S variable ***************************************
    lambda           = 1/(dtfactor * qdt * D1)
    call helmholtz_set_lambda(lambda)

    !This is mg_copy_to_tree from psb state
    !!!  replaces::  call mg_copy_to_tree(su_, mg_irhs, factor=-lambda)
    iw_from=su_
    iw_to=mg_irhs
    fac=-lambda
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)
       ! Include one layer of ghost cells on grid leaves
       {^IFONED
       mg%boxes(id)%cc(0:nc+1, iw_to) = fac * &
            psb(igrid)%w(ixMlo1-1:ixMhi1+1, iw_from)
       }
       {^IFTWOD
       mg%boxes(id)%cc(0:nc+1, 0:nc+1, iw_to) = fac * &
            psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, iw_from)
       }
       {^IFTHREED
       mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iw_to) = fac * &
            psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, &
            ixMlo3-1:ixMhi3+1, iw_from)
       }
    end do

    !This is mg_copy_to_tree from psb state
    !!!  replaces::  call mg_copy_to_tree(su_, mg_iphi)
    iw_from=su_
    iw_to=mg_iphi
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)
       ! Include one layer of ghost cells on grid leaves
       {^IFONED
       mg%boxes(id)%cc(0:nc+1, iw_to) = fac * &
            psb(igrid)%w(ixMlo1-1:ixMhi1+1, iw_from)
       }
       {^IFTWOD
       mg%boxes(id)%cc(0:nc+1, 0:nc+1, iw_to) = fac * &
            psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, iw_from)
       }
       {^IFTHREED
       mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iw_to) = fac * &
            psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, &
            ixMlo3-1:ixMhi3+1, iw_from)
       }
    end do

    call mg_fas_fmg(mg, .true., max_res=res)
    do n = 1, 10
       call mg_fas_vcycle(mg, max_res=res)
       if (res < max_residual) exit
    end do

    !This is mg_copy_from_tree_gc for psa state
    !!! replaces:: call mg_copy_from_tree_gc(mg_iphi, su_)
    iw_from=mg_iphi
    iw_to=su_
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)
       ! Include one layer of ghost cells on grid leaves
       {^IFONED
       psa(igrid)%w(ixMlo1-1:ixMhi1+1, iw_to) = &
            mg%boxes(id)%cc(0:nc+1, iw_from)
       }
       {^IFTWOD
       psa(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, iw_to) = &
            mg%boxes(id)%cc(0:nc+1, 0:nc+1, iw_from)
       }
       {^IFTHREED
       psa(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, &
            ixMlo3-1:ixMhi3+1, iw_to) = &
            mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iw_from)
       }
    end do
    ! Done with the S variable ***************************************

    ! Next handle the I variable ***************************************
    lambda = 1/(dtfactor * qdt * D2)
    call helmholtz_set_lambda(lambda)

    !This is mg_copy_to_tree from psb state
    !!!  replaces::  call mg_copy_to_tree(in_, mg_irhs, factor=-lambda)
    iw_from=in_
    iw_to=mg_irhs
    fac=-lambda
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)
       ! Include one layer of ghost cells on grid leaves
       {^IFONED
       mg%boxes(id)%cc(0:nc+1, iw_to) = fac * &
            psb(igrid)%w(ixMlo1-1:ixMhi1+1, iw_from)
       }
       {^IFTWOD
       mg%boxes(id)%cc(0:nc+1, 0:nc+1, iw_to) = fac * &
            psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, iw_from)
       }
       {^IFTHREED
       mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iw_to) = fac * &
            psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, &
            ixMlo3-1:ixMhi3+1, iw_from)
       }
    end do

    !This is mg_copy_to_tree from psb state
    !!!  replaces::  call mg_copy_to_tree(in_, mg_iphi)
    iw_from=in_
    iw_to=mg_iphi
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)
       ! Include one layer of ghost cells on grid leaves
       {^IFONED
       mg%boxes(id)%cc(0:nc+1, iw_to) = fac * &
            psb(igrid)%w(ixMlo1-1:ixMhi1+1, iw_from)
       }
       {^IFTWOD
       mg%boxes(id)%cc(0:nc+1, 0:nc+1, iw_to) = fac * &
            psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, iw_from)
       }
       {^IFTHREED
       mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iw_to) = fac * &
            psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, &
            ixMlo3-1:ixMhi3+1, iw_from)
       }
    end do

    call mg_fas_fmg(mg, .true., max_res=res)
    do n = 1, 10
       call mg_fas_vcycle(mg, max_res=res)
       if (res < max_residual) exit
    end do

    !This is mg_copy_from_tree_gc for psa state
    !!! replaces:: call mg_copy_from_tree_gc(mg_iphi, in_)
    iw_from=mg_iphi
    iw_to=in_
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)
       ! Include one layer of ghost cells on grid leaves
       {^IFONED
       psa(igrid)%w(ixMlo1-1:ixMhi1+1, iw_to) = &
            mg%boxes(id)%cc(0:nc+1, iw_from)
       }
       {^IFTWOD
       psa(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, iw_to) = &
            mg%boxes(id)%cc(0:nc+1, 0:nc+1, iw_from)
       }
       {^IFTHREED
       psa(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, &
           ixMlo3-1:ixMhi3+1, iw_to) = &
            mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iw_from)
       }
    end do
    ! Done with the I variable ***************************************

    ! Next handle the R variable ***************************************
    lambda = 1/(dtfactor * qdt * D3)
    call helmholtz_set_lambda(lambda)

    !This is mg_copy_to_tree from psb state
    !!!  replaces::  call mg_copy_to_tree(re_, mg_irhs, factor=-lambda)
    iw_from=re_
    iw_to=mg_irhs
    fac=-lambda
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)
       ! Include one layer of ghost cells on grid leaves
       {^IFONED
       mg%boxes(id)%cc(0:nc+1, iw_to) = fac * &
            psb(igrid)%w(ixMlo1-1:ixMhi1+1, iw_from)
       }
       {^IFTWOD
       mg%boxes(id)%cc(0:nc+1, 0:nc+1, iw_to) = fac * &
            psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, iw_from)
       }
       {^IFTHREED
       mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iw_to) = fac * &
            psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, &
            ixMlo3-1:ixMhi3+1, iw_from)
       }
    end do

    !This is mg_copy_to_tree from psb state
    !!!  replaces::  call mg_copy_to_tree(re_, mg_iphi)
    iw_from=re_
    iw_to=mg_iphi
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)
       ! Include one layer of ghost cells on grid leaves
       {^IFONED
       mg%boxes(id)%cc(0:nc+1, iw_to) = fac * &
            psb(igrid)%w(ixMlo1-1:ixMhi1+1, iw_from)
       }
       {^IFTWOD
       mg%boxes(id)%cc(0:nc+1, 0:nc+1, iw_to) = fac * &
            psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, iw_from)
       }
       {^IFTHREED
       mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iw_to) = fac * &
            psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, &
            ixMlo3-1:ixMhi3+1, iw_from)
       }
    end do

    call mg_fas_fmg(mg, .true., max_res=res)
    do n = 1, 10
       call mg_fas_vcycle(mg, max_res=res)
       if (res < max_residual) exit
    end do

    !This is mg_copy_from_tree_gc for psa state
    !!! replaces:: call mg_copy_from_tree_gc(mg_iphi, re_)
    iw_from=mg_iphi
    iw_to=re_
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)
       ! Include one layer of ghost cells on grid leaves
       {^IFONED
       psa(igrid)%w(ixMlo1-1:ixMhi1+1, iw_to) = &
            mg%boxes(id)%cc(0:nc+1, iw_from)
       }
       {^IFTWOD
       psa(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, iw_to) = &
            mg%boxes(id)%cc(0:nc+1, 0:nc+1, iw_from)
       }
       {^IFTHREED
       psa(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, &
           ixMlo3-1:ixMhi3+1, iw_to) = &
            mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iw_from)
       }
    end do
    ! Done with the R variable ***************************************

  end subroutine sir_implicit_update

end module mod_sir_phys
