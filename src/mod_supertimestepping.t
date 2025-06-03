!> Generic supertimestepping method
!> 1) in amrvac.par in sts_list set the following parameters which have the default values: 
!> sts_dtpar=0.9,sts_ncycles=1000,sts_method=1,sourcetype_sts=2
!> These parametes are general for all the methods added TODO: check if there is any need
!> to have terms implemented with different sets of parameters, and these cannot be general anymore
!> 2) then add programatically in the code a term with the subroutine
!> add_sts_method
!> This method takes as parameters a function which calculated the explicit timestep
!> associated with the term, a subroutine which sets the source term 
!> types for the BC and the BC are generated from  the variables startVar:endVar
!> flux conservation (fixconserve) is done for the variables specified by  ixChangeStart, ixChangeN, ixChangeFixC
!> The following two steps are done in this way as in fortran it is not allowed to pass null function pointers as parameters:
!> 3)in order to  to have hooks before_first_cycle, after_last_cycle (e.g. conversion from e_tot to e_int before first sts cycle
!> and back from  e_int to e_tot after the last STS cycle  for the thermal conductivity module) add them just afterwards with the subroutine
!> set_conversion_methods_to_head
!> 4) to add the hook for error handling (e.g check small values in the thermal conductivity module ) 
!> call set_error_handling_to_head which takes as parameter a subroutine
!> the error handling subroutine is called before setting BC
module mod_supertimestepping
  use mod_geometry
  use mod_comm_lib, only: mpistop
  implicit none
  private

  public :: is_sts_initialized
  public :: sts_init
  public :: add_sts_method,set_conversion_methods_to_head,set_error_handling_to_head
  public :: sts_add_source
  public :: set_dt_sts_ncycles
  public :: sourcetype_sts, sourcetype_sts_prior, sourcetype_sts_after, sourcetype_sts_split

  ! input parameters from parameter file
  !> the coefficient that multiplies the sts dt
  double precision :: sts_dtpar=0.8d0

  !The following is used only for method 2, not input parameter TODO check if we want as input parameter
  double precision,parameter :: nu_sts = 0.5d0
  !> the maximum number of subcycles
  integer :: sts_ncycles=1000
  integer :: sts_method = 1
  integer, parameter :: sourcetype_sts_prior =0
  integer, parameter :: sourcetype_sts_after =1
  integer, parameter :: sourcetype_sts_split =2
  integer :: sourcetype_sts = sourcetype_sts_split
  !> Whether to conserve fluxes at the current partial step
  logical :: fix_conserve_at_step = .true.
  logical :: sts_initialized = .false.

  abstract interface

    !>interface for setting sources in the derived type
    subroutine subr1(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux)
      use mod_global_parameters
      integer, intent(in) :: ixI^L, ixO^L, igrid, nflux
      double precision, intent(in) ::  x(ixI^S,1:ndim)
      double precision, intent(inout) :: wres(ixI^S,1:nw), w(ixI^S,1:nw)
      double precision, intent(in) :: my_dt
      logical, intent(in) :: fix_conserve_at_step
    end subroutine subr1

    !>interface for the function which gets the timestep in dtnew in the derived type
    function subr2(w,ixG^L,ix^L,dx^D,x) result(dtnew)
      use mod_global_parameters
      integer, intent(in) :: ixG^L, ix^L
      double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
      double precision, intent(in) :: w(ixG^S,1:nw) 
      double precision :: dtnew
    end function subr2

    !>interface for error handling subroutine in the derived type
    subroutine subr_e(w, x, ixI^L, ixO^L, step)
      use mod_global_parameters
      use mod_small_values
      integer, intent(in)             :: ixI^L,ixO^L
      double precision, intent(inout) :: w(ixI^S,1:nw)
      double precision, intent(in)    :: x(ixI^S,1:ndim)
      integer, intent(in)    :: step
    end subroutine subr_e

    !>interface for the subroutines before_first_cycle and after_last_cycle in the derived type
    subroutine subr5(ixI^L, ixO^L, w, x)
      use mod_global_parameters
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(inout) :: w(ixI^S,1:nw) 
    end subroutine subr5

    !>for the subroutines in this module, which do not depend on the term, but
    !>on the parameter sts_method = 1/2 in the parameter file
    !>sts_add_source
    subroutine subr3(dt)
      double precision,intent(in) :: dt
    end subroutine subr3

    !>sts_get_ncycles 
    function subr4(dt,dtnew,dt_modified) result(s)
      double precision,intent(in) :: dtnew
      double precision,intent(inout) :: dt
      logical,intent(inout) :: dt_modified
      integer :: s
    end function subr4

  end interface

  type sts_term

    double precision :: dt_expl
    integer, public :: s

    !>types used for send/recv ghosts, see mod_ghostcells_update
    integer, dimension(-1:1^D&) :: type_send_srl_sts_1, type_recv_srl_sts_1
    integer, dimension(-1:1^D&) :: type_send_r_sts_1
    integer, dimension( 0:3^D&) :: type_recv_r_sts_1
    integer, dimension( 0:3^D&) :: type_recv_p_sts_1, type_send_p_sts_1

    integer, dimension(-1:1^D&) :: type_send_srl_sts_2, type_recv_srl_sts_2
    integer, dimension(-1:1^D&) :: type_send_r_sts_2
    integer, dimension( 0:3^D&) :: type_recv_r_sts_2
    integer, dimension( 0:3^D&) :: type_recv_p_sts_2, type_send_p_sts_2

    integer :: startVar
    integer :: endVar
    !> number of flux involved in STS update
    integer :: nflux
    integer :: startwbc
    integer :: nwbc
    logical :: types_initialized
    logical :: evolve_magnetic_field
    procedure (subr1), pointer, nopass :: sts_set_sources
    procedure (subr2), pointer, nopass :: sts_getdt
    procedure (subr5), pointer, nopass :: sts_before_first_cycle, sts_after_last_cycle
    procedure (subr_e), pointer, nopass :: sts_handle_errors
    type(sts_term), pointer :: next

  end type sts_term

  type(sts_term), pointer :: head_sts_terms
  !The following two subroutine/function  pointers 
  !make the difference between the two STS methods implemented
  procedure (subr3), pointer :: sts_add_source
  procedure (subr4), pointer :: sts_get_ncycles

contains

  !> Initialize sts module
  subroutine sts_init()
    use mod_global_parameters
    use mod_physics
    if(.not. sts_initialized) then  
      nullify(head_sts_terms)
      call sts_params_read(par_files)
      sts_dtpar=sts_dtpar/dble(ndim)
      sts_initialized = .true.
      if(sts_method .eq. 1) then
        sts_add_source => sts_add_source1
        sts_get_ncycles => sts_get_ncycles1
      else if(sts_method .eq. 2) then
        sts_add_source => sts_add_source2
        sts_get_ncycles => sts_get_ncycles2
      else
        call mpistop("Unknown sts method")
      end if
    endif

  end subroutine sts_init

  pure function is_sts_initialized() result(res)
    logical :: res
    if (sts_initialized) then
      res = associated(head_sts_terms)
    else
      res = .false.
    endif
  end function is_sts_initialized

  !> Read module parameters from par file
  subroutine sts_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /sts_list/ sts_dtpar,sts_ncycles,sts_method,sourcetype_sts

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, sts_list, end=111)
111    close(unitpar)
    end do

  end subroutine sts_params_read

  !> subroutine which added programatically a term to be calculated using STS
  !> Params: 
  !> sts_getdt function calculates the explicit timestep for this term
  !> sts_set_sources subroutine sets the source term
  !> startVar, endVar, nflux indices of start, end, and number of the variables that need fix conservation
  !> startwbc, nwbc indices of start and number of the variables that need ghost cell exchange
  !> These terms implemented by an element of the derived type sts_term are put in a linked list
  subroutine add_sts_method(sts_getdt, sts_set_sources, startVar, nflux, startwbc, nwbc, evolve_B)
    use mod_global_parameters
    use mod_ghostcells_update

    integer, intent(in) :: startVar, nflux, startwbc, nwbc
    logical, intent(in) :: evolve_B

    interface

      subroutine sts_set_sources(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux)
        use mod_global_parameters
        use mod_fix_conserve
        integer, intent(in) :: ixI^L, ixO^L, igrid, nflux
        double precision, intent(in) ::  x(ixI^S,1:ndim)
        double precision, intent(inout) :: wres(ixI^S,1:nw), w(ixI^S,1:nw)
        double precision, intent(in) :: my_dt
        logical, intent(in) :: fix_conserve_at_step
      end subroutine sts_set_sources
  
      function sts_getdt(w,ixG^L,ix^L,dx^D,x) result(dtnew)
        use mod_global_parameters
        integer, intent(in) :: ixG^L, ix^L
        double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
        double precision, intent(in) :: w(ixG^S,1:nw) 
        double precision :: dtnew
      end function sts_getdt

    end interface

    type(sts_term), pointer  :: temp
    allocate(temp)

    temp%sts_getdt => sts_getdt
    temp%sts_set_sources => sts_set_sources
    temp%sts_before_first_cycle => null()
    temp%sts_after_last_cycle => null()
    temp%sts_handle_errors => null()
    temp%startVar = startVar
    temp%endVar= startVar+nflux-1
    temp%nflux = nflux
    temp%startwbc = startwbc
    temp%nwbc = nwbc
    temp%types_initialized = .false.
    temp%evolve_magnetic_field=evolve_B

    temp%next => head_sts_terms
    head_sts_terms => temp

  end subroutine add_sts_method

  !> Set the hooks called before the first cycle and after the last cycle in the STS update
  !> This method should be called after add_sts_method. The hooks are added to the last term added with this subroutine
  !> Params: sts_before_first_cycle, sts_after_last_cycle subroutines which implement the hooks called before first cycle and after last cycle
  subroutine set_conversion_methods_to_head(sts_before_first_cycle, sts_after_last_cycle)
    interface
      subroutine sts_before_first_cycle(ixI^L, ixO^L, w, x)
        use mod_global_parameters
        integer, intent(in) :: ixI^L, ixO^L
        double precision, intent(in) :: x(ixI^S,1:ndim)
        double precision, intent(inout) :: w(ixI^S,1:nw) 
      end subroutine sts_before_first_cycle

      subroutine sts_after_last_cycle(ixI^L, ixO^L, w, x)
        use mod_global_parameters
        integer, intent(in) :: ixI^L, ixO^L
        double precision, intent(in) :: x(ixI^S,1:ndim)
        double precision, intent(inout) :: w(ixI^S,1:nw)
      end subroutine sts_after_last_cycle
    end interface

    head_sts_terms%sts_before_first_cycle => sts_before_first_cycle
    head_sts_terms%sts_after_last_cycle => sts_after_last_cycle

  end subroutine set_conversion_methods_to_head

  !> Set the hook of error handling in the STS update. This method is called before updating the BC.
  !> This method should be called after add_sts_method. The hook is added to the last term added with this subroutine.
  !> Param: sts_error_handing the subroutine which handles the errors
  subroutine set_error_handling_to_head(sts_error_handling)
    interface
      subroutine sts_error_handling(w, x, ixI^L, ixO^L, step)
        use mod_global_parameters
        use mod_small_values
        integer, intent(in)             :: ixI^L,ixO^L
        double precision, intent(inout) :: w(ixI^S,1:nw)
        double precision, intent(in)    :: x(ixI^S,1:ndim)
        integer, intent(in)    :: step
      end subroutine sts_error_handling
    end interface
    head_sts_terms%sts_handle_errors => sts_error_handling

  end subroutine set_error_handling_to_head

  !> method used to set the number of cycles for the STS1 method
  function sts_get_ncycles1(dt,dtnew,dt_modified) result(is)
    double precision,intent(in) :: dtnew
    double precision,intent(inout) :: dt
    logical,intent(inout) :: dt_modified
    integer :: is

    double precision    :: ss

    !!ss is now limit of dt because of sts_ncycles
    ss = dtnew*((2.d0*sts_ncycles+1)**2-9.d0)/16.d0
    if(dt>ss) then
      dt_modified = .true.
      dt = ss
      is = sts_ncycles
    else
      ss = dt/dtnew
      ! get number of sub-steps of supertime stepping (Meyer 2012 MNRAS 422,2102)
      if(ss .le. 1.d0) then
        is=1
      else
        is=ceiling((dsqrt(9.d0+16.d0*ss)-1.d0)*0.5d0)
        is=is/2*2+1
      end if
    end if

  end function sts_get_ncycles1

  !> method used to set the number of cycles for the STS2 method
  function sts_get_ncycles2(dt,dtnew,dt_modified) result(is)
    double precision,intent(in) :: dtnew
    double precision,intent(inout) :: dt
    logical,intent(inout) :: dt_modified
    integer :: is

    double precision    :: ss,rr
    integer:: ncycles

    rr = dt/dtnew
    !print*, dt, " --DTEXPL-- ", dtnew, ", rr ",rr
    ncycles = sts_ncycles 
    !print*, "NCYCLES BEFORE ",ncycles
    ss=sum_chev(nu_sts,ncycles,rr)
    !print*, "NCYCLES AFTER ",ncycles
    is = ncycles
    !print*, "SUMCHEV ", ss, " NCYCLES ", is
    if(ss < rr) then
      dt_modified = .true.
      dt = ss *  dtnew
    endif

  end function sts_get_ncycles2

  !> This sets the explicit dt and calculates the number of cycles for each of the terms implemented with STS.
  function set_dt_sts_ncycles(my_dt) result(dt_modified)
    use mod_global_parameters

    double precision,intent(inout) :: my_dt
    double precision :: my_dt1
    logical :: dt_modified, dt_modified1, dt_modified2

    double precision :: dtnew,dtmin_mype
    double precision    :: dx^D, ss
    integer:: iigrid, igrid, ncycles
    type(sts_term), pointer  :: temp,oldTemp
    nullify(oldTemp)
    temp => head_sts_terms
    dt_modified = .false.
    do while(associated(temp))
      dt_modified2 = .false.
      dtmin_mype=bigdouble
      !$OMP PARALLEL DO PRIVATE(igrid,dx^D) REDUCTION(min:dtmin_mype)
      do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
         ! maybe the following global variables are needed in get_dt!
         ! next few lines ensure correct usage of routines like divvector etc
         ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
         block=>ps(igrid)
         ! end maybe the following global variables are needed in get_dt!!!!!!!
         dx^D=rnode(rpdx^D_,igrid);
         dtmin_mype=min(dtmin_mype, sts_dtpar * temp%sts_getdt(ps(igrid)%w,ixG^LL,ixM^LL,dx^D,ps(igrid)%x))
      end do
      !$OMP END PARALLEL DO
      call MPI_ALLREDUCE(dtmin_mype,dtnew,1,MPI_DOUBLE_PRECISION,MPI_MIN,icomm,ierrmpi)
      temp%s = sts_get_ncycles(my_dt,dtnew,dt_modified2)

      !print*, "NCYCLES ", temp%s, dt_modified2, my_dt, dtnew
      temp%dt_expl = dtnew
 
      ! Note that as for some term it may happen that the dt is modified: it may be reduced if the
      ! number of cycles is overpassed, the list has to be reiterated to update ncycles for previous 
      ! terms which did not modify dt TODO add pointer to previous and loop backward to update
      if(dt_modified2) then
        dt_modified = .true.
        !reiterate all the other sts elements and recalculate s 
        oldTemp => head_sts_terms
        my_dt1 = my_dt
        dt_modified1 = .false.
        do while(.not. associated(oldTemp,temp))  
         oldTemp%s = sts_get_ncycles(my_dt1,oldTemp%dt_expl,dt_modified1) 
         !check dt is not modified again, and this should not happen, except for bug in sts_get_ncycles1,2
         if(dt_modified1) call mpistop("sts dt modified twice")
         oldTemp=>oldTemp%next
        end do
      end if
      temp=>temp%next

    end do

  end function set_dt_sts_ncycles

  pure FUNCTION chev(j,nu,N)
    use mod_constants

    double precision, INTENT(IN) :: nu
    INTEGER, INTENT(IN)       :: j, N
    double precision             :: chev

    chev = 1d0 / ((-1d0 + nu)*cos(((2d0*j - 1d0) / N)* (dpi/2d0)) + 1d0 + nu)

  END FUNCTION chev

  FUNCTION sum_chev(nu,N,limMax)
    double precision, intent(in) :: nu,limmax
    integer, intent(inout)       ::  N
    double precision             :: sum_chev, tmp

    integer :: j

    j=1
    sum_chev = 0d0
    do while (j < N .and. sum_chev < limMax)
      sum_chev = sum_chev + chev(j,nu,N)
      j=j+1
    enddo
   N=j-1 
  END FUNCTION sum_chev

  !TODO the following not used
!  PURE FUNCTION total_chev(nu,N)
!    double precision, INTENT(IN) :: nu
!    INTEGER, INTENT(IN)       :: N
!    double precision             :: total_chev
!
!    total_chev = N/(2d0*dsqrt(nu)) * ( (1d0 + dsqrt(nu))**(2d0*N) - (1d0 - dsqrt(nu))**(2d0*N) ) / &
!         ( (1d0 + dsqrt(nu))**(2d0*N) + (1d0 - dsqrt(nu))**(2d0*N) )
!
!  END FUNCTION total_chev

  !> Iterates all the terms implemented with STS and adds the sources
  !> STS method 2 implementation
  subroutine sts_add_source2(my_dt)
  ! Turlough Downes 2006,2007
    use mod_ghostcells_update
    use mod_global_parameters
    use mod_fix_conserve
    use mod_physics

    double precision, intent(in) :: my_dt
    double precision, allocatable :: bj(:)
    double precision :: sumbj,dtj  

    integer:: iigrid, igrid, j, ixC^L
    logical :: stagger_flag=.false., prolong_flag=.false., coarsen_flag=.false.
    type(sts_term), pointer  :: temp

    ! do not fill physical boundary conditions
    bcphys=.false.

    fix_conserve_at_step = time_advance .and. levmax>levmin

    temp => head_sts_terms
    do while(associated(temp))

      if(.not.temp%evolve_magnetic_field) then
        ! not do fix conserve and getbc for staggered values
        stagger_flag=stagger_grid
        stagger_grid=.false.
      else if(stagger_grid) then
        ixCmax^D=ixMhi^D;
        ixCmin^D=ixMlo^D-1;
      end if

      call init_comm_fix_conserve(1,ndim,temp%nflux)

      if(associated(temp%sts_before_first_cycle)) then
        prolong_flag=prolongprimitive
        coarsen_flag=coarsenprimitive
        prolongprimitive=.false.
        coarsenprimitive=.false.
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
          block=>ps(igrid)
          call temp%sts_before_first_cycle(ixG^LL,ixG^LL,ps(igrid)%w,ps(igrid)%x)  
        end do
      end if

      allocate(bj(1:temp%s))
      do j=1,temp%s
        bj(j) = chev(j,nu_sts,sts_ncycles)
      end do

      type_send_srl=>temp%type_send_srl_sts_1
      type_recv_srl=>temp%type_recv_srl_sts_1
      type_send_r=>temp%type_send_r_sts_1
      type_recv_r=>temp%type_recv_r_sts_1
      type_send_p=>temp%type_send_p_sts_1
      type_recv_p=>temp%type_recv_p_sts_1

      if(.not. temp%types_initialized) then
        call create_bc_mpi_datatype(temp%startwbc,temp%nwbc)
        if(temp%nflux>temp%nwbc) then
          ! prepare types for the changed no-need-ghost-update variables in the last getbc
          type_send_srl=>temp%type_send_srl_sts_2
          type_recv_srl=>temp%type_recv_srl_sts_2
          type_send_r=>temp%type_send_r_sts_2
          type_recv_r=>temp%type_recv_r_sts_2
          type_send_p=>temp%type_send_p_sts_2
          type_recv_p=>temp%type_recv_p_sts_2
          call create_bc_mpi_datatype(temp%startVar,temp%nflux)
          type_send_srl=>temp%type_send_srl_sts_1
          type_recv_srl=>temp%type_recv_srl_sts_1
          type_send_r=>temp%type_send_r_sts_1
          type_recv_r=>temp%type_recv_r_sts_1
          type_send_p=>temp%type_send_p_sts_1
          type_recv_p=>temp%type_recv_p_sts_1
        end if
        temp%types_initialized = .true.
      end if

      sumbj=0.d0
      do j=1,temp%s
        if(j .eq. temp%s .and. (sumbj + bj(j)) * temp%dt_expl > my_dt) then
          dtj = my_dt - sumbj * temp%dt_expl
        else
          dtj = bj(j)* temp%dt_expl
        end if
        sumbj = sumbj + bj(j)
        if(stagger_grid) then
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
            ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
            block=>ps(igrid)
            call temp%sts_set_sources(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x,ps1(igrid)%w,fix_conserve_at_step,dtj,igrid,temp%nflux)
            if(temp%nflux>ndir) then
              ps(igrid)%w(ixM^T,temp%startVar)=ps(igrid)%w(ixM^T,temp%startVar)+dtj*ps1(igrid)%w(ixM^T,temp%startVar)
            end if
            ps(igrid)%ws(ixC^S,1:nws)=ps(igrid)%ws(ixC^S,1:nws)+dtj*ps1(igrid)%w(ixC^S,iw_mag(1:nws))
            call phys_face_to_center(ixM^LL,ps(igrid))
          end do
          !$OMP END PARALLEL DO
        else
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
            ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
            block=>ps(igrid)
            call temp%sts_set_sources(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x,ps1(igrid)%w,fix_conserve_at_step,dtj,igrid,temp%nflux)
            ps(igrid)%w(ixM^T,temp%startVar:temp%endVar)=ps(igrid)%w(ixM^T,temp%startVar:temp%endVar)+&
              dtj*ps1(igrid)%w(ixM^T,temp%startVar:temp%endVar)
          end do
          !$OMP END PARALLEL DO
        end if
        !fix conserve the fluxes set in the STS method
        if(fix_conserve_at_step) then
          call recvflux(1,ndim)
          call sendflux(1,ndim)
          call fix_conserve(ps,1,ndim,temp%startVar,temp%nflux)
          if(stagger_grid) then
            call fix_edges(ps,1,ndim)
            ! fill the cell-center values from the updated staggered variables
            !$OMP PARALLEL DO PRIVATE(igrid)
            do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
              call phys_face_to_center(ixG^LL,ps(igrid))
            end do
            !$OMP END PARALLEL DO
          end if
        end if
        if(associated(temp%sts_handle_errors)) then
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
            ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
            block=>ps(igrid)
            call temp%sts_handle_errors(ps(igrid)%w,ps(igrid)%x,ixG^LL,ixM^LL,j)
          end do
          !$OMP END PARALLEL DO
        end if

        if(temp%nflux>temp%nwbc.and.temp%s==j) then
          ! include the changed no-need-ghost-update variables in the last getbc
          type_send_srl=>temp%type_send_srl_sts_2
          type_recv_srl=>temp%type_recv_srl_sts_2
          type_send_r=>temp%type_send_r_sts_2
          type_recv_r=>temp%type_recv_r_sts_2
          type_send_p=>temp%type_send_p_sts_2
          type_recv_p=>temp%type_recv_p_sts_2
          call getbc(global_time,0.d0,ps,temp%startVar,temp%nflux)
        else
          call getbc(global_time,0.d0,ps,temp%startwbc,temp%nwbc)
        end if
      end do

      if(associated(temp%sts_after_last_cycle)) then 
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
          block=>ps(igrid)
          call temp%sts_after_last_cycle(ixG^LL,ixG^LL,ps(igrid)%w,ps(igrid)%x)
        end do
        prolongprimitive  = prolong_flag
        coarsenprimitive = coarsen_flag
      end if
      deallocate(bj)

      if(.not.temp%evolve_magnetic_field) then
        ! restore stagger_grid value
        stagger_grid=stagger_flag
      end if

      temp=>temp%next
    end do

    if(associated(head_sts_terms)) then
      ! point bc mpi data type back to full type for (M)HD
      type_send_srl=>type_send_srl_f
      type_recv_srl=>type_recv_srl_f
      type_send_r=>type_send_r_f
      type_recv_r=>type_recv_r_f
      type_send_p=>type_send_p_f
      type_recv_p=>type_recv_p_f
    end if

    bcphys=.true.

    if(phys_partial_ionization) then
      ! update temperature variable in w
      !$OMP PARALLEL DO PRIVATE(igrid)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        call phys_update_temperature(ixG^LL,ixG^LL,ps(igrid)%w,ps(igrid)%w,ps(igrid)%x)
      end do
      !$OMP END PARALLEL DO
    end if

  end subroutine sts_add_source2

  !> Iterates all the terms implemented with STS and adds the sources
  !> STS method 1 implementation
  subroutine sts_add_source1(my_dt)
  ! Meyer 2012 MNRAS 422,2102
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_fix_conserve
    use mod_physics
    use mod_amr_solution_node, only: alloc_state

    double precision, intent(in) :: my_dt
    double precision :: dtj
    double precision :: omega1,cmu,cmut,cnu,cnut,one_mu_nu
    double precision, allocatable :: bj(:)
    integer:: iigrid, igrid, j, ixC^L, ixGext^L
    logical :: evenstep, stagger_flag=.false., prolong_flag=.false., coarsen_flag=.false., total_energy_flag=.true.
    type(sts_term), pointer  :: temp
    type(state), dimension(:), pointer :: tmpPs1, tmpPs2

    ! do not fill physical boundary conditions
    bcphys=.false.

    fix_conserve_at_step = time_advance .and. levmax>levmin

    temp => head_sts_terms
    do while(associated(temp))

      if(.not.temp%evolve_magnetic_field) then
        ! not do fix conserve and getbc for staggered values
        stagger_flag=stagger_grid
        stagger_grid=.false.
      else if(stagger_grid) then
        ixCmax^D=ixMhi^D;
        ixCmin^D=ixMlo^D-1;
      end if

      call init_comm_fix_conserve(1,ndim,temp%nflux)

      if(associated(temp%sts_before_first_cycle)) then
        prolong_flag = prolongprimitive
        coarsen_flag = coarsenprimitive
        prolongprimitive = .false.
        coarsenprimitive = .false.
        total_energy_flag=phys_total_energy
        phys_total_energy=.false.
        !$OMP PARALLEL DO PRIVATE(igrid)
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
          block=>ps(igrid)
          call temp%sts_before_first_cycle(ixG^LL,ixG^LL,ps(igrid)%w,ps(igrid)%x)
          if(.not. allocated(ps2(igrid)%w)) allocate(ps2(igrid)%w(ixG^T,1:nw))
          if(.not. allocated(ps3(igrid)%w)) allocate(ps3(igrid)%w(ixG^T,1:nw))
          if(.not. allocated(ps4(igrid)%w)) allocate(ps4(igrid)%w(ixG^T,1:nw))
          ps1(igrid)%w(ixG^T,1:nw)=ps(igrid)%w(ixG^T,1:nw)
          ps2(igrid)%w(ixG^T,1:nw)=ps(igrid)%w(ixG^T,1:nw)
        end do
        !$OMP END PARALLEL DO
      else
        if(stagger_grid) then
          ixGext^L=ixG^LL^LADD1;
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail; igrid=igrids(iigrid);
            if(.not. allocated(ps2(igrid)%w)) then
              call alloc_state(igrid, ps2(igrid), ixG^LL, ixGext^L, .false.)
            end if
            if(.not. allocated(ps3(igrid)%w)) allocate(ps3(igrid)%w(ixG^T,1:nw))
            if(.not. allocated(ps4(igrid)%w)) allocate(ps4(igrid)%w(ixG^T,1:nw))
            ps1(igrid)%w=ps(igrid)%w
            ps2(igrid)%w=ps(igrid)%w
            ps1(igrid)%ws=ps(igrid)%ws
            ps2(igrid)%ws=ps(igrid)%ws
          end do
          !$OMP END PARALLEL DO
        else
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail; igrid=igrids(iigrid);
            if(.not. allocated(ps2(igrid)%w)) allocate(ps2(igrid)%w(ixG^T,1:nw))
            if(.not. allocated(ps3(igrid)%w)) allocate(ps3(igrid)%w(ixG^T,1:nw))
            if(.not. allocated(ps4(igrid)%w)) allocate(ps4(igrid)%w(ixG^T,1:nw))
            ps1(igrid)%w(ixG^T,1:nw)=ps(igrid)%w(ixG^T,1:nw)
            ps2(igrid)%w(ixG^T,1:nw)=ps(igrid)%w(ixG^T,1:nw)
          end do
          !$OMP END PARALLEL DO
        end if
      end if

      allocate(bj(0:temp%s))
      bj(0)=1.d0/3.d0
      bj(1)=bj(0)
      if(temp%s>1) then
        omega1=4.d0/dble(temp%s**2+temp%s-2)
        cmut=omega1/3.d0
      else
        omega1=0.d0
        cmut=1.d0
      end if

      type_send_srl=>temp%type_send_srl_sts_1
      type_recv_srl=>temp%type_recv_srl_sts_1
      type_send_r=>temp%type_send_r_sts_1
      type_recv_r=>temp%type_recv_r_sts_1
      type_send_p=>temp%type_send_p_sts_1
      type_recv_p=>temp%type_recv_p_sts_1

      if(.not. temp%types_initialized) then 
        call create_bc_mpi_datatype(temp%startwbc,temp%nwbc)
        if(temp%nflux>temp%nwbc) then
          ! prepare types for the changed no-need-ghost-update variables in the last getbc
          type_send_srl=>temp%type_send_srl_sts_2
          type_recv_srl=>temp%type_recv_srl_sts_2
          type_send_r=>temp%type_send_r_sts_2
          type_recv_r=>temp%type_recv_r_sts_2
          type_send_p=>temp%type_send_p_sts_2
          type_recv_p=>temp%type_recv_p_sts_2
          call create_bc_mpi_datatype(temp%startVar,temp%nflux)
          type_send_srl=>temp%type_send_srl_sts_1
          type_recv_srl=>temp%type_recv_srl_sts_1
          type_send_r=>temp%type_send_r_sts_1
          type_recv_r=>temp%type_recv_r_sts_1
          type_send_p=>temp%type_send_p_sts_1
          type_recv_p=>temp%type_recv_p_sts_1
        end if
        temp%types_initialized = .true.
      end if
      dtj = cmut*my_dt
      if(stagger_grid) then
        !$OMP PARALLEL DO PRIVATE(igrid)
        do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
          block=>ps(igrid)
          ps4(igrid)%w=zero
          call temp%sts_set_sources(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x,ps4(igrid)%w,fix_conserve_at_step,dtj,igrid,temp%nflux)
          !!!eq solved: dU/dt = S, ps3 is stored S^n
            ps3(igrid)%w(ixC^S,temp%startVar:temp%endVar) = my_dt * ps4(igrid)%w(ixC^S,temp%startVar:temp%endVar)
          if(temp%nflux>ndir) then
            ps1(igrid)%w(ixM^T,temp%startVar) = ps1(igrid)%w(ixM^T,temp%startVar) + cmut * ps3(igrid)%w(ixM^T,temp%startVar)
          end if
          ps1(igrid)%ws(ixC^S,1:nws) = ps1(igrid)%ws(ixC^S,1:nws) + cmut * ps3(igrid)%w(ixC^S,iw_mag(1:nws))
          call phys_face_to_center(ixM^LL,ps1(igrid))
        end do
        !$OMP END PARALLEL DO
      else
        !$OMP PARALLEL DO PRIVATE(igrid)
        do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
          block=>ps(igrid)
          call temp%sts_set_sources(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x,ps4(igrid)%w,fix_conserve_at_step,dtj,igrid,temp%nflux)
          !!!eq solved: dU/dt = S, ps3 is stored S^n
          ps3(igrid)%w(ixM^T,temp%startVar:temp%endVar) = my_dt * ps4(igrid)%w(ixM^T,temp%startVar:temp%endVar)
          ps1(igrid)%w(ixM^T,temp%startVar:temp%endVar) = ps1(igrid)%w(ixM^T,temp%startVar:temp%endVar) + &
            cmut * ps3(igrid)%w(ixM^T,temp%startVar:temp%endVar)
        end do
        !$OMP END PARALLEL DO
      end if
      if(fix_conserve_at_step) then
        call recvflux(1,ndim)
        call sendflux(1,ndim)
        call fix_conserve(ps1,1,ndim,temp%startVar,temp%nflux)
        if(stagger_grid) then
          call fix_edges(ps1,1,ndim)
          ! fill the cell-center values from the updated staggered variables
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
            call phys_face_to_center(ixG^LL,ps1(igrid))
          end do
          !$OMP END PARALLEL DO
        end if
      end if
      ! fix conservation of AMR grid by replacing flux from finer neighbors
      if(associated(temp%sts_handle_errors)) then
        !$OMP PARALLEL DO PRIVATE(igrid)
        do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
          block=>ps(igrid)
          call temp%sts_handle_errors(ps1(igrid)%w,ps1(igrid)%x,ixG^LL,ixM^LL,1)
        end do
        !$OMP END PARALLEL DO
      end if
      if(temp%nflux>temp%nwbc.and.temp%s==1) then
        ! include the changed no-need-ghost-update variables in the last getbc
        type_send_srl=>temp%type_send_srl_sts_2
        type_recv_srl=>temp%type_recv_srl_sts_2
        type_send_r=>temp%type_send_r_sts_2
        type_recv_r=>temp%type_recv_r_sts_2
        type_send_p=>temp%type_send_p_sts_2
        type_recv_p=>temp%type_recv_p_sts_2
        call getbc(global_time,0.d0,ps1,temp%startVar,temp%nflux)
      else
        call getbc(global_time,0.d0,ps1,temp%startwbc,temp%nwbc)
      end if
      !!first step end

      evenstep=.true.

      tmpPs2=>ps1

      do j=2,temp%s
        bj(j)=dble(j**2+j-2)/dble(2*j*(j+1))
        cmu=dble(2*j-1)/dble(j)*bj(j)/bj(j-1)
        cmut=omega1*cmu
        cnu=dble(1-j)/dble(j)*bj(j)/bj(j-2)
        cnut=(bj(j-1)-1.d0)*cmut
        one_mu_nu=1.d0-cmu-cnu
        if(evenstep) then
          tmpPs1=>ps1
          tmpPs2=>ps2
        else
          tmpPs1=>ps2
          tmpPs2=>ps1
        end if

        dtj = cmut*my_dt
        if(stagger_grid) then
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
            ! maybe the following global variables are needed in set_sources
            ! next few lines ensure correct usage of routines like divvector etc
            ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
            block=>ps(igrid)
            ! end maybe the following global variables are needed in set_sources
            call temp%sts_set_sources(ixG^LL,ixM^LL,tmpPs1(igrid)%w,ps(igrid)%x,ps4(igrid)%w,fix_conserve_at_step,dtj,igrid,temp%nflux)
            if(temp%nflux>ndir) then
              tmpPs2(igrid)%w(ixM^T,temp%startVar)=cmu*tmpPs1(igrid)%w(ixM^T,temp%startVar)+&
                cnu*tmpPs2(igrid)%w(ixM^T,temp%startVar)+one_mu_nu*ps(igrid)%w(ixM^T,temp%startVar)+&
                dtj*ps4(igrid)%w(ixM^T,temp%startVar)+cnut*ps3(igrid)%w(ixM^T,temp%startVar)
            end if
            tmpPs2(igrid)%ws(ixC^S,1:nws)=cmu*tmpPs1(igrid)%ws(ixC^S,1:nws)+&
              cnu*tmpPs2(igrid)%ws(ixC^S,1:nws)+one_mu_nu*ps(igrid)%ws(ixC^S,1:nws)+&
              dtj*ps4(igrid)%w(ixC^S,iw_mag(1:nws))+cnut*ps3(igrid)%w(ixC^S,iw_mag(1:nws))
            call phys_face_to_center(ixM^LL,tmpPs2(igrid))
          end do
          !$OMP END PARALLEL DO
        else
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
            ! maybe the following global variables are needed in set_sources
            ! next few lines ensure correct usage of routines like divvector etc
            ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
            block=>ps(igrid)
            ! end maybe the following global variables are needed in set_sources
            call temp%sts_set_sources(ixG^LL,ixM^LL,tmpPs1(igrid)%w,ps(igrid)%x,ps4(igrid)%w,fix_conserve_at_step,dtj,igrid,temp%nflux)
            tmpPs2(igrid)%w(ixM^T,temp%startVar:temp%endVar)=cmu*tmpPs1(igrid)%w(ixM^T,temp%startVar:temp%endVar)+&
              cnu*tmpPs2(igrid)%w(ixM^T,temp%startVar:temp%endVar)+one_mu_nu*ps(igrid)%w(ixM^T,temp%startVar:temp%endVar)+&
              dtj*ps4(igrid)%w(ixM^T,temp%startVar:temp%endVar)+cnut*ps3(igrid)%w(ixM^T,temp%startVar:temp%endVar)
          end do
          !$OMP END PARALLEL DO
        end if
        if(fix_conserve_at_step) then
          call recvflux(1,ndim)
          call sendflux(1,ndim)
          call fix_conserve(tmpPs2,1,ndim,temp%startVar,temp%nflux)
          if(stagger_grid) then
            call fix_edges(tmpPs2,1,ndim)
            ! fill the cell-center values from the updated staggered variables
            !$OMP PARALLEL DO PRIVATE(igrid)
            do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
              call phys_face_to_center(ixG^LL,tmpPs2(igrid))
            end do
            !$OMP END PARALLEL DO
          end if
        end if
        if(associated(temp%sts_handle_errors)) then
        !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
            ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
            block=>ps(igrid)
            call temp%sts_handle_errors(tmpPs2(igrid)%w,ps(igrid)%x,ixG^LL,ixM^LL,j)
          end do
        !$OMP END PARALLEL DO
        end if
        if(temp%nflux>temp%nwbc.and.temp%s==j) then
          ! include the changed no-need-ghost-update variables in the last getbc
          type_send_srl=>temp%type_send_srl_sts_2
          type_recv_srl=>temp%type_recv_srl_sts_2
          type_send_r=>temp%type_send_r_sts_2
          type_recv_r=>temp%type_recv_r_sts_2
          type_send_p=>temp%type_send_p_sts_2
          type_recv_p=>temp%type_recv_p_sts_2
          call getbc(global_time,0.d0,tmpPs2,temp%startVar,temp%nflux)
        else
          call getbc(global_time,0.d0,tmpPs2,temp%startwbc,temp%nwbc)
        end if
        evenstep=.not.evenstep
      end do

      if(associated(temp%sts_after_last_cycle)) then
        !$OMP PARALLEL DO PRIVATE(igrid)
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
          block=>ps(igrid)
          ps(igrid)%w(ixG^T,temp%startVar:temp%endVar)=tmpPs2(igrid)%w(ixG^T,temp%startVar:temp%endVar)
          call temp%sts_after_last_cycle(ixG^LL,ixG^LL,ps(igrid)%w,ps(igrid)%x)
        end do
        !$OMP END PARALLEL DO
        phys_total_energy=total_energy_flag
        prolongprimitive  = prolong_flag
        coarsenprimitive = coarsen_flag
      else
        if(stagger_grid) then
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail; igrid=igrids(iigrid);
            ps(igrid)%w(ixG^T,temp%startVar:temp%endVar)=tmpPs2(igrid)%w(ixG^T,temp%startVar:temp%endVar)
            ps(igrid)%ws=tmpPs2(igrid)%ws
          end do
          !$OMP END PARALLEL DO
        else
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail; igrid=igrids(iigrid);
            ps(igrid)%w(ixG^T,temp%startVar:temp%endVar)=tmpPs2(igrid)%w(ixG^T,temp%startVar:temp%endVar)
          end do
          !$OMP END PARALLEL DO
        end if
      end if

      deallocate(bj)

      if(.not.temp%evolve_magnetic_field) then
        ! restore stagger_grid value
        stagger_grid=stagger_flag
      end if

      temp=>temp%next
    end do

    if(associated(head_sts_terms)) then
      ! point bc mpi data type back to full type for (M)HD
      type_send_srl=>type_send_srl_f
      type_recv_srl=>type_recv_srl_f
      type_send_r=>type_send_r_f
      type_recv_r=>type_recv_r_f
      type_send_p=>type_send_p_f
      type_recv_p=>type_recv_p_f
    end if

    bcphys=.true.

    if(phys_partial_ionization) then
      ! update temperature variable in w
      !$OMP PARALLEL DO PRIVATE(igrid)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        call phys_update_temperature(ixG^LL,ixG^LL,ps(igrid)%w,ps(igrid)%w,ps(igrid)%x)
      end do
      !$OMP END PARALLEL DO
    end if

  end subroutine sts_add_source1

end module mod_supertimestepping
