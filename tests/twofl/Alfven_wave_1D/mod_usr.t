module mod_usr
  use mod_twofl


  implicit none

  ! background
  double precision :: pn0, rhon0, pc0, rhoc0, b0

  !!user defined
  double precision:: ampl=1d-3
  double precision :: nn = 10
  type arrayptr
    double precision,dimension(:), pointer :: p
    character(len=30) ::  namevar
  end type arrayptr

contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ ampl,nn

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine usr_params_read





  subroutine usr_init()
    use mod_variables
    
    ! default units = 1
    ! background
    pn0 = 2d1    
    pc0 = 1d1 
    rhon0 = 2d1   
    rhoc0 = 1d1   
    b0 = 0.4

    call usr_params_read(par_files)

    if(mype .eq. 0) then
      print*, "Amplitude ", ampl
      print*, "NX= ", nn
    endif

    usr_init_one_grid => initonegrid_usr
    usr_set_equi_vars => special_set_equi_vars
    usr_set_B0              => specialset_B0
    usr_set_J0              => specialset_J0

    usr_set_parameters  => init_params_usr
  
    call set_coordinate_system("Cartesian_1.75D")

    call twofl_activate()


  end subroutine usr_init
  
    !> Here one can add a steady (time-independent) equi vars
    subroutine special_set_equi_vars(ixI^L,ixO^L,x,w0)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L,ixO^L
      double precision, intent(in)    :: x(ixI^S,1:ndim)
      double precision, intent(inout) :: w0(ixI^S,1:number_equi_vars)

      w0(ixO^S,equi_pe_n0_) = pn0
      w0(ixO^S,equi_pe_c0_) = pc0
      w0(ixO^S,equi_rho_n0_) = rhon0
      w0(ixO^S,equi_rho_c0_) = rhoc0


    end subroutine special_set_equi_vars


  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here add a time-independent background magnetic field
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)
      wB0(ixO^S,1)=b0
      wB0(ixO^S,2:3)=0d0

  end subroutine specialset_B0

  subroutine specialset_J0(ixI^L,ixO^L,x,wJ0)
  ! Here add a time-independent background current density 
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wJ0(ixI^S,7-2*ndir:ndir)
    wJ0(ixO^S,:)=zero


  end subroutine specialset_J0




  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    complex, parameter :: ic = dcmplx(0,1)
    double complex :: wave(ixO^S), B1, Vc, Vn,omega
    double precision :: va0, k


    k = 2*dpi*nn/(xprobmax1-xprobmin1)

    va0 = b0/dsqrt(rhoc0)

    if(nn .ne. 5) then
      if(mype .eq. 0) print* , "nn=",nn, &
          " and VALUES INTRODUCED only for nn=5",&
          " and twofl_alpha_coll 1d-3,1d-2,1d-1,1d0,1d2,1d3,1d4.",&
          " twofl_alpha_coll set to 0"
      twofl_alpha_coll = 0d0
    else
      if(twofl_alpha_coll == 1d-3) then
        omega = dcmplx(2.483586707670307,0.009999837878226345)
      elseif(twofl_alpha_coll == 1d-2) then
        omega = dcmplx(2.4776002530940726, 0.09983709383324003)
      elseif(twofl_alpha_coll == 1d-1) then
        omega = dcmplx(1.8476978931291896,0.7139437294531664)
      elseif(twofl_alpha_coll == 1d0) then
        omega = dcmplx(1.435586475340743, 0.06869586285089642)
      elseif(twofl_alpha_coll == 1d2) then
        omega = dcmplx(1.4339344875543152, 0.0006853893715324451)
      elseif(twofl_alpha_coll == 1d3) then
        omega = dcmplx(1.4339343253916603,0.00006853892165121452)
      elseif(twofl_alpha_coll == 1d4) then
        omega = dcmplx(1.4339343237700348,6.853892149619427D-6)
      else
        if(mype .eq. 0) print* , "twofl_alpha_coll=",twofl_alpha_coll, &
          " and VALUES INTRODUCED only for nn=5",&
          " and twofl_alpha_coll 1d-3,1d-2,1d-1,1d0,1d2,1d3,1d4.",&
          " twofl_alpha_coll set to 0"
        twofl_alpha_coll = 0d0
      endif  
    endif

    if(twofl_alpha_coll == 0d0) then
      if(mype .eq. 0) print*, "no coupling"
      omega = va0 * k
    endif


    Vc = ampl * va0
    Vn = Vc * (twofl_alpha_coll * rhon0 *rhoc0)/(ic * rhon0 * omega + twofl_alpha_coll * rhon0 *rhoc0)
    B1 = - k * Vc * b0 / omega

    if(mype .eq. 0) then
      print*, "Vc, Vn, B1", Vc, Vn, B1
    endif 
    wave(ixO^S) = exp(ic *  (- k*(x(ixO^S,1) -xprobmin1))) 
    w(ixO^S,mom_n(3)) = real(Vn * wave(ixO^S))
    w(ixO^S,mom_c(3)) = real(Vc * wave(ixO^S))
    w(ixO^S,mag(3)) = real(B1 * wave(ixO^S))
    call twofl_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine initonegrid_usr

  


  subroutine init_params_usr()
  end subroutine init_params_usr





  function rargsort(a) result(b)
    ! Returns the indices that would sort an array.
    !
    ! Arguments
    ! ---------
    !
    double precision, intent(in):: a(:)   ! array of numbers
    integer :: b(size(a))         ! indices into the array 'a' that sort it
    !
    ! Example
    ! -------
    !
    ! rargsort([4.1_dp, 2.1_dp, 2.05_dp, -1.5_dp, 4.2_dp]) ! Returns [4, 3, 2, 1, 5]
    
    integer :: N                           ! number of numbers/vectors
    integer :: i,imin                      ! indices: i, i of smallest
    integer :: temp1                       ! temporary
    double precision:: temp2
    double precision:: a2(size(a))
    a2 = a
    N=size(a)
    do i = 1, N
        b(i) = i
    end do
    do i = 1, N-1
        ! find ith smallest in 'a'
        imin = minloc(a2(i:),1) + i - 1
        ! swap to position i in 'a' and 'b', if not already there
        if (imin /= i) then
            temp2 = a2(i); a2(i) = a2(imin); a2(imin) = temp2
            temp1 = b(i); b(i) = b(imin); b(imin) = temp1
        end if
    end do
  end function






end module mod_usr
