module mod_usr
  use mod_global_parameters
  use mod_srmhd_parameters
  use mod_srmhd  
  implicit none
  double precision :: rhoRight, pRight, v1Right, v2Right, v3Right, &
                        b1Right,b2Right,b3Right, &
                        rhoLeft, pLeft, v1Left, v2Left, v3Left, &
                        b1Left,b2Left,b3Left, &
                        xsplit1
  character(len=30):: coordinate_system
contains

  !> Read this module s parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ rhoRight, pRight, v1Right, v2Right, v3Right, &
                        b1Right,b2Right,b3Right, &
                        rhoLeft, pLeft, v1Left, v2Left, v3Left, &
                        b1Left,b2Left,b3Left, &
                        xsplit1,coordinate_system


    write(*,*)'Reading usr_list'
    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine usr_params_read

  subroutine usr_init()
    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    call usr_params_read(par_files)
    call set_coordinate_system(trim(coordinate_system))
    call srmhd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_variables

    call usr_params_read(par_files) 
    if(dabs(b1Left-b1Right)>smalldouble)then
      call mpistop(' it should be b1left=b1Rigth')
    end if
  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ixO^L,w,x)
    ! initialize one grid 
    integer, intent(in) :: ixG^L,ixO^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)


    logical,save     :: first=.true.
    double precision :: wold(ixG^S,1:nw) 



    where(x(ixG^S,1)<xsplit1)
       w(ixG^S,rho_)     = rholeft
       w(ixG^S,mom(1))   = v1left
       w(ixG^S,mom(2))   = v2left
       w(ixG^S,mom(3))   = v3left
       w(ixG^S,p_ )      = pleft
       w(ixG^S,mag(1) )  = b1left
       w(ixG^S,mag(2) )  = b2left
       w(ixG^S,mag(3) )  = b3left
    elsewhere
       w(ixG^S,rho_)     = rhoright
       w(ixG^S,mom(1))   = v1right
       w(ixG^S,mom(2))   = v2right
       w(ixG^S,mom(3))   = v3right
       w(ixG^S,p_  )     = pright
       w(ixG^S,mag(1) )  = b1right
       w(ixG^S,mag(2) )  = b2right
       w(ixG^S,mag(3) )  = b3right
    endwhere
    if (srmhd_glm) w(ixG^S,psi_)     = 0.0d0



    call srmhd_get_4u_from_3v(ixG^L,ixG^L,w(ixG^S,mom(1):mom(ndir)))
    call srmhd_to_conserved(ixG^L,ixG^L,w,x)



    wold=w

 !    PRINT*,' well see',w(ixG^S,lfac_),mom(1),mag(1)

    call srmhd_to_primitive(ixG^L,ixG^L,w,x)
! print*,'byyye' ,saveigrid 
    call srmhd_to_conserved(ixG^L,ixG^L,w,x)
!    if(maxval(dabs(w(ixG^S,1:nwflux)-wold(ixG^S,1:nwflux)))>1d-7)then
!      PRINT*,' the max user', (w(ixG^S,1)-wold(ixG^S,1)),all(x(ixG^S,1)<xsplit1)
!      PRINT*,' the max is ', maxval(dabs(w(ixG^S,1:nwflux)-wold(ixG^S,1:nwflux)))
!      PRINT*,'the position', maxloc(dabs(w(ixG^S,1:nwflux)-wold(ixG^S,1:nwflux)))
!PRINT*,' is you',w(ixG^S,lfac_),wold(ixG^S,lfac_)
!      call mpistop('is wrong ')
!    end if


  end subroutine initonegrid_usr

end module mod_usr
