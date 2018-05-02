module mod_usr
  use mod_mhd
  implicit none

contains

  subroutine usr_init()
    usr_init_one_grid => initonegrid_usr

    call set_coordinate_system('Cartesian')
    call mhd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    ! initialize one grid
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision:: mpol,k1,Rjet,width,qv,B0,dv
    double precision, dimension(ixG^S) :: r,phi
    logical, save :: first=.true.

    ! KH in 3D, Keppens & Toth
    r(ixG^S)=dsqrt(x(ixG^S,2)**2+x(ixG^S,3)**2)
    phi(ixG^S)=datan2(x(ixG^S,3),x(ixG^S,2))

    mpol=one
    k1=two*dpi/(xprobmax1-xprobmin1)
    Rjet=0.5d0
    width=0.1d0*Rjet
    qv=0.645d0
    B0=0.129d0
    dv=0.01d0

    w(ix^S,rho_)=one
    w(ix^S,p_)=one
    w(ix^S,mag(1))=B0
    w(ix^S,mag(2))=zero
    w(ix^S,mag(3))=zero
    w(ix^S,mom(1))=qv*dtanh((r(ix^S)-Rjet)/width)
    w(ix^S,mom(2))=dv*dexp(-((r(ix^S)-Rjet)/(4.0d0*width))**2) &
                   *dcos(mpol*phi(ix^S))*dsin(k1*x(ix^S,1))*dcos(phi(ix^S))
    w(ix^S,mom(3))=dv*dexp(-((r(ix^S)-Rjet)/(4.0d0*width))**2) &
                   *dcos(mpol*phi(ix^S))*dsin(k1*x(ix^S,1))*dsin(phi(ix^S))

    call mhd_to_conserved(ixG^L,ix^L,w,x)

    if(first)then
      if(mype==0)then
       write(*,*)'Doing 3D ideal MHD, KH jet problem'
      endif
      first=.false.
    endif

  end subroutine initonegrid_usr

end module mod_usr
