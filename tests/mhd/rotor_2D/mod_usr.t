module mod_usr
  use mod_mhd
  implicit none

contains

  subroutine usr_init()
    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 
    usr_set_B0        => specialset_B0
    usr_init_vector_potential=>initvecpot_usr

    call mhd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: rhodisk,rr0,rr1,x1disk,x2disk,v0
    double precision :: r(ixI^S),fslope(ixI^S),rinv(ixI^S)
    double precision :: Bloc(ixI^S,1:ndir)
    logical, save :: first=.true.

    rhodisk=10.0d0
    rr0=0.1d0
    rr1=0.115d0
    x1disk=xprobmin1+half*(xprobmax1-xprobmin1)
    x2disk=xprobmin2+half*(xprobmax2-xprobmin2)
    v0=2.0d0
    if(first .and. mype==0)then
      write(*,*)'Doing 2D ideal MHD, rotor problem'
      write(*,*)'rr0,rr1:'
      write(*,*)rr0,rr1
      write(*,*)'gamma,v0:'
      write(*,*)mhd_gamma,v0
      first=.false.
    endif
    w(ixO^S,p_)=one
    if(B0field) then
      w(ixO^S,mag(:))=0.d0
    else if(stagger_grid) then
      call b_from_vector_potential(block%ixGs^L,ixI^L,ixO^L,block%ws,x)
      call mhd_face_to_center(ixO^L,block)
    else
      call get_B(ixI^L,ixO^L,Bloc,x)
      w(ixO^S,mag(:))=Bloc(ixO^S,:)
    end if
    r(ixO^S)=dsqrt((x(ixO^S,1)-x1disk)**2+(x(ixO^S,2)-x2disk)**2)
    where(r(ixO^S)>rr0)
      fslope(ixO^S)=(rr1-r(ixO^S))/(rr1-rr0)
      rinv(ixO^S)=one/r(ixO^S)
    elsewhere
      fslope(ixO^S)=one
      rinv(ixO^S)=one/rr0
    endwhere
    where(r(ixO^S)<rr1)
      w(ixO^S,rho_)=one+(rhodisk-one)*fslope(ixO^S)
      w(ixO^S,mom(1))=-v0*fslope(ixO^S)*(x(ixO^S,2)-x2disk)*rinv(ixO^S)
      w(ixO^S,mom(2))=v0*fslope(ixO^S)*(x(ixO^S,1)-x1disk)*rinv(ixO^S)
    elsewhere
      w(ixO^S,rho_)=one
      w(ixO^S,mom(1))=zero
      w(ixO^S,mom(2))=zero
    endwhere

    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixC^L,idir
    double precision, intent(in)       :: xC(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)

    if (idir==3) then
      A(ixC^S) =5.0d0/dsqrt(4.0d0*dpi)*xC(ixC^S,2)
    else
      A(ixC^S) = 0.d0
    end if

  end subroutine initvecpot_usr

  subroutine get_B(ixI^L,ixO^L,B,x)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out) :: B(ixI^S,1:ndir)

    B(ixO^S,1)=5.0d0/dsqrt(4.0d0*dpi)
    B(ixO^S,2)=0.d0

  end subroutine get_B

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    ! this subroutine can be used in convert, to add auxiliary variables to the
    ! converted output file, for further analysis using tecplot, paraview, ....
    ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
    !
    ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with 
    ! corresponding normalization values (default value 1)
    use mod_mhd_phys
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    ! .. local ..
    double precision:: divb(ixI^S)
    double precision :: current(ixI^S,7-2*ndir:3)
    integer          :: idirmin
    
    call get_divb(w,ixI^L,ixO^L,divb)
    w(ixO^S,nw+1)=divb(ixO^S)

    call get_normalized_divb(w,ixI^L,ixO^L,divb)
    w(ixO^S,nw+2)=divb(ixO^S)
    
    call get_current(w,ixI^L,ixO^L,idirmin,current)
    w(ixO^S,nw+3)=current(ixO^S,3)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames
    varnames='divb divbnorm jz'

  end subroutine specialvarnames_output

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here one can add a steady (time-independent) potential background field
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    double precision :: Bloc(ixI^S,1:ndir)

    call get_B(ixI^L,ixO^L,Bloc,x)
    wB0(ixO^S,:)=wB0(ixO^S,:)+Bloc(ixO^S,:)

  end subroutine specialset_B0

end module mod_usr
