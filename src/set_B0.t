subroutine set_B0_grid(igrid)
  use mod_global_parameters

  integer, intent(in) :: igrid

  call set_B0_cell(pw(igrid)%B0(:^D&,:,0),pw(igrid)%x,ixG^LL,ixG^LL)
  call set_J0_cell(igrid,pw(igrid)%J0,ixG^LL,ixM^LL^LADD1)
  call set_B0_face(igrid,pw(igrid)%x,ixG^LL,ixM^LL)

end subroutine set_B0_grid

subroutine set_B0_cell(wB0,x,ixI^L,ix^L)
  use mod_usr_methods, only: usr_set_B0
  use mod_global_parameters
  
  integer, intent(in):: ixI^L,ix^L
  double precision, intent(inout) :: wB0(ixI^S,1:ndir)
  double precision, intent(in) :: x(ixI^S,1:ndim)

  wB0(ix^S,1:ndir)=zero

  ! approximate cell-averaged B0 as cell-centered B0
  select case (typeaxial)
  case ("spherical")
     {^NOONED
     if (dabs(Bdip)>smalldouble) then
        wB0(ix^S,1)=2.0d0*Bdip*dcos(x(ix^S,2))/x(ix^S,1)**3
        wB0(ix^S,2)=Bdip*dsin(x(ix^S,2))/x(ix^S,1)**3
     end if
  
     if (abs(Bquad)>smalldouble) then
        wB0(ix^S,1)=wB0(ix^S,1) &
             +Bquad*0.5d0*(1.0d0+3.0d0*dcos(2.0d0*x(ix^S,2)))/x(ix^S,1)**4
        wB0(ix^S,2)=wB0(ix^S,2)+Bquad*dsin(2.0d0*x(ix^S,2))/x(ix^S,1)**4
     end if
     if (abs(Boct)>smalldouble) then
        wB0(ix^S,1)=wB0(ix^S,1) &
                     +Boct*(10.0d0*dcos(2.0d0*x(ix^S,2))-2.0d0) &
                          *dcos(x(ix^S,2))/x(ix^S,1)**5
        wB0(ix^S,2)=wB0(ix^S,2) &
                     +Boct*1.5d0*(3.0d0+5.0d0*dcos(2.0d0*x(ix^S,2))) &
                          *dsin(x(ix^S,2))/x(ix^S,1)**5
     end if
     }
  end select
  
  if (associated(usr_set_B0)) call usr_set_B0(ixI^L,ix^L,x,wB0)

end subroutine set_B0_cell

subroutine set_J0_cell(igrid,wJ0,ixI^L,ix^L)
  use mod_usr_methods, only: usr_set_J0
  use mod_global_parameters

  integer, intent(in):: igrid,ixI^L,ix^L
  double precision, intent(inout) :: wJ0(ixI^S,7-2*ndir:3)
  integer :: idirmin0, idirmin

  if(associated(usr_set_J0)) then
    call usr_set_J0(ixI^L,ix^L,pw(igrid)%x,wJ0)
  else
    idirmin0 = 7-2*ndir
    call curlvector(pw(igrid)%B0(:^D&,:,0),ixI^L,ix^L,wJ0,idirmin,idirmin0,ndir)
  end if

end subroutine set_J0_cell

subroutine set_B0_face(igrid,x,ixI^L,ix^L)
  use mod_global_parameters

  integer, intent(in) :: igrid, ixI^L, ix^L
  double precision, intent(in) :: x(ixI^S,1:ndim)

  double precision :: xC(ixI^S,1:ndim),dxu(ndim),xshift(ndim)
  integer :: idims, ixC^L, idims2

  ^D&dxu(^D)=rnode(rpdx^D_,igrid);

  do idims=1,ndim
     ixCmin^D=ixmin^D-kr(^D,idims); ixCmax^D=ixmax^D;
     xshift(1:ndim)=half*(one-kr(1:ndim,idims))
     do idims2=1,ndim
       xC(ixC^S,idims2)=x(ixC^S,idims2)+(half-xshift(idims2))*dxu(idims2)
     end do
     if(stretched_grid) then
       if(slab_stretched.and.idims==^ND) then
         xC(ixC^S,^ND)=x(ixC^S,^ND)+0.5d0*pw(igrid)%dx(ixC^S,^ND)
       else if(.not.slab_stretched.and.idims==1) then
         xC(ixC^S,1)=x(ixC^S,1)+0.5d0*pw(igrid)%dx(ixC^S,1)
       end if
     end if
     call set_B0_cell(pw(igrid)%B0(:^D&,:,idims),xC,ixI^L,ixC^L)
  end do

end subroutine set_B0_face

subroutine alloc_B0_grid(igrid)
  use mod_global_parameters

  integer, intent(in) :: igrid

  if(.not. allocated(pw(igrid)%B0)) then
    allocate(pw(igrid)%B0(ixG^T,1:ndir,0:ndim))
    allocate(pw(igrid)%J0(ixG^T,7-2*ndir:3))
  end if

end subroutine alloc_B0_grid

subroutine dealloc_B0_grid(igrid)
  use mod_global_parameters

  integer, intent(in) :: igrid

  deallocate(pw(igrid)%B0)
  deallocate(pw(igrid)%J0)

end subroutine dealloc_B0_grid
