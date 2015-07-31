!=============================================================================
subroutine set_B0_grid(igrid)

include 'amrvacdef.f'

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------

call set_B0_cell(pB0_cell(igrid)%w,px(igrid)%x,ixGlo1,ixGhi1)
call set_B0_face(igrid,px(igrid)%x)

end subroutine set_B0_grid
!=============================================================================
subroutine set_B0_cell(wB0,x,ixmin1,ixmax1)

include 'amrvacdef.f'

double precision, intent(inout) :: wB0(ixmin1:ixmax1,1:ndir)
double precision, intent(in) :: x(ixGlo1:ixGhi1,1:ndim)
integer, intent(in):: ixmin1,ixmax1
!-----------------------------------------------------------------------------
wB0(ixmin1:ixmax1,1:ndir)=zero

! approximate cell-averaged B0 as cell-centered B0
select case (typeaxial)
case ("spherical")
   
end select

if(dabs(Busr)/=zero) then
   call specialset_B0(ixmin1,ixmax1,ixmin1,ixmax1,x,wB0)
end if

end subroutine set_B0_cell
!=============================================================================
subroutine set_B0_face(igrid,x)

include 'amrvacdef.f'

integer, intent(in) :: igrid
double precision, intent(in) :: x(ixGlo1:ixGhi1,1:ndim)

double precision :: xC(ixGlo1:ixGhi1,1:ndim),dx1,xmin1,xshift1
integer :: idims, ixCmin1,ixCmax1, ix, idims2
!-----------------------------------------------------------------------------
dx1=rnode(rpdx1_,igrid);
xmin1=rnode(rpxmin1_,igrid);

do idims=1,ndim
   ixCmin1=ixMlo1-kr(1,idims); ixCmax1=ixMhi1;
   xshift1=half*(one-kr(1,idims));
   do idims2=1,ndim
      select case(idims2)
      case(1)
        do ix = ixCmin1,ixCmax1
          xC(ix,1)=xmin1+(dble(ix-dixB)-xshift1)*dx1
        end do
      end select
   end do

   select case(idims)
   case (1)
      call set_B0_cell(pB0_face1(igrid)%w,xC,ixCmin1,ixCmax1) 
   end select
end do

end subroutine set_B0_face
!=============================================================================
subroutine alloc_B0_grid(igrid)

include 'amrvacdef.f'

integer, intent(in) :: igrid

integer :: ixCmin1,ixCmax1
!-----------------------------------------------------------------------------

allocate(pB0_cell(igrid)%w(ixGlo1:ixGhi1,1:ndir))
ixCmin1=ixMlo1-kr(1,1); ixCmax1=ixMhi1;
allocate(pB0_face1(igrid)%w(ixCmin1:ixCmax1,1:ndir))

end subroutine alloc_B0_grid
!=============================================================================
subroutine dealloc_B0_grid(igrid)

include 'amrvacdef.f'

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------

deallocate(pB0_cell(igrid)%w)
deallocate(pB0_face1(igrid)%w)

end subroutine dealloc_B0_grid
!=============================================================================
