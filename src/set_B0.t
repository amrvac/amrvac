!=============================================================================
subroutine set_B0_grid(igrid)

include 'amrvacdef.f'

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------

call set_B0_cell(pB0_cell(igrid)%w,px(igrid)%x,ixG^LL)
call set_B0_face(igrid,px(igrid)%x)

end subroutine set_B0_grid
!=============================================================================
subroutine set_B0_cell(wB0,x,ix^L)

include 'amrvacdef.f'

double precision, intent(inout) :: wB0(ix^S,1:ndir)
double precision, intent(in) :: x(ixG^T,1:ndim)
integer, intent(in):: ix^L
!-----------------------------------------------------------------------------
wB0(ix^S,1:ndir)=zero

! approximate cell-averaged B0 as cell-centered B0
select case (typeaxial)
case ("spherical")
   {^NOONED
   if (abs(Bdip)>smalldouble) then
      wB0(ix^S,1)=2.0d0*Bdip*cos(x(ix^S,2))/x(ix^S,1)**3.0d0
      wB0(ix^S,2)=Bdip*sin(x(ix^S,2))/x(ix^S,1)**3.0d0
   end if

   if (abs(Bquad)>smalldouble) then
      wB0(ix^S,1)=wB0(ix^S,1) &
           +Bquad*0.5d0*(1.0d0+3.0d0*cos(2.0d0*x(ix^S,2)))/x(ix^S,1)**4
      wB0(ix^S,2)=wB0(ix^S,2)+Bquad*sin(2.0d0*x(ix^S,2))/x(ix^S,1)**4
   end if
   if (abs(Boct)>smalldouble) then
      wB0(ix^S,1)=wB0(ix^S,1) &
                   +Boct*(10.0d0*cos(2.0d0*x(ix^S,2))-2.0d0) &
                        *cos(x(ix^S,2))/x(ix^S,1)**5
      wB0(ix^S,2)=wB0(ix^S,2) &
                   +Boct*1.5d0*(3.0d0+5.0d0*cos(2.0d0*x(ix^S,2))) &
                        *sin(x(ix^S,2))/x(ix^S,1)**5
   end if
   }
end select

if(dabs(Busr)/=zero) then
   call specialset_B0(ix^L,ix^L,x,wB0)
end if

end subroutine set_B0_cell
!=============================================================================
subroutine set_B0_face(igrid,x)

include 'amrvacdef.f'

integer, intent(in) :: igrid
double precision, intent(in) :: x(ixG^T,1:ndim)

double precision :: xC(ixG^T,1:ndim),dx^D,xmin^D,xshift^D
integer :: idims, ixC^L, ix, idims2
!-----------------------------------------------------------------------------
dx^D=rnode(rpdx^D_,igrid);
xmin^D=rnode(rpxmin^D_,igrid);

do idims=1,ndim
   ixCmin^D=ixMlo^D-kr(^D,idims); ixCmax^D=ixMhi^D;
   xshift^D=half*(one-kr(^D,idims));
   do idims2=1,ndim
      select case(idims2)
      {case(^D)
        do ix = ixC^LIM^D
          xC(ix^D%ixC^S,^D)=xmin^D+(dble(ix-dixB)-xshift^D)*dx^D
        end do\}
      end select
   end do

   select case(idims)
   {case (^D)
      call set_B0_cell(pB0_face^D(igrid)%w,xC,ixC^L) \}
   end select
end do

end subroutine set_B0_face
!=============================================================================
subroutine alloc_B0_grid(igrid)

include 'amrvacdef.f'

integer, intent(in) :: igrid

integer :: ixC^L
!-----------------------------------------------------------------------------

allocate(pB0_cell(igrid)%w(ixG^T,1:ndir))
{ixCmin^DD=ixMlo^DD-kr(^DD,^D); ixCmax^DD=ixMhi^DD;
allocate(pB0_face^D(igrid)%w(ixC^S,1:ndir))\}

end subroutine alloc_B0_grid
!=============================================================================
subroutine dealloc_B0_grid(igrid)

include 'amrvacdef.f'

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------

deallocate(pB0_cell(igrid)%w)
{deallocate(pB0_face^D(igrid)%w)\}

end subroutine dealloc_B0_grid
!=============================================================================
