!=============================================================================
subroutine set_B0_grid(igrid)
use mod_global_parameters

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------

call set_B0_cell(pB0_cell(igrid)%w,px(igrid)%x,ixG^LL,ixG^LL)
{#IFDEF FCT
call set_B0_face(igrid,px(igrid)%x,ixG^LL,ixM^LL^LADD1)
}{#IFNDEF FCT
call set_B0_face(igrid,px(igrid)%x,ixG^LL,ixM^LL)
}

end subroutine set_B0_grid
!=============================================================================
subroutine set_B0_cell(wB0,x,ixI^L,ix^L)
use mod_usr_methods, only: usr_set_B0
use mod_global_parameters

integer, intent(in):: ixI^L,ix^L
double precision, intent(inout) :: wB0(ixI^S,1:ndir)
double precision, intent(in) :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------
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

if (dabs(Busr)/=zero) then
   if (.not. associated(usr_set_B0)) then
      call mpistop("usr_set_B0 not defined")
   else
      call usr_set_B0(ixI^L,ix^L,x,wB0)
   end if
end if

end subroutine set_B0_cell
!=============================================================================
subroutine set_B0_face(igrid,x,ixI^L,ix^L)

use mod_global_parameters

integer, intent(in) :: igrid, ixI^L, ix^L
double precision, intent(in) :: x(ixI^S,1:ndim)

double precision :: xC(ixI^S,1:ndim),dx^D,xmin^D,xshift^D
integer :: idims, ixC^L, ix, idims2
!-----------------------------------------------------------------------------
dx^D=rnode(rpdx^D_,igrid);
xmin^D=rnode(rpxmin^D_,igrid);
{#IFDEF STRETCHGRID
logG=logGs(node(plevel_,igrid))
qst=qsts(node(plevel_,igrid))
}

do idims=1,ndim
   ixCmin^D=ixmin^D-kr(^D,idims); ixCmax^D=ixmax^D;
   xshift^D=half*(one-kr(^D,idims));
   do idims2=1,ndim
      select case(idims2)
      {case(^D)
        do ix = ixC^LIM^D
          xC(ix^D%ixC^S,^D)=xmin^D+(dble(ix-dixB)-xshift^D)*dx^D
        end do\}
      end select
   end do
{#IFDEF STRETCHGRID
   do ix = ixCmin1,ixCmax1
     if(xshift1==0.d0) then
       xC(ix^%1ixC^S,1)=xmin1*qst**(ix-dixB)
     else
       xC(ix^%1ixC^S,1)=xmin1/(one-half*logG)*qst**(ix-dixB-1)
     end if
   end do
}

   select case(idims)
   {case (^D)
      call set_B0_cell(pB0_face^D(igrid)%w,xC,ixI^L,ixC^L) \}
   end select
end do

end subroutine set_B0_face
!=============================================================================
subroutine alloc_B0_grid(igrid)

use mod_global_parameters

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------

allocate(pB0_cell(igrid)%w(ixG^T,1:ndir))
{allocate(pB0_face^D(igrid)%w(ixG^T,1:ndir))\}

end subroutine alloc_B0_grid
!=============================================================================
subroutine dealloc_B0_grid(igrid)

use mod_global_parameters

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------

deallocate(pB0_cell(igrid)%w)
{deallocate(pB0_face^D(igrid)%w)\}

end subroutine dealloc_B0_grid
!=============================================================================
