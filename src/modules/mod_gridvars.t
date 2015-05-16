{#IFDEF PARTICLES
module mod_gridvars
!
! Contains gridvars to be used with particles
!
use mod_physicaldata
implicit none

type(walloc), save,dimension(ngridshi) :: gridvars
type(walloc), save,dimension(ngridshi) :: gridvars_old

! I use this variables to set the current igrid and ipart for the particle integrator:
integer                                :: igrid_working, ipart_working

{#IFDEF PARTICLES_LORENTZ
integer,parameter                      :: ngridvars=2*^NC
integer,parameter                      :: {bp^C_=^C},{^IFSRMHD {up^C_=bp^NC_+^C} } {^IFMHD {ep^C_=bp^NC_+^C} }
}
{#IFDEF PARTICLES_ADVECT
integer,parameter                      :: ngridvars=^NC
integer,parameter                      :: {vp^C_=^C}
}
{#IFDEF PARTICLES_GCA
integer,parameter                      :: ngridvars=7*^NC
integer,parameter                      :: {bp^C_=^C},{ep^C_=bp^NC_+^C},{grad_kappa_B^C_=ep^NC_+^C}
integer,parameter                      :: {b_dot_grad_b^C_=grad_kappa_B^NC_+^C}
integer,parameter                      :: {ue_dot_grad_b^C_=b_dot_grad_b^NC_+^C}
integer,parameter                      :: {b_dot_grad_ue^C_=ue_dot_grad_b^NC_+^C}
integer,parameter                      :: {ue_dot_grad_ue^C_=b_dot_grad_ue^NC_+^C}
}
{#IFDEF PARTICLES_USER
! specify parameters here
}

contains

!=============================================================================
subroutine init_gridvars()

include 'amrvacdef.f'
integer                                   :: igrid, iigrid
!-----------------------------------------------------------------------------

do igrid=1,ngridshi
   nullify(gridvars(igrid)%w)
   if (time_advance) nullify(gridvars_old(igrid)%w)
end do

do iigrid=1,igridstail; igrid=igrids(iigrid);
   allocate(gridvars(igrid)%w(ixG^T,1:ngridvars))
   if (time_advance)  allocate(gridvars_old(igrid)%w(ixG^T,1:ngridvars))
end do

! now fill the gridvars:
call fill_gridvars(gridvars,pw)
if (time_advance) call fill_gridvars(gridvars_old,pwold)

end subroutine init_gridvars
!=============================================================================
subroutine finish_gridvars()

include 'amrvacdef.f'

integer             :: iigrid, igrid
!-----------------------------------------------------------------------------
do iigrid=1,igridstail; igrid=igrids(iigrid);

   deallocate(gridvars(igrid)%w)
   if (time_advance) deallocate(gridvars_old(igrid)%w)

end do

end subroutine finish_gridvars
!=============================================================================
subroutine interpolate_var(igrid,ixI^L,ixO^L,gf,x,xloc,gfloc)

include 'amrvacdef.f'
integer, intent(in)                   :: igrid,ixI^L, ixO^L
double precision, intent(in)          :: gf(ixI^S)
double precision, intent(in)          :: x(ixI^S,1:ndim)
double precision, intent(in)          :: xloc(1:ndir)
double precision, intent(out)         :: gfloc
integer                               :: ic^D, ic1^D, ic2^D, idir
double precision                      :: xd^D
{#IFDEF D2
double precision                      :: c00, c10
}
{#IFDEF D3
double precision                      :: c0, c1, c00, c10, c01, c11
}
character(len=1024)                   :: line
!-----------------------------------------------------------------------------

! flat interpolation:
{ic^D = int((xloc(^D)-rnode(rpxmin^D_,igrid))/rnode(rpdx^D_,igrid)) + 1 + dixB \}
!gfloc = gf(ic^D)


! linear interpolation:
{
if (x({ic^DD},^D) .lt. xloc(^D)) then
   ic1^D = ic^D
else
   ic1^D = ic^D -1
end if
ic2^D = ic1^D + 1
\}

{^D& 
if (ic1^D.lt.ixGlo^D+1 .or. ic2^D.gt.ixGhi^D-1) then
    line = ''
    write(line,"(a)")'Trying to interpolate from out of grid!'
    write(line,"(a,a,i3.2)")trim(line),' direction: ',^D
    write(line,"(a,a,^NDes14.6)")trim(line),' position: ',xloc(1:ndim)
    write(line,"(a,a,^NDi3.2,^NDi3.2)"),trim(line),' indices:', ic1^D,ic2^D
    call mpistop(line)
end if
\}


{#IFDEF D1
xd1 = (xloc(1)-x(ic11,1)) / (x(ic21,1) - x(ic11,1))
gfloc  = gf(ic11) * (1.0d0 - xd1) + gf(ic21) * xd1
}
{#IFDEF D2
xd1 = (xloc(1)-x(ic11,ic12,1)) / (x(ic21,ic12,1) - x(ic11,ic12,1))      
xd2 = (xloc(2)-x(ic11,ic12,2)) / (x(ic11,ic22,2) - x(ic11,ic12,2))
c00 = gf(ic11,ic12) * (1.0d0 - xd1) + gf(ic21,ic12) * xd1
c10 = gf(ic11,ic22) * (1.0d0 - xd1) + gf(ic21,ic22) * xd1
gfloc  = c00 * (1.0d0 - xd2) + c10 * xd2
}
{#IFDEF D3
xd1 = (xloc(1)-x(ic11,ic12,ic13,1)) / (x(ic21,ic12,ic13,1) - x(ic11,ic12,ic13,1))      
xd2 = (xloc(2)-x(ic11,ic12,ic13,2)) / (x(ic11,ic22,ic13,2) - x(ic11,ic12,ic13,2))      
xd3 = (xloc(3)-x(ic11,ic12,ic13,3)) / (x(ic11,ic12,ic23,3) - x(ic11,ic12,ic13,3))    

c00 = gf(ic11,ic12,ic13) * (1.0d0 - xd1) + gf(ic21,ic12,ic13) * xd1
c10 = gf(ic11,ic22,ic13) * (1.0d0 - xd1) + gf(ic21,ic22,ic13) * xd1
c01 = gf(ic11,ic12,ic23) * (1.0d0 - xd1) + gf(ic21,ic12,ic23) * xd1
c11 = gf(ic11,ic22,ic23) * (1.0d0 - xd1) + gf(ic21,ic22,ic23) * xd1

c0  = c00 * (1.0d0 - xd2) + c10 * xd2
c1  = c01 * (1.0d0 - xd2) + c11 * xd2

gfloc = c0 * (1.0d0 - xd3) + c1 * xd3
}

end subroutine interpolate_var
!=============================================================================
subroutine get_vec(igrid,x,tloc,var,ibeg,iend)

include 'amrvacdef.f'

integer,intent(in)                                   :: igrid, ibeg, iend
double precision,dimension(^NC), intent(in)          :: x
double precision, intent(in)                         :: tloc
double precision,dimension(iend-ibeg+1), intent(out) :: var
double precision,dimension(iend-ibeg+1)              :: e1, e2
integer                                              :: ivar, iloc
double precision                                     :: td
!-----------------------------------------------------------------------------

if (.not.time_advance) then

   do ivar=ibeg,iend
      iloc = ivar-ibeg+1

      call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ivar),px(igrid)%x(ixG^T,1:ndim),x,var(iloc))

   end do
   
else
   
   td = (tloc/(UNIT_LENGTH/UNIT_VELOCITY) - t) / dt
   
   do ivar=ibeg,iend
      iloc = ivar-ibeg+1

      call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars_old(igrid)%w(ixG^T,ivar),px(igrid)%x(ixG^T,1:ndim),x,e1(iloc))
      call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ivar),px(igrid)%x(ixG^T,1:ndim),x,e2(iloc))
      
      var(iloc) = e1(iloc) * (1.0d0 - td) + e2(iloc) * td

   end do

end if !.not.time_advance

end subroutine get_vec
!=============================================================================
end module mod_gridvars
!=============================================================================
}
