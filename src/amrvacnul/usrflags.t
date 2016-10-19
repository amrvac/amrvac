subroutine special_process_usr
use mod_global_parameters
use mod_input_output

{#IFDEF MAGNETOFRICTION
if(itmaxmf>0) then
   time_in=MPI_WTIME()
   call magnetofriction
   if(mype==0) write(*,*) 'Magnetofriction phase took : ',MPI_WTIME()-time_in,' sec'
endif
}
if(nwtf>0 .and. neqpartf>0) then
  call write_snapshot_tf
end if

end subroutine special_process_usr
!=============================================================================
!> regenerate w and eqpar arrays to output into *tf.dat
subroutine transformw_usr(w,wtf,eqpar_tf,ixI^L,ixO^L)
use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: w(ixI^S,1:nw)
double precision, intent(out):: wtf(ixI^S,1:nwtf)
double precision, intent(out):: eqpar_tf(neqpartf)

double precision :: gamma_usr
integer :: iwup_,diw,e_usr
!-----------------------------------------------------------------------------
{^IFMHD
gamma_usr=5.d0/3.d0
e_usr=m^NC_+1
wtf(ixO^S,1:m^NC_)=w(ixO^S,1:m^NC_)
wtf(ixO^S,e_usr)=w(ixO^S,rho_)*1.d0/(gamma_usr-one)+&
      half*((^C&w(ixO^S,m^C_)**2.0d0+)/w(ixO^S,rho_)+(^C&w(ixO^S,b^C_)**2.0d0+))
wtf(ixO^S,e_usr+1:e_usr+^NC)=w(ixO^S,b1_:b^NC_)
iwup_=b^NC_
if(iwup_<nw .and. nw<nwtf) then
  diw=nw-iwup_
  wtf(ixO^S,e_usr+^NC+1:e_usr+^NC+diw)=w(ixO^S,iwup_+1:iwup_+diw)
endif
eqpar_tf(1:4)=eqpar(1:4)
}
                                                                                                                                                          
end subroutine transformw_usr
!=============================================================================
subroutine fixp_usr(ixI^L,ixO^L,w,x)
use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(inout)    :: w(ixI^S,1:nw)
double precision, intent(in)       :: x(ixI^S,1:ndim)
!----------------------------------------------------------------------------


end subroutine fixp_usr
!=============================================================================
subroutine flag_grid_usr(qt,ixG^L,ixO^L,w,x,flag)

use mod_global_parameters

integer, intent(in)             :: ixG^L, ixO^L
integer, intent(inout)          :: flag
double precision, intent(in)    :: qt
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision, intent(in)    :: x(ixG^S,1:ndim)

! flag=-1 : Treat all cells active, omit deactivation (onentry, default)
! flag=0  : Treat as normal domain
! flag=1  : Treat as passive, but reduce by safety belt
! flag=2  : Always treat as passive

!-----------------------------------------------------------------------------
      
end subroutine flag_grid_usr
!=============================================================================
