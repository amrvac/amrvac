subroutine correctaux_usr(ixI^L,ixO^L,w,x,patchierror,subname)

include 'amrvacdef.f'

integer, intent(in)            :: ixI^L, ixO^L
integer, intent(inout)         :: patchierror(ixG^T)
character(len=*), intent(in)   :: subname
double precision, intent(inout):: w(ixI^S,1:nw)
double precision, intent(in)   :: x(ixI^S,1:ndim)

! correct solution from analytic case

end subroutine correctaux_usr
!==========================================================================================
