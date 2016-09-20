!##############################################################################
! module amrvacphys.correctauxsrhd.t version july 2009 
!=============================================================================
subroutine correctaux(ixI^L,ixO^L,w,x,patchierror,subname)

use mod_global_parameters

integer, intent(in)         :: ixI^L, ixO^L
integer, intent(in)         :: patchierror(ixG^T)
character(len=*), intent(in)   :: subname
double precision, intent(inout):: w(ixI^S,1:nw)
double precision, intent(in)      :: x(ixI^S,1:ndim)

integer        :: iw, kxO^L, ix^D, i
logical        :: patchw(ixG^T)
!-----------------------------------------------------------------------------

{do ix^D= ixO^LIM^D\}
    ! point with local failure identified by patchierror
    ! patchierror=-1 then from smallvalues call
    ! patchierror=1  then from getaux      call
    if (patchierror(ix^D)/=0) then
        ! verify in cube with border width nflatgetaux the presence
        ! of cells where all went ok
        do i=1,nflatgetaux
           {kxOmin^D= max(ix^D-i,ixOmin^D);
            kxOmax^D= min(ix^D+i,ixOmax^D);\}
           ! in case cells are fine within smaller cube than 
           ! the userset nflatgetaux: use that smaller cube
           if (any(patchierror(kxO^S)==0)) exit
        end do
        if (any(patchierror(kxO^S)==0))then
           ! within surrounding cube, cells without problem were found
           if (useprimitive) then
              patchw(kxO^S)=(patchierror(kxO^S)/=0)
              call primitiven(ixI^L,kxO^L,w,patchw)
           end if
           ! faulty cells are corrected by averaging here
           ! only average those which were ok and replace faulty cells
           do iw = 1,nw
              w(ix^D,iw)=sum(w(kxO^S,iw),patchierror(kxO^S)==0)&
                           /count(patchierror(kxO^S)==0)
           end do
           if (useprimitive) then
              ! in addition to those switched to primitive variables
              ! above, also switch the corrected variables
              patchw(ix^D)=.false.
              call conserven(ixI^L,kxO^L,w,patchw)
           end if
        else
           ! no cells without error were found in cube of size nflatgetaux
           ! --> point of no recovery
           write(*,*)'Getaux error:',patchierror(ix^D),'ix^D=',ix^D
           !write(*,*)'New ','d=',w(ix^D,d_),'s=', &
           !        {^C&w(ix^D,s^C_)},'tau=',w(ix^D,tau_)
           !write(*,*)'position  ', px(saveigrid)%x(ix^D, 1:ndim)
           write(*,*)'Called from: ',subname
           if (patchierror(ix^D)<0) then
               call mpistop("---------correctaux from smallvalues-----")
           else
               call mpistop("---------correctaux from getaux----------")
           end if
        end if
    end if
{enddo^D&\}

end subroutine correctaux
!=============================================================================
subroutine smallvalues(w,x,ixI^L,ixO^L,subname)

use mod_global_parameters

integer, intent(in)             :: ixI^L,ixO^L
double precision, intent(inout) :: w(ixI^S,1:nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
character(len=*), intent(in)    ::subname

!!integer                         :: posvec(ndim)
integer, dimension(ixG^T)       :: patchierror
!-----------------------------------------------------------------------------

if(any(w(ixO^S,d_) < minrho) .or. any(w(ixO^S,tau_) < smalltau))then
  if(strictsmall)then
     write(*,*)'smallvalues :: for tau = ',minval(w(ixO^S,tau_)), &
         'With limit=',smalltau,' For D =', minval(w(ixO^S,d_)),&
         ' With limit=',smallrho
     write(*,*)'From::  ', subname,' iteration ', it
     call mpistop("Smallvalues with strictsmall=T failed")
  else
     if(strictgetaux)then
       where(w(ixO^S,d_) < minrho .or. w(ixO^S,tau_) < smalltau)
         w(ixO^S,d_)  = 2.0*(1.0d0 + 10.0d0 * minrho)*minrho
         w(ixO^S,tau_)= 2.0*(1.0d0 + 10.0d0 * minp)*smalltau
         {^C&w(ixO^S,s^C_) =zero;}
         w(ixO^S,lfac_)=one
       end where
     else
       where(w(ixO^S,d_) < minrho .or. w(ixO^S,tau_) < smalltau)
         patchierror(ixO^S)=-1
       else where
         patchierror(ixO^S)=0
       end where
       call correctaux(ixI^L,ixO^L,w,x,patchierror,subname)
     end if
  end if ! strict
end if

{#IFDEF TRACER
if(any(w(ixO^S,Dtr1_) < minrho) .or. &
   any(w(ixO^S,Dtr1_) > 10.d10 ))then
	where(w(ixO^S,Dtr1_) < minrho .or. &
              w(ixO^S,Dtr1_) > 10.d10) 
		w(ixO^S,Dtr1_) = 0.d0
	end where
end if
if(any(w(ixO^S,Dtr1_) .NE. w(ixO^S,Dtr1_)))then
	where(w(ixO^S,Dtr1_) .NE. w(ixO^S,Dtr1_)) 
		w(ixO^S,Dtr1_) = 0.d0
	end where
end if\}

{#IFDEF EPSINF 
! fix to floor:
if(any(w(ixO^S,Dne_)<eqpar(rho1e_)) &
     .or. any(w(ixO^S,Dne0_)<eqpar(rho1e_)*eqpar(rho0floor_))) then 
   where(w(ixO^S,Dne_)<eqpar(rho1e_).or. w(ixO^S,Dne0_)<eqpar(rho1e_)*eqpar(rho0floor_))
      w(ixO^S,Dne_)=eqpar(rho1e_)
      w(ixO^S,Dne0_)=w(ixO^S,Dne_)*eqpar(rho0floor_)
      w(ixO^S,Depsinf_) = eqpar(epsfloor_)*w(ixO^S,Dne_)**(2.0d0/3.0d0)
   end where
end if\}

end subroutine smallvalues
!=============================================================================
subroutine getaux(clipping,w,x,ixI^L,ixO^L,subname)

! Calculate auxilary variables ixO^L from non-auxiliary entries in w
! clipping can be set to .true. to e.g. correct unphysical pressures,
! densities, v>c,  etc.

use mod_global_parameters

logical, intent(in)             :: clipping
integer, intent(in)             :: ixI^L, ixO^L
double precision, intent(inout) :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
character(len=*), intent(in)    :: subname

integer          :: err,ix^D
double precision :: dold,tauold,{^C&sold^C_},pold,lfacold
integer          :: patchierror(ixG^T)
!-----------------------------------------------------------------------------

! artificial replacement of small D and tau values
! with corresponding smallrho/smallp settings,
! together with nullifying momenta
call smallvalues(w,x,ixI^L,ixO^L,subname)

! we compute auxiliaries p, lfac from D,S,tau
! put the p and lfac in the auxiliary fields lfac_ and p_
! on entry: p_ field may contain previous value for iteration
! however, for filling ghost cells, this does not exist, so we no longer
! use it

{do ix^D= ixO^LIM^D\}
     dold=w(ix^D,d_)
     tauold=w(ix^D,tau_)
     pold=w(ix^D,p_)
     lfacold=w(ix^D,lfac_)
    { ^C&sold^C_=w(ix^D,s^C_);}

    call con2prim(w(ix^D,p_),w(ix^D,lfac_), &
            w(ix^D,d_),{^C&w(ix^D,s^C_)},w(ix^D,tau_),err)
    patchierror(ix^D)=err
    if (err/=0.and.strictgetaux) then
       write(*,*)'Getaux error:',err,'ix^D=',ix^D,' global timestep it=',it
       write(*,*)'Postion: ',x(ix^D,1:ndim)
       write(*,*)'start value for p=',pold
       write(*,*)'start value for lfac=',lfacold
       write(*,*)'end value for p=',w(ix^D,p_)
       write(*,*)'end value for lfac=',w(ix^D,lfac_)
       write(*,*)'exit  d=',w(ix^D,d_),'s=',{^C&w(ix^D,s^C_)},'tau=',w(ix^D,tau_)
       write(*,*)'input d=',dold,'s=',{^C&sold^C_},'tau=',tauold
       write(*,*)'Called from: ',subname
       call mpistop("problem in getaux")
   endif
{enddo^D&\}

if(.not.strictgetaux.and.any(patchierror(ixO^S)/=0)) &
     call correctaux(ixI^L,ixO^L,w,x,patchierror,subname)


end subroutine getaux
!=============================================================================
! end module amrvacphys.correctauxsrhd.t
!##############################################################################
