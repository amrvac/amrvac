!#############################################################################
! module amrvacphys.correctauxsrmhdiso .t version april 2013 (F.Casse)
!=============================================================================
subroutine smallvalues(w,x,ixI^L,ixO^L,patchw,subname)

! difference with srhd case only in padding array patchw

include 'amrvacdef.f'

integer, intent(in)          :: ixI^L,ixO^L
double precision             :: w(ixI^S,1:nw)
double precision, intent(in) :: x(ixI^S,1:ndim)
logical, intent(in)          :: patchw(ixG^T)
character(len=*), intent(in) :: subname

integer                      :: patchierror(ixG^T)
!-----------------------------------------------------------------------------
if(any((w(ixO^S,d_) < minrho)))then
  if(strictsmall)then
     print*,'smallvalues :: ', 'for tau = ',minval(w(ixO^S,tau_))&
        ,' With limit=',smalltau ,' For D = ', minval(w(ixO^S,d_))&
        ,' With limit=', smallrho
     print*,'From::  ', subname,' iteration ', it
     call mpistop("Smallvalues with strictsmall=T failed")
  else
     if(strictgetaux)then
       ! (optional) artificial replacement of small D and tau values
       ! above corresponding smallrho/smallp settings,
       ! together with nullifying momenta
       where((w(ixO^S,d_) < minrho) &
           .and. .not.patchw(ixO^S))
         w(ixO^S,d_)   = 2.0d0*(1.0d0 + 10.0d0 * minrho)*minrho
         w(ixO^S,tau_) = 2.0d0*(1.0d0 + 10.0d0 * minp)*smalltau&
                      + ({^C&w(ixO^S,b^C_)**2.0d0+})/2.0d0
         {^C&w(ixO^S,s^C_) =zero;}
         w(ixO^S,lfac_)=one
       end where
     else
       where((w(ixO^S,d_) < minrho)&
           .and. .not.patchw(ixO^S))
          patchierror(ixO^S)=-1
       else where
          patchierror(ixO^S)=0
       end where
       call correctaux_usr(ixI^L,ixO^L,w,x,patchierror,subname)
       call correctaux(ixI^L,ixO^L,w,x,patchierror,subname)
     end if
  end if ! strict
end if
{#IFDEF EPSINF 
if(any(w(ixO^S,Drho1_)<eqpar(rho1e_) &
     .or. w(ixO^S,Dn_)<eqpar(rho1e_))) then 
   where(w(ixO^S,Drho1_)<eqpar(rho1e_) .and. .not.patchw(ixO^S))
      w(ixO^S,Drho1_)=eqpar(rho1e_)
      w(ixO^S,Drho0_)=w(ixO^S,Drho1_)*eqpar(rho0floor_)
      w(ixO^S,Depsinf_) = eqpar(epsfloor_)*w(ixO^S,Drho1_)**(2.0d0/3.0d0)
   end where
   where(w(ixO^S,Dn_)<eqpar(rho1e_) .and. .not.patchw(ixO^S))
      w(ixO^S,Dn_)=eqpar(rho1e_)
      w(ixO^S,Dn0_)=w(ixO^S,Dn_)*eqpar(rho0floor_)
   end where
end if\}
end subroutine smallvalues
!=============================================================================
subroutine getaux(clipping,w,x,ixI^L,ixO^L,subname)

! Calculate the auxiliary variables lfac and xi within ixO^L
! clipping is not used (yet) 

include 'amrvacdef.f'

logical, intent(in)            :: clipping
integer, intent(in)            :: ixI^L, ixO^L
double precision, intent(inout):: w(ixI^S,1:nw)
double precision, intent(in)   :: x(ixI^S,1:ndim)
character(len=*), intent(in)   :: subname

double precision:: divb(ixG^T), lorleft, b2left
integer::          err,ix^D
double precision :: dold,tauold,{^C&sold^C_},{^C&bold^C_}
integer          :: patchierror(ixG^T)
logical          :: patchw(ixG^T)
!-----------------------------------------------------------------------------
patchw(ixO^S)=.false.
if (fixsmall) call smallvalues(w,x,ixI^L,ixO^L,patchw,subname)

! we compute auxiliaries lfac,xi from D,S,tau,B
! put the lfac and xi in the auxiliary fields lfac_ and xi_
{do ix^D= ixO^LIM^D\}
    dold=w(ix^D,d_)
    tauold=w(ix^D,tau_)
    { ^C&sold^C_=w(ix^D,s^C_);}
    { ^C&bold^C_=w(ix^D,b^C_);}

    call con2prim(w(ix^D,lfac_),w(ix^D,xi_), &
             w(ix^D,d_),{^C&w(ix^D,s^C_)},w(ix^D,tau_),{^C&w(ix^D,b^C_)},err)
    patchierror(ix^D)=err

    if (err/=0.and.strictgetaux) then
       print*,'Getaux error:',err,'ix^D=',ix^D
       print*,'input d=',dold,'s=',{^C&sold^C_},'tau=',tauold,'b=',{^C&bold^C_}
       print*,'New ','d=',w(ix^D,d_),'s=',{^C&w(ix^D,s^C_)},'tau=',w(ix^D,tau_)
       print*,'B=',{^C&w(ix^D,b^C_)}
       print*,'position  ', x(ix^D, 1:ndim)
       print*,'index ', ix^D
       print*,'Called from: ',subname
       call mpistop("problem in getaux: retry with strictgetaux=F")
    end if

{enddo^D&\}
! first try user-defined method:
if(.not.strictgetaux.and.any(patchierror(ixO^S)/=0)) &
     call correctaux_usr(ixI^L,ixO^L,w,x,patchierror,subname)
! then, if any patchierror was not reset to 0, use default:
if(.not.strictgetaux.and.any(patchierror(ixO^S)/=0)) &
     call correctaux(ixI^L,ixO^L,w,x,patchierror,subname)

!if(any(patchierror(ixO^S)/=0)) call mpistop('getaux final error')


end subroutine getaux
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!=============================================================================
subroutine correctaux(ixI^L,ixO^L,w,x,patchierror,subname)

include 'amrvacdef.f'

integer, intent(in)            :: ixI^L, ixO^L
integer, intent(in)            :: patchierror(ixG^T)
character(len=*), intent(in)   :: subname
double precision, intent(inout):: w(ixI^S,1:nw)
double precision, intent(in)   :: x(ixI^S,1:ndim)

integer        :: iw, kxO^L, ix^D, i
logical        :: patchw(ixG^T)
!-----------------------------------------------------------------------------

{do ix^D= ixO^LIM^D\}
    if (patchierror(ix^D)/=0) then
        do i=1,nflatgetaux
           {kxOmin^D= max(ix^D-i,ixOmin^D);
            kxOmax^D= min(ix^D+i,ixOmax^D);\}
           if (any(patchierror(kxO^S)==0)) exit
        end do
        if (any(patchierror(kxO^S)==0))then
            ! in contrast to srhd case: always switch to primitive
            ! as error can become large when this is not done
            patchw(kxO^S)=(patchierror(kxO^S)/=0)
            call primitiven(ixI^L,kxO^L,w,patchw)
            do iw = 1,nw
               ! in contrast to srhd: do not alter magnetic field
               if (iw/=b^C_|.and.) then
                   w(ix^D,iw)=sum(w(kxO^S,iw),patchierror(kxO^S)==0)&
                           /count(patchierror(kxO^S)==0)
               end if  
            end do
            patchw(ix^D)=.false.
            call conserven(ixI^L,kxO^L,w,patchw)
        else
           ! no cells without error were found in cube of size nflatgetaux
           ! --> point of no recovery
           print*,'Getaux error:',patchierror(ix^D),'ix^D=',ix^D
           print*,'New ','d=',w(ix^D,d_),'s=',{^C&w(ix^D,s^C_)},'tau=',w(ix^D,tau_)
           print*,'B=',{^C&w(ix^D,b^C_)},w(ix^D,lfac_),w(ix^D,xi_)
           print*,'position  ', x(ix^D, 1:ndim)
           print*,'index ', ix^D
           print*,'Called from: ',subname
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
! end module amrvacphys.correctauxsrmhdeos.t
!#############################################################################
