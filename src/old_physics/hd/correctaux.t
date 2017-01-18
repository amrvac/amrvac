!##############################################################################
! module hd/correctaux.t version july 2009 
! merged with mhdiso, September 2012 by Oliver Porth
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

{do ix^DB= ixO^LIM^DB\}
    ! point with local failure identified by patchierror
    ! patchierror=-1 then from smallvalues call
    ! patchierror=1  then from getpthermal or primitive call
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
           if (useprimitive.and.subname/='primitive') then
              patchw(kxO^S)=(patchierror(kxO^S)/=0)
              call primitiven(ixI^L,kxO^L,w,patchw)
           end if
           ! faulty cells are corrected by averaging here
           ! only average those which were ok and replace faulty cells
           do iw = 1,nw
              w(ix^D,iw)=sum(w(kxO^S,iw),patchierror(kxO^S)==0)&
                           /count(patchierror(kxO^S)==0)
           end do
           if (useprimitive.and.subname/='primitive') then
              ! in addition to those switched to primitive variables
              ! above, also switch the corrected variables
              patchw(ix^D)=.false.
              call conserven(ixI^L,kxO^L,w,patchw)
           end if
        else
           ! no cells without error were found in cube of size nflatgetaux
           ! --> point of no recovery
           print*,'Getaux error:',patchierror(ix^D),'ix^D=',ix^D
           !print*,'New ','rho=',w(ix^D,rho_),'m=', &
           !        {^C&w(ix^D,m^C_)},'e=',w(ix^D,e_)
           !print*,'position  ', px(saveigrid)%x(ix^D, 1:ndim)
           print*,'Called from: ',subname
           if (patchierror(ix^D)<0) then
               call mpistop("-correctaux from smallvalues-----")
           else
               call mpistop("-correctaux from primitive or getpthermal--")
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
character(len=*), intent(in)    :: subname
!.. local ..
integer                         :: posvec(ndim)
integer, dimension(ixG^T)       :: patchierror
double precision                :: pth(ixG^T)
!-----------------------------------------------------------------------------

{#IFDEF ENERGY
pth(ixO^S)=(eqpar(gamma_)-one)*(w(ixO^S,e_)- &
       half*(({^C&w(ixO^S,m^C_)**2+})/w(ixO^S,rho_)))
if(any(w(ixO^S,rho_) < minrho).or.any(w(ixO^S,e_) < smalle)&
  .or.any(pth(ixO^S) < minp))then
  if(strictsmall)then
     write(*,*)'SMALLVALUES under strictsmall problem From:  ', &
             subname,' iteration ', it,' time ',t
     if(any(w(ixO^S,rho_) < minrho)) then
       posvec(1:ndim)=minloc(w(ixO^S,rho_))
       ^D&posvec(^D)=posvec(^D)+ixOmin^D-1;
       write(*,*)'minimum rho =', minval(w(ixO^S,rho_)),' with limit=',minrho,&
       ' at x=', x(^D&posvec(^D),1:ndim),' array index=',posvec,' where E_k=',&
       half*(^C&w(^D&posvec(^D),m^C_)**2+)/w(^D&posvec(^D),rho_),&
       ' E_total=',w(^D&posvec(^D),e_),&
       ' w(1:nwflux)=',w(^D&posvec(^D),1:nwflux)
     endif
     if(any(w(ixO^S,e_) < smalle)) then
       posvec(1:ndim)=minloc(w(ixO^S,e_))
       ^D&posvec(^D)=posvec(^D)+ixOmin^D-1;
       write(*,*)'minimum e =', minval(w(ixO^S,e_)),' with limit=',smalle,&
       ' at x=', x(^D&posvec(^D),1:ndim),' array index=',posvec,' where E_k=',&
       half*(^C&w(^D&posvec(^D),m^C_)**2+)/w(^D&posvec(^D),rho_),&
       ' E_total=',w(^D&posvec(^D),e_),&
       ' w(1:nwflux)=',w(^D&posvec(^D),1:nwflux)
     endif
     if(any(pth(ixO^S) < minp)) then
       posvec(1:ndim)=minloc(pth(ixO^S))
       ^D&posvec(^D)=posvec(^D)+ixOmin^D-1;
       write(*,*)'minimum pressure = ', minval(pth(ixO^S)),' with limit=',minp,&
       ' at x=',x(^D&posvec(^D),1:ndim),' array index=',posvec,' where rho=',&
       w({^D&posvec(^D)},rho_),', velocity v=',&
       ^C&w({^D&posvec(^D)},m^C_)/w({^D&posvec(^D)},rho_),&
      ' w(1:nwflux)=',w(^D&posvec(^D),1:nwflux)
     endif
     call mpistop("Smallvalues with strictsmall=T failed")
  else
     if(strictgetaux)then
       where(w(ixO^S,rho_) < minrho .or. w(ixO^S,e_) < smalle&
             .or. pth(ixO^S)<minp)
         w(ixO^S,rho_)  = 2.0*(1.0d0 + 10.0d0 * minrho)*minrho
         w(ixO^S,e_)    = 2.0*(1.0d0 + 10.0d0 * minp)*smalle
         {^C&w(ixO^S,m^C_) =zero;}
       end where
     else
       where(w(ixO^S,rho_) < minrho .or. w(ixO^S,e_) < smalle&
             .or. pth(ixO^S)<minp)
         patchierror(ixO^S)=-1
       elsewhere
         patchierror(ixO^S)=0
       end where
       call correctaux(ixI^L,ixO^L,w,x,patchierror,subname)
     end if
  end if ! strict
end if
}
{#IFNDEF ENERGY
if(any(w(ixO^S,rho_) < minrho))then
  if(strictsmall)then
     write(*,*)'SMALLVALUES under strictsmall problem From:  ', &
             subname,' iteration ', it,' time ',t
     posvec(1:ndim)=minloc(w(ixO^S,rho_))
     ^D&posvec(^D)=posvec(^D)+ixOmin^D-1;
     write(*,*)'minimum rho =', minval(w(ixO^S,rho_)),' with limit=',minrho,&
     ' at x=', x(^D&posvec(^D),1:ndim),' array index=',posvec,' where E_k=',&
     half*(^C&w(^D&posvec(^D),m^C_)**2+)/w(^D&posvec(^D),rho_),&
     ' w(1:nwflux)=',w(^D&posvec(^D),1:nwflux)
     call mpistop("Smallvalues with strictsmall=T failed")
  else
     if(strictgetaux)then
       where(w(ixO^S,rho_) < minrho )
         w(ixO^S,rho_)  = 2.0*(1.0d0 + 10.0d0 * minrho)*minrho
         {^C&w(ixO^S,m^C_) =zero;}
       end where
     else
       where(w(ixO^S,rho_) < minrho)
         patchierror(ixO^S)=-1
       elsewhere
         patchierror(ixO^S)=0
       end where
       call correctaux(ixI^L,ixO^L,w,x,patchierror,subname)
     end if
  end if ! strict
end if
}
end subroutine smallvalues
!=============================================================================
! end module hd/correctaux.t
!##############################################################################
