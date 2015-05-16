!##############################################################################
! module mhd/correctaux.t version July 2009 
! merged with mhdiso, September 2012 by Oliver Porth
!=============================================================================
subroutine correctaux(ixI^L,ixO^L,w,x,patchierror,subname)

include 'amrvacdef.f'

integer, intent(in)         :: ixI^L, ixO^L
integer, intent(in)         :: patchierror(ixG^T)
character(len=*), intent(in)   :: subname
double precision, intent(inout):: w(ixI^S,1:nw)
double precision, intent(in)      :: x(ixI^S,1:ndim)
integer        :: iw, kxO^L, ix^D, i
logical        :: patchw(ixI^S)
!-----------------------------------------------------------------------------

{do ix^DB= ixO^LIM^DB\}
    ! point with local failure identified by patchierror
    ! patchierror=-1 then from smallvalues call
    ! patchierror=1  then from getpthermal/total or primitive call
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
           if (subname/='primitive') then
              patchw(kxO^S)=(patchierror(kxO^S)/=0)
              call primitiven(ixI^L,kxO^L,w,patchw)
           endif
           ! faulty cells are corrected by averaging here
           ! only average those which were ok and replace faulty cells
           do iw = 1,nw
              if (iw/=b^C_|.and.) then
                   w(ix^D,iw)=sum(w(kxO^S,iw),patchierror(kxO^S)==0)&
                              /count(patchierror(kxO^S)==0)
              end if
           end do
           if (subname/='primitive') then
              ! in addition to those switched to primitive variables
              ! above, also switch the corrected variables
              patchw(ix^D)=.false.
              call conserven(ixI^L,kxO^L,w,patchw)
           endif
        else
           ! no cells without error were found in cube of size nflatgetaux
           ! --> point of no recovery
           print*,'Getaux error:',patchierror(ix^D),'ix^D=',ix^D
           print*,'New ','rho=',w(ix^D,rho_),'m=', &
                   {^C&w(ix^D,m^C_)}
           print*,'B=',{^C&w(ix^D,b^C_)}
{#IFDEF ENERGY
           print*,'pressure ', (eqpar(gamma_)-one)*(w(ix^D,e_)- &
                           half*(({^C&w(ix^D,m^C_)**2+})/w(ix^D,rho_)&
                                 +{^C&w(ix^D,b^C_)**2+}))
           print*,'e=',w(ix^D,e_)
}
           print*,'position ', x(ix^D, 1:ndim),' time ',t,it
           print*,'Called from: ',subname
           if (patchierror(ix^D)<0) then
               call mpistop("-correctaux from smallvalues-----")
           else
               call mpistop("-correctaux from primitive-getpthermal-total-")
           end if
        end if
    end if
{enddo^D&\}

end subroutine correctaux
!=============================================================================
subroutine smallvalues(w,x,ixI^L,ixO^L,subname)

include 'amrvacdef.f'

integer, intent(in)             :: ixI^L,ixO^L
double precision, intent(inout) :: w(ixI^S,1:nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
character(len=*), intent(in)    ::subname
!.. local ..
integer                         :: posvec(ndim)
integer, dimension(ixG^T)       :: patchierror
double precision                :: pth(ixG^T), Te(ixG^T)
!-----------------------------------------------------------------------------

{#IFDEF ENERGY
pth(ixO^S)=(eqpar(gamma_)-one)*(w(ixO^S,e_)- &
       half*(({^C&w(ixO^S,m^C_)**2+})/w(ixO^S,rho_)&
       +{ ^C&w(ixO^S,b^C_)**2+}))
if(smallT>0.d0) Te(ixO^S)=pth(ixO^S)/w(ixO^S,rho_)
if(strictsmall) then
  if(smallT>0.d0 .and. any(Te(ixO^S) <=smallT)) then
    print *,'SMALLVALUES of temperature under strictsmall problem From:  ', &
            subname,' iteration=', it,' time=',t
    posvec(1:ndim)=minloc(Te(ixO^S))
    ^D&posvec(^D)=posvec(^D)+ixOmin^D-1;
    write(*,*)'minimum temperature= ', minval(Te(ixO^S)),' with limit=',smallT,&
    ' at x=', x(^D&posvec(^D),1:ndim),' array index=',posvec,' where rho=',&
    w({^D&posvec(^D)},rho_),', velocity v=',&
    ^C&w({^D&posvec(^D)},m^C_)/w({^D&posvec(^D)},rho_),&
    ', and magnetic field B=',^C&w({^D&posvec(^D)},b^C_),&
    ' w(1:nwflux)=',w(^D&posvec(^D),1:nwflux)
    call mpistop("Smallvalues of temperature with strictsmall=T failed")
  endif
  if(any(pth(ixO^S) <=minp)) then
    print *,'SMALLVALUES of pressure under strictsmall problem From:  ', &
            subname,' iteration=', it,' time=',t
    posvec(1:ndim)=minloc(pth(ixO^S))
    ^D&posvec(^D)=posvec(^D)+ixOmin^D-1;
    write(*,*)'minimum pressure = ', minval(pth(ixO^S)),' with limit=',minp,&
    ' at x=', x(^D&posvec(^D),1:ndim),' array index=',posvec,' where rho=',&
    w({^D&posvec(^D)},rho_),', velocity v=',&
    ^C&w({^D&posvec(^D)},m^C_)/w({^D&posvec(^D)},rho_),&
    ', and magnetic field B=',^C&w({^D&posvec(^D)},b^C_),&
    ' w(1:nwflux)=',w(^D&posvec(^D),1:nwflux)
    call mpistop("Smallvalues of pressure with strictsmall=T failed")
  endif
  if(any(w(ixO^S,e_) <=smalle)) then
    print *,'SMALLVALUES of energy under strictsmall problem From:  ', &
            subname,' iteration=', it,' time=',t
    posvec(1:ndim)=minloc(w(ixO^S,e_))
    ^D&posvec(^D)=posvec(^D)+ixOmin^D-1;
    write(*,*)'minimum e =', minval(w(ixO^S,e_)),' with limit=',smalle,&
    ' at x=', x(^D&posvec(^D),1:ndim),' array index=',posvec,' where E_k=',&
    half*(^C&w(^D&posvec(^D),m^C_)**2+)/w(^D&posvec(^D),rho_),&
    ' E_total=',w(^D&posvec(^D),e_),&
    ' w(1:nwflux)=',w(^D&posvec(^D),1:nwflux)
    call mpistop("Smallvalues of energy with strictsmall=T failed")
  endif
  if(any(w(ixO^S,rho_) <=minrho)) then
    print *,'SMALLVALUES of density under strictsmall problem From:  ', &
            subname,' iteration=', it,' time=',t
    posvec(1:ndim)=minloc(w(ixO^S,rho_))
    ^D&posvec(^D)=posvec(^D)+ixOmin^D-1;
    write(*,*)'minimum rho =', minval(w(ixO^S,rho_)),' with limit=',minrho,&
    ' at x=', x(^D&posvec(^D),1:ndim),' array index=',posvec,' where E_k=',&
    half*(^C&w(^D&posvec(^D),m^C_)**2+)/w(^D&posvec(^D),rho_),&
    ' E_total=',w(^D&posvec(^D),e_),&
    ' w(1:nwflux)=',w(^D&posvec(^D),1:nwflux)
    call mpistop("Smallvalues of density with strictsmall=T failed")
  endif
else
  if(strictgetaux)then
    where(w(ixO^S,rho_) < minrho)
      w(ixO^S,rho_)=minrho
      {^C&w(ixO^S,m^C_) =zero;}
    end where
    where(pth(ixO^S) < minp)
      w(ixO^S,e_)=minp/(eqpar(gamma_)-one)+&
       (({^C&w(ixO^S,m^C_)**2+})/w(ixO^S,rho_)+{^C&w(ixO^S,b^C_)**2+})*half
    end where
    where(Te(ixO^S) < smallT)
      w(ixO^S,e_)=smallT*w(ixO^S,rho_)/(eqpar(gamma_)-one)+&
       (({^C&w(ixO^S,m^C_)**2+})/w(ixO^S,rho_)+{^C&w(ixO^S,b^C_)**2+})*half
    end where
  else
    where(w(ixO^S,rho_) < minrho .or. w(ixO^S,e_) < smalle&
          .or. pth(ixO^S) < minp .or. Te(ixO^S) < smallT)
      patchierror(ixO^S)=-1
    elsewhere
      patchierror(ixO^S)=0
    end where
    call correctaux(ixI^L,ixO^L,w,x,patchierror,subname)
  end if
end if
}

{#IFDEF ISO
if(any(w(ixO^S,rho_) < minrho)) then
  if(strictsmall)then
     write(*,*)'SMALLVALUES of density under strictsmall problem From:  ', &
             subname,' iteration ', it,' time ',t
     posvec(1:ndim)=minloc(w(ixO^S,rho_))
     ^D&posvec(^D)=posvec(^D)+ixOmin^D-1;
     write(*,*)'minimum rho =', minval(w(ixO^S,rho_)),' with limit=',minrho,&
     ' at x=', x(^D&posvec(^D),1:ndim),' array index=',posvec,' where E_k=',&
     half*(^C&w(^D&posvec(^D),m^C_)**2+)/w(^D&posvec(^D),rho_),&
     ' w(1:nwflux)=',w(^D&posvec(^D),1:nwflux)
     call mpistop("Smallvalues of density with strictsmall=T failed")
  else
     if(strictgetaux)then
       where(w(ixO^S,rho_) < minrho)
         w(ixO^S,rho_)  = 2.0*(1.0d0 + 10.0d0 * minrho)*minrho
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
! end module mhd/correctaux.t
!##############################################################################
