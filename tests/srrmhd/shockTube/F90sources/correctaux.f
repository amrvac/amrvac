!#############################################################################
! module correctaux.t version September 2012
! SRMHD physics module
!=============================================================================
subroutine smallvalues(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,patchw,subname)

! difference with srhd case only in padding array patchw

include 'amrvacdef.f'

integer, intent(in)          :: ixImin1,ixImax1,ixOmin1,ixOmax1
double precision             :: w(ixImin1:ixImax1,1:nw)
double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
logical, intent(in)          :: patchw(ixGlo1:ixGhi1)
character(len=*), intent(in) :: subname

integer                      :: patchierror(ixGlo1:ixGhi1), ixmin(ndim)
!-----------------------------------------------------------------------------
if(any((w(ixOmin1:ixOmax1,d_) < minrho .or. w(ixOmin1:ixOmax1,tau_) &
   < smalltau))) then
  if(strictsmall)then
      ixmin=minloc(w(ixOmin1:ixOmax1,d_))
     write(*,*),'smallvalues :: ', 'for tau = ',minval(w(ixOmin1:ixOmax1,&
        tau_)),' With limit=',smalltau ,' For D = ', minval(w(ixOmin1:ixOmax1,&
        d_)),' With limit=', minrho
     write(*,*),'From::  ', subname,' iteration ', it
        write(*,*),'position  ', x(ixmin(1), 1:ndim)
     call mpistop("Smallvalues with strictsmall=T failed")
  else
     if(strictgetaux)then
       ! (optional) artificial replacement of small D and tau values
       ! above corresponding smallrho/smallp settings,
       ! together with nullifying momenta
       where((w(ixOmin1:ixOmax1,d_) < minrho .or. w(ixOmin1:ixOmax1,tau_) &
          < smalltau).and. .not.patchw(ixOmin1:ixOmax1))
         w(ixOmin1:ixOmax1,d_)   = 2.0d0*(1.0d0 + 10.0d0 * minrho)*minrho
         w(ixOmin1:ixOmax1,tau_) = 2.0d0*(1.0d0 + 10.0d0 * minp)*smalltau&
            + (w(ixOmin1:ixOmax1,b1_)**2.0d0+w(ixOmin1:ixOmax1,b2_)**2.0d0&
            +w(ixOmin1:ixOmax1,b3_)**2.0d0)/2.0d0
         w(ixOmin1:ixOmax1,s1_) =zero;w(ixOmin1:ixOmax1,s2_) =zero
         w(ixOmin1:ixOmax1,s3_) =zero;
         w(ixOmin1:ixOmax1,lfac_)=one
       end where
     else
       where((w(ixOmin1:ixOmax1,d_) < minrho).or.(w(ixOmin1:ixOmax1,tau_) &
          < smalltau).and. .not.patchw(ixOmin1:ixOmax1))
          patchierror(ixOmin1:ixOmax1)=-1
       elsewhere
          patchierror(ixOmin1:ixOmax1)=0
       end where
       call correctaux_usr(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,patchierror,&
          subname)
       call correctaux(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,patchierror,&
          subname)
     end if
  end if ! strict
end if


end subroutine smallvalues
!=============================================================================
subroutine getaux(clipping,w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,subname)

! Calculate the auxiliary variables lfac and xi within ixO^L
! clipping is not used (yet) 

include 'amrvacdef.f'

logical, intent(in)            :: clipping
integer, intent(in)            :: ixImin1,ixImax1, ixOmin1,ixOmax1
double precision, intent(inout):: w(ixImin1:ixImax1,1:nw)
double precision, intent(in)   :: x(ixImin1:ixImax1,1:ndim)
character(len=*), intent(in)   :: subname

double precision:: divb(ixGlo1:ixGhi1)
integer::          err,ix1, i1
double precision :: sold1_,sold2_,sold3_,bold1_,bold2_,bold3_, lfacold, xiold
integer          :: patchierror(ixGlo1:ixGhi1)=0
logical          :: patchw(ixGlo1:ixGhi1)

!-----------------------------------------------------------------------------

patchw(ixOmin1:ixOmax1)=.false.
if (fixsmall) call smallvalues(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,patchw,&
   subname)

! we compute auxiliaries lfac,xi from D,S,tau,B
! put the lfac and xi in the auxiliary fields lfac_ and xi_
do ix1= ixOmin1,ixOmax1
    lfacold = w(ix1,lfac_)
    xiold = w(ix1,xi_)
     sold1_=w(ix1,s1_); sold2_=w(ix1,s2_); sold3_=w(ix1,s3_);
     bold1_=w(ix1,b1_); bold2_=w(ix1,b2_); bold3_=w(ix1,b3_);

    call con2prim(w(ix1,lfac_),w(ix1,xi_), &
             w(ix1,d_),w(ix1,s1_),w(ix1,s2_),w(ix1,s3_),w(ix1,tau_),w(ix1,&
                b1_),w(ix1,b2_),w(ix1,b3_),w(ix1,b1_),w(ix1,b2_),w(ix1,b3_),&
                err)


    patchierror(ix1)=err
    if (err/=0.and.strictgetaux) then
       write(*,*),'Getaux error:',err,'ix^D=',ix1
       write(*,*),'input lfac=',lfacold,'s=',sold1_,sold2_,sold3_,'xi=',xiold,&
          'b=',bold1_,bold2_,bold3_
       write(*,*),'New ','lfac=',w(ix1,lfac_),'s=',w(ix1,s1_),w(ix1,s2_),&
          w(ix1,s3_),'xi=',w(ix1,xi_)
       write(*,*),'B=',w(ix1,b1_),w(ix1,b2_),w(ix1,b3_)
       write(*,*),'E=',w(ix1,e1_),w(ix1,e2_),w(ix1,e3_)
       write(*,*),'igrid =', saveigrid
       write(*,*),'level:', node(plevel_,saveigrid)
       write(*,*),'position  ', x(ix1, 1:ndim)
       write(*,*),'Called from: ',subname
       call mpistop("problem in getaux: retry with strictgetaux=F")
    end if
enddo

     ! first try user-defined method:
    if(.not.strictgetaux.and.any(patchierror(ixOmin1:ixOmax1)&
       /=0)) call correctaux_usr(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,&
       patchierror,subname)
     ! then, if any patchierror was not reset to 0, use default:
    if(.not.strictgetaux.and.any(patchierror(ixOmin1:ixOmax1)&
       /=0)) call correctaux(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,patchierror,&
       subname)

end subroutine getaux
!=============================================================================
subroutine correctaux(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,patchierror,subname)

include 'amrvacdef.f'

integer, intent(in)            :: ixImin1,ixImax1, ixOmin1,ixOmax1
integer, intent(in)            :: patchierror(ixGlo1:ixGhi1)
character(len=*), intent(in)   :: subname
double precision, intent(inout):: w(ixImin1:ixImax1,1:nw)
double precision, intent(in)   :: x(ixImin1:ixImax1,1:ndim)

integer        :: iw, kxOmin1,kxOmax1, ix1, i
logical        :: patchw(ixGlo1:ixGhi1)
!-----------------------------------------------------------------------------

do ix1= ixOmin1,ixOmax1
    if (patchierror(ix1)/=0) then
        do i=1,nflatgetaux
           kxOmin1= max(ix1-i,ixOmin1);
            kxOmax1= min(ix1+i,ixOmax1);
           if (any(patchierror(kxOmin1:kxOmax1)==0)) exit
        end do
        if (any(patchierror(kxOmin1:kxOmax1)==0))then
            ! in contrast to srhd case: always switch to primitive
            ! as error can become large when this is not done
            patchw(kxOmin1:kxOmax1)=(patchierror(kxOmin1:kxOmax1)/=0)
            call primitiven(ixImin1,ixImax1,kxOmin1,kxOmax1,w,patchw)
            do iw = 1,nw
               ! in contrast to srhd: do not alter magnetic field
               if (iw/=b1_.and.iw/=b2_.and.iw/=b3_ .and. iw/=e1_.and.iw&
                  /=e2_.and.iw/=e3_) then
                   w(ix1,iw)=sum(w(kxOmin1:kxOmax1,iw),patchierror&
                      (kxOmin1:kxOmax1)==0)/count(patchierror&
                      (kxOmin1:kxOmax1)==0)
               end if  
            end do
            patchw(ix1)=.false.
            call conserven(ixImin1,ixImax1,kxOmin1,kxOmax1,w,patchw)
        else
           ! no cells without error were found in cube of size nflatgetaux
           ! --> point of no recovery
           write(*,*),'Getaux error:',patchierror(ix1),'ix^D=',ix1
           write(*,*),'Called from: ',subname
           write(*,*),'New ','d=',w(ix1,d_),'s=',w(ix1,s1_),w(ix1,s2_),w(ix1,&
              s3_),'tau=',w(ix1,tau_)
           write(*,*),'B=',w(ix1,b1_),w(ix1,b2_),w(ix1,b3_)
           write(*,*),'E=',w(ix1,e1_),w(ix1,e2_),w(ix1,e3_)
           write(*,*),'lfac=',w(ix1,lfac_),'xi=',w(ix1,xi_)


           write(*,*),'position  ', x(ix1, 1:ndim)
           if (patchierror(ix1)<0) then
              call mpistop("---------correctaux from smallvalues-----")
           else
              call mpistop("---------correctaux from getaux----------")
           end if
        end if
    end if
enddo

end subroutine correctaux
!=============================================================================
! end module correctaux.t
!#############################################################################
