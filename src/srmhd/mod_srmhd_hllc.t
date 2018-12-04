module mod_srmhd_hllc
  use mod_srmhd_parameters
  use mod_srmhd_phys 
  implicit none
  private

  public :: srmhd_hllc_init

contains

  subroutine srmhd_hllc_init()
    use mod_physics_hllc

    phys_diffuse_hllcd => srmhd_diffuse_hllcd
    phys_get_lCD       => srmhd_get_lCD
    phys_get_wCD       => srmhd_get_wCD

  end subroutine srmhd_hllc_init

  subroutine srmhd_diffuse_hllcd(ixI^L,ixO^L,idim,wLC,wRC,fLC,fRC,patchf)
  ! made by Z. MELIANI 14/02/2018
  ! when method is hllcd or hllcd1 then: 
  ! this subroutine is to enforce regions where we AVOID HLLC
  ! and use TVDLF instead: this is achieved by setting patchf to 4 in
  ! certain regions. An additional input parameter is nxdiffusehllc
  ! which sets the size of the fallback region.
    use mod_global_parameters
    
    integer, intent(in)                                     :: ixI^L,ixO^L,idim
    double precision, dimension(ixI^S,1:nw), intent(in)     :: wRC,wLC
    double precision, dimension(ixI^S,1:nwflux),intent(in)  :: fLC, fRC
    integer, dimension(ixI^S), intent(inout) :: patchf
    
    integer :: ixOO^D,TxOO^L

    
    ! In a user-controlled region around any point with flux sign change between
    ! left and right, ensure fallback to TVDLF
    {do ixOO^D= ixO^LIM^D\}
      {
      TxOOmin^D= max(ixOO^D - nxdiffusehllc*kr(idim,^D), ixOmin^D);
      TxOOmax^D= min(ixOO^D + nxdiffusehllc*kr(idim,^D), ixOmax^D);
      \}
      if(abs(patchf(ixOO^D)) == 1 .or. abs(patchf(ixOO^D)) == 4)Then
         if(any(fRC(ixOO^D,1:nwflux)*fLC(ixOO^D,1:nwflux)<-smalldouble))Then
           where(Abs(patchf(TxOO^S))==1)
             patchf(TxOO^S) = 4
           endwhere
         endif
      endif
    {enddo^D&\}
  
  end subroutine srmhd_diffuse_hllcd

  subroutine srmhd_get_lCD(wLC,wRC,fLC,fRC,cmin,cmax,idim,ixI^L,ixO^L, &
                    whll,Fhll,lambdaCD,patchf)
  ! made by Z. MELIANI 14/02/2018
  
  ! Calculate lambda at CD and set the patchf to know the orientation
  ! of the riemann fan and decide on the flux choice
  ! We also compute here the HLL flux and w value, for fallback strategy
  
    use mod_global_parameters
    implicit none    

    integer, intent(in)                                     :: ixI^L,ixO^L,idim
    double precision, dimension(ixI^S,1:nw), intent(in)     :: wLC,wRC
    double precision, dimension(ixI^S,1:nwflux), intent(in) :: fLC,fRC
    double precision, dimension(ixI^S), intent(in)          :: cmax,cmin
    integer         , dimension(ixI^S), intent(inout)       :: patchf
    double precision, dimension(ixI^S,1:nwflux), intent(out):: Fhll,whll
    double precision, dimension(ixI^S), intent(out)         :: lambdaCD
    

    double precision, dimension(ixG^T,mom(1):mom(ndir))     :: vCD
    double precision, dimension(ixG^T)     :: Aco,Bco,Cco,Delta, &
                                              BpdotfBp, Bp2,fBp2
    logical         , dimension(ixI^S)     :: Cond_patchf, Cond_Bidimhll
    double precision                       :: Epsilon
    integer                                :: iw

    ! on entry, patch is preset to contain values from -2,1,2,4
    !      -2: take left flux, no computation here
    !      +2: take right flux, no computation here
    !      +4: take TVDLF flux, no computation here
    !       1: compute the characteristic speed for the CD
    
    Cond_patchf(ixO^S)=abs(patchf(ixO^S))==1
    
    do iw=1,nwflux
      where(Cond_patchf(ixO^S))
      !============= compute HLL flux ==============!
      Fhll(ixO^S,iw)= (cmax(ixO^S)*fLC(ixO^S,iw)-cmin(ixO^S)*fRC(ixO^S,iw) &
                       + cmin(ixO^S)*cmax(ixO^S)*(wRC(ixO^S,iw)-wLC(ixO^S,iw)))&
                       /(cmax(ixO^S)-cmin(ixO^S))
      !======== compute intermediate HLL state =======!
      whll(ixO^S,iw) = (cmax(ixO^S)*wRC(ixO^S,iw)-cmin(ixO^S)*wLC(ixO^S,iw)&
                        +fLC(ixO^S,iw)-fRC(ixO^S,iw))/(cmax(ixO^S)-cmin(ixO^S))
      endwhere
    enddo


!Calculate the Characteristic speed at the contact
! part specific to SRMHD and HLLC
! Eq. 41, Mignone & Bodo, MNRAS 2006
!  write this eq as:  Aco lambda^2 - Bco lambda + Cco = 0
where(Cond_patchf(ixO^S))
  ! only store the HD part first, sufficient for case where normal B vanishes
  Aco(ixO^S)    = Fhll(ixO^S,e_)+Fhll(ixO^S,d_)
  Bco(ixO^S)    = Fhll(ixO^S,mom(idim))+whll(ixO^S,e_)+whll(ixO^S,d_)
  Cco(ixO^S)    = whll(ixO^S,mom(idim))
endwhere

! condition on the normal magnetic field
Cond_Bidimhll(ixO^S) = (dabs(whll(ixO^S,mag(idim)))<=smalldouble)

! Case With Normal Magnetic field
if(any(.not.Cond_Bidimhll(ixO^S).and.Cond_patchf(ixO^S)))then
  where(.not.Cond_Bidimhll(ixO^S).and.Cond_patchf(ixO^S))
   !---- Initialisation ----!
    BpdotfBp(ixO^S) = zero
    Bp2(ixO^S) = zero
    fBp2(ixO^S) = zero
  endwhere

  !The calculation of the Transverse components part
  do iw=mag(1),mag(ndir)
    if(iw /= mag(idim))then
      where(.not.Cond_Bidimhll(ixO^S) .and. Cond_patchf(ixO^S))
        BpdotfBp(ixO^S) = BpdotfBp(ixO^S) + whll(ixO^S,iw)*Fhll(ixO^S,iw)
        Bp2(ixO^S) = Bp2(ixO^S) + whll(ixO^S,iw)**2.0d0
        fBp2(ixO^S) = fBp2(ixO^S) + Fhll(ixO^S,iw)**2.0d0
      endwhere
     endif
    enddo

    where(.not.Cond_Bidimhll(ixO^S) .and. Cond_patchf(ixO^S))
     Aco(ixO^S)    = Aco(ixO^S) - BpdotfBp(ixO^S)
     Bco(ixO^S)    = Bco(ixO^S)- Bp2(ixO^S) - fBp2(ixO^S)
     Cco(ixO^S)    = Cco(ixO^S) - BpdotfBp(ixO^S)
    endwhere
   endif


   where(Cond_patchf(ixO^S))
    Delta(ixO^S) = Bco(ixO^S)**2.0d0- 4.0d0*Aco(ixO^S) * Cco(ixO^S)
    where(Aco(ixO^S)/=zero .and. Delta(ixO^S)>=zero)
     !Calculate the Characteristic speed at the contact
     ! only the minus sign is between [-1,1]
     lambdaCD(ixO^S) = (Bco(ixO^S) - dsqrt(Delta(ixO^S)))/(2.0d0*Aco(ixO^S))
    elsewhere(Aco(ixO^S)==zero .and.  Bco(ixO^S)/=zero)
     lambdaCD(ixO^S) = Cco(ixO^S)/Bco(ixO^S)
    elsewhere(Delta(ixO^S)<zero)
     lambdaCD(ixO^S) = zero
     ! we will fall back to HLL flux case in this degeneracy
     patchf(ixO^S) =  3
    endwhere
   endwhere

   where(patchf(ixO^S)==3)
    Cond_patchf(ixO^S)=.false.
   end where

   where(Cond_patchf(ixO^S))
    ! double check whether obtained speed is in between min and max speeds given
    ! and identify in which part of the Riemann fan the time-axis is
    where(cmin(ixO^S)<zero.and.lambdaCD(ixO^S)>zero&
          .and.lambdaCD(ixO^S)<cmax(ixO^S))
     patchf(ixO^S) = -1
    elsewhere(cmax(ixO^S)>zero.and.lambdaCD(ixO^S)<zero&
              .and.lambdaCD(ixO^S)>cmin(ixO^S))
     patchf(ixO^S) =  1
    elsewhere(lambdaCD(ixO^S)>=cmax(ixO^S).or.lambdaCD(ixO^S) <= cmin(ixO^S))
     lambdaCD(ixO^S) = zero
     ! we will fall back to HLL flux case in this degeneracy
     patchf(ixO^S) =  3
    endwhere
   endwhere

   where(patchf(ixO^S)== 3)
    Cond_patchf(ixO^S)=.false.
   end where


   ! handle the specific case where the time axis is exactly on the CD
   if(any(lambdaCD(ixO^S)==zero.and.Cond_patchf(ixO^S)))then
    !determine which sector (forward or backward) of the Riemann fan is smallest
    ! and select left or right flux accordingly
    where(lambdaCD(ixO^S)==zero.and.Cond_patchf(ixO^S))
     where(-cmin(ixO^S)>=cmax(ixO^S))
      patchf(ixO^S) =  1
     elsewhere
      patchf(ixO^S) = -1
     endwhere
    endwhere
   endif



   ! eigenvalue lambda for contact is near zero: decrease noise by this trick
   if(flathllc)then
    Epsilon=1.0d-6
    where(Cond_patchf(ixO^S).and. &
     dabs(lambdaCD(ixO^S))/max(cmax(ixO^S),Epsilon)< Epsilon  .and. &
     dabs(lambdaCD(ixO^S))/max(dabs(cmin(ixO^S)),Epsilon)< Epsilon)
     lambdaCD(ixO^S) =  zero
    end where
   end if
    

   if(any(dabs(lambdaCD(ixO^S))>1.0d0 .and. Cond_patchf(ixO^S)))then
    call mpistop("problems with lambdaCD>1")
   endif

   ! next should never happen
   if(any(patchf(ixO^S)==0))then
    call mpistop("patchf=0")
   endif


  end subroutine srmhd_get_lCD

  subroutine srmhd_get_wCD(x,wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,&
                    cmin,cmax,&
                    ixI^L,ixO^L,idim,f)
  ! made by Z. MELIANI 14/02/2018
  ! compute the intermediate state U*
  ! only needed where patchf=-1/1
  
  ! reference Li S., JCP, 203, 2005, 344-357
  ! reference T. Miyoski, Kusano JCP, 2008, 2005.
    use mod_global_parameters
    
    integer, intent(in)                                      :: ixI^L,ixO^L,idim
    double precision, dimension(ixI^S,1:nw), intent(in)      :: wRC,wLC
    double precision, dimension(ixI^S,1:ndim), intent(in)    :: x
    double precision, dimension(ixI^S,1:nwflux), intent(in)  :: whll, Fhll
    double precision, dimension(ixI^S), intent(in)           :: lambdaCD
    double precision, dimension(ixI^S), intent(in)           :: cmax,cmin
    double precision, dimension(ixI^S,1:nwflux), intent(in)  :: fRC,fLC
    double precision, dimension(ixI^S,1:nwflux),intent(out)  :: f
    double precision, dimension(ixI^S,1:nw)         :: wCD,wSub
    double precision, dimension(ixI^S,1:nwflux)     :: fSub
    double precision, dimension(ixI^S)              :: vSub,cspeed,pCD,VdotBCD
    integer         , dimension(ixI^S),intent(inout):: patchf

    double precision, dimension(ixG^T,mom(1):mom(ndir))      :: vCD
    double precision, dimension(ixI^S)             :: Ratio_CD
    double precision, dimension(ixI^S)             :: vRC,vLC
    logical         , dimension(ixI^S)             :: Cond_Bidimhll
    integer                                        :: n, iw, idir,ix^D
    double precision, dimension(ixO^S)             :: VdotB, B2

    !-------------- auxiliary Speed and array-------------!
    call srmhd_get_v_idim(wRC,x,ixI^L,ixO^L,idim,vRC)
    call srmhd_get_v_idim(wLC,x,ixI^L,ixO^L,idim,vLC)


    where(patchf(ixO^S) == 1)
     cspeed(ixO^S) = cmax(ixO^S)
     vSub(ixO^S)   = vRC(ixO^S)
    elsewhere(patchf(ixO^S) == -1)
     cspeed(ixO^S) = cmin(ixO^S)
     vSub(ixO^S)   = vLC(ixO^S)
    endwhere

    do iw=1,nwflux
     if(iw /= mag(idim))then
      where(patchf(ixO^S) == 1)
       wSub(ixO^S,iw) =  wRC(ixO^S,iw)
       fSub(ixO^S,iw) =  fRC(ixO^S,iw)
      elsewhere(patchf(ixO^S) == -1)
       wSub(ixO^S,iw) =  wLC(ixO^S,iw)
       fSub(ixO^S,iw) =  fLC(ixO^S,iw)
      endwhere
     endif
    enddo

    wSub(ixO^S,mag(idim)) = whll(ixO^S,mag(idim))
    Cond_Bidimhll(ixO^S) = (dabs(wSub(ixO^S,mag(idim)))<=smalldouble)

    where(abs(patchf(ixO^S))==1)
     Ratio_CD(ixO^S) = (cspeed(ixO^S)-vSub(ixO^S))&
                       /(cspeed(ixO^S)-lambdaCD(ixO^S))
     wCD(ixO^S,d_)   = wSub(ixO^S,d_)*Ratio_CD(ixO^S)
    end where
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      if(abs(patchf(ix^D))==1) then
       do n=1,srmhd_n_tracer
          iw = tracer(n)
          wCD(ix^D,iw) = wSub(ix^D,iw)*Ratio_CD(ix^D)
        end do
      end if
    {end do\}







     ! in case we have somewhere a normal component, need to distinguish
     bxnozero : if(any(.not.Cond_Bidimhll(ixO^S)))then

     !==== Magnetic field ====!
     do iw =mag(1),mag(ndir)
      if(iw /= mag(idim))then
       ! Transverse components
       where(.not.Cond_Bidimhll(ixO^S)  .and. abs(patchf(ixO^S))==1)
         ! case from eq 37
         wCD(ixO^S,iw) = whll(ixO^S,iw)
       elsewhere(Cond_Bidimhll(ixO^S) .and. abs(patchf(ixO^S))==1)
         ! case from eq 53
         wCD(ixO^S,iw) = wSub(ixO^S,iw) * Ratio_CD(ixO^S)
       endwhere
      else  !  Normal component
       wCD(ixO^S,iw) = wSub(ixO^S,iw)
      endif
     enddo

     !====== velocities ========!
     do idir = 1,ndir
      if(idir /=idim)then
       ! Transverse components
       where(.not.Cond_Bidimhll(ixO^S)  .and. abs(patchf(ixO^S))==1)
        ! case from eq 38
        vCD(ixO^S,mom(idir))=(wCD(ixO^S,mag(idir))*lambdaCD(ixO^S)&
                       -Fhll(ixO^S,mag(idir)))/wCD(ixO^S,mag(idim))
       elsewhere(Cond_Bidimhll(ixO^S)  .and. abs(patchf(ixO^S))==1)
        ! unused case
        vCD(ixO^S,mom(idir))=zero
       endwhere
      else ! Normal component
       where(abs(patchf(ixO^S))==1)
        vCD(ixO^S,mom(idir))=lambdaCD(ixO^S)
       endwhere
      endif
     enddo
     ! enforce fallback strategy for case where characteristic velocity 
     !  at contact
     ! is unphysical
     where(.not. Cond_Bidimhll(ixO^S)  &
           .and. sum(vCD(ixO^S,mom(:))**2.0d0,dim=ndim+1) > 1.0d0 &
           .and. abs(patchf(ixO^S))==1)
      patchf(ixO^S) = 3
     endwhere

     where(.not.Cond_Bidimhll(ixO^S).and. abs(patchf(ixO^S))==1)
      wCD(ixO^S,lfac_) = 1.0d0/dsqrt(1.0d0-sum(vCD(ixO^S,mom(:))**2.0d0,dim=ndim+1))
      VdotB(ixO^S) = sum(vCD(ixO^S,mom(:))*wCD(ixO^S,mag(:)),dim=ndim+1)
      !--- total Pressure from eq 40 ---!
      pCD(ixO^S)  = -lambdaCD(ixO^S)&
           *(Fhll(ixO^S,e_)+Fhll(ixO^S,d_)-VdotB(ixO^S)*wCD(ixO^S,mag(idim)))&
           +Fhll(ixO^S,mom(idim))+(wCD(ixO^S,mag(idim))/wCD(ixO^S,lfac_))**2.0d0

     elsewhere(Cond_Bidimhll(ixO^S) .and. abs(patchf(ixO^S))==1)
      !------ total Pressure from 48 ------!
      pCD(ixO^S)  = -(Fhll(ixO^S,e_)+Fhll(ixO^S,d_))*lambdaCD(ixO^S) &
                  + Fhll(ixO^S,mom(idim))
     endwhere

     ! enforce fallback strategy for case where total pressure at contact
     ! is unphysical (should never happen?)
     where(pCD(ixO^S) <= zero.and. abs(patchf(ixO^S))==1)
      patchf(ixO^S) = 3
     endwhere
     !------- Momentum ------!
     do iw =mom(1),mom(ndir)
      if(iw /= mom(idim))then
       where(.not.Cond_Bidimhll(ixO^S)   .and.  abs(patchf(ixO^S))==1)
        ! eq. 44 45
        wCD(ixO^S,iw) = (cspeed(ixO^S)*wSub(ixO^S,iw)-fSub(ixO^S,iw)&
                         -wCD(ixO^S,mag(idim))*&
                   (wCD(ixO^S,iw-mom(1)+mag(1))/wCD(ixO^S,lfac_)**2.0d0&
           +VdotB(ixO^S) * vCD(ixO^S,iw))) /(cspeed(ixO^S)-lambdaCD(ixO^S))
       elsewhere(Cond_Bidimhll(ixO^S) .and. abs(patchf(ixO^S))==1)
        ! eq. 51
        wCD(ixO^S,iw) = wSub(ixO^S,iw) * Ratio_CD(ixO^S)
       endwhere
      endif
     enddo

     where(.not.Cond_Bidimhll(ixO^S)   .and. abs(patchf(ixO^S))==1)
      !---- Tau Right combine 46 43----!
      wCD(ixO^S,e_) = (cspeed(ixO^S) * wSub(ixO^S,e_) &
                 + vSub(ixO^S) * wSub(ixO^S,d_) - wSub(ixO^S,mom(idim))&
                 +lambdaCD(ixO^S)*pCD(ixO^S)-VdotB(ixO^S)*wCD(ixO^S,mag(idim)))&
                 /(cspeed(ixO^S)-lambdaCD(ixO^S))
      !--- Sidim Right from eq 33 ---!
      wCD(ixO^S,mom(idim)) = (wCD(ixO^S,e_)+wCD(ixO^S,d_)+pCD(ixO^S))&
                  *lambdaCD(ixO^S)-VdotB(ixO^S)*wCD(ixO^S,mag(idim))
     elsewhere(Cond_Bidimhll(ixO^S)  .and. abs(patchf(ixO^S))==1)
      !---- Tau Right combine 50 52----!
      wCD(ixO^S,e_) = (cspeed(ixO^S) * wSub(ixO^S,e_) &
                  + vSub(ixO^S) * wSub(ixO^S,d_) - wSub(ixO^S,mom(idim))&
                  +lambdaCD(ixO^S)*pCD(ixO^S))&
                  /(cspeed(ixO^S)-lambdaCD(ixO^S))
      !--- sidim Right from 49 ---!
      wCD(ixO^S,mom(idim)) = (wCD(ixO^S,e_)+wCD(ixO^S,d_)+pCD(ixO^S))&
                  *lambdaCD(ixO^S)
      !-- Not real lfac, but fill it anyway --!
      wCD(ixO^S,lfac_) = 1.0d0
      !-------------------!
     endwhere
else ! bxnozero
  ! in case we have everywhere NO normal B component, no need to distinguish
  !---- Magnetic field ----!
  do iw =mag(1),mag(ndir)
    if(iw /= mag(idim))then
      !Transverse components
      where(abs(patchf(ixO^S))==1)
        wCD(ixO^S,iw) = wSub(ixO^S,iw) * Ratio_CD(ixO^S)
      endwhere
    else   !Normal components
      where(abs(patchf(ixO^S))==1)
        wCD(ixO^S,iw) = wSub(ixO^S,iw)
      endwhere
    endif
  enddo

  !---- momentum ----!
  do iw =mom(1),mom(ndir)
    if(iw /= mom(idim))Then
      where(abs(patchf(ixO^S))==1)
        wCD(ixO^S,iw) = wSub(ixO^S,iw) * Ratio_CD(ixO^S)
      endwhere
    endif
  enddo

  where(abs(patchf(ixO^S))==1)
   !------- Pressure -------!
   pCD(ixO^S)      = Fhll(ixO^S,mom(idim))&
                    - (Fhll(ixO^S,e_)+Fhll(ixO^S,d_))* lambdaCD(ixO^S)
   !---- Tau Right ---!
   wCD(ixO^S,e_) = (cspeed(ixO^S) * wSub(ixO^S,e_) &
                    + vSub(ixO^S) * wSub(ixO^S,d_) - wSub(ixO^S,mom(idim))&
                    +lambdaCD(ixO^S)*pCD(ixO^S))/(cspeed(ixO^S)-lambdaCD(ixO^S))
   !---- Sx_ Right ---!
   wCD(ixO^S,mom(idim)) =(wCD(ixO^S,e_)+wCD(ixO^S,d_)+pCD(ixO^S))&
                         *lambdaCD(ixO^S)
   !-- Not real lfac, but fill it anyway --!
   wCD(ixO^S,lfac_) = 1.0d0
  endwhere

 ! In case of B_idim = 0, we need only the normal velocity
  do iw = mom(1),mom(ndir)
    if(iw /= mom(idim))Then
      where(abs(patchf(ixO^S))==1)
        vCD(ixO^S,iw) = zero
      endwhere
    else
      where(abs(patchf(ixO^S))==1)
        vCD(ixO^S,iw) = lambdaCD(ixO^S)
      endwhere
    endif
  enddo

endif bxnozero

do iw=1,nwflux
 if(iw /= mag(idim).and. iw/=psi_)then
   where(abs(patchf(ixO^S))==1)
       ! f_i=fSub+lambda (wCD-wSub)
       f(ixO^S,iw)=fSub(ixO^S,iw)+cspeed(ixO^S)*(wCD(ixO^S,iw)-wsub(ixO^S,iw))
   endwhere
 else
       f(ixO^S,iw)=zero
 end if
end do

  end subroutine srmhd_get_wCD

end module mod_srmhd_hllc
