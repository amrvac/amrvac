module mod_mhd_hllc
  use mod_mhd_phys

  implicit none
  private

  public :: mhd_hllc_init

contains

  subroutine mhd_hllc_init()
    use mod_physics_hllc

    phys_diffuse_hllcd => mhd_diffuse_hllcd
    phys_get_lCD => mhd_get_lCD
    phys_get_wCD => mhd_get_wCD

  end subroutine mhd_hllc_init

  subroutine mhd_diffuse_hllcd(ixI^L,ixO^L,idim,wLC,wRC,fLC,fRC,patchf)
  ! when method is hllcd or hllcd1 then: 
  ! this subroutine is to enforce regions where we AVOID HLLC
  ! and use TVDLF instead: this is achieved by setting patchf to 4 in
  ! certain regions. An additional input parameter is nxdiffusehllc
  ! which sets the size of the fallback region.
    use mod_global_parameters
    
    integer, intent(in)                                      :: ixI^L,ixO^L,idim
    double precision, dimension(ixI^S,1:nw), intent(in)      :: wRC,wLC
    double precision, dimension(ixI^S,1:nwflux),intent(in) :: fLC, fRC
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
  
  end subroutine mhd_diffuse_hllcd

  subroutine mhd_get_lCD(wLC,wRC,fLC,fRC,cmin,cmax,idim,ixI^L,ixO^L, &
                    whll,Fhll,lambdaCD,patchf)
  
  ! Calculate lambda at CD and set the patchf to know the orientation
  ! of the riemann fan and decide on the flux choice
  ! We also compute here the HLL flux and w value, for fallback strategy
  
    use mod_global_parameters
    
    integer, intent(in)                                      :: ixI^L,ixO^L,idim
    double precision, dimension(ixI^S,1:nw), intent(in)      :: wLC,wRC
    double precision, dimension(ixI^S,1:nwflux), intent(in):: fLC,fRC
    double precision, dimension(ixI^S), intent(in)           :: cmax,cmin
    integer         , dimension(ixI^S), intent(inout)        :: patchf
    double precision, dimension(ixI^S,1:nwflux), intent(out) :: Fhll,whll
    double precision, dimension(ixI^S), intent(out)            :: lambdaCD
    
    logical         , dimension(ixI^S)     :: Cond_patchf
    double precision                       :: Epsilon
    integer                                :: iw

    ! on entry, patch is preset to contain values from -2,1,2,4
    !      -2: take left flux, no computation here
    !      +2: take right flux, no computation here
    !      +4: take TVDLF flux, no computation here
    !       1: compute the characteristic speed for the CD
    
    Cond_patchf(ixO^S)=(abs(patchf(ixO^S))==1)
    
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
    
    ! deduce the characteristic speed at the CD
    where(Cond_patchf(ixO^S))
      lambdaCD(ixO^S)=whll(ixO^S,mom(idim))/whll(ixO^S,rho_)
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
      ! determine which sector (forward or backward) of the Riemann fan is smallest
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

  end subroutine mhd_get_lCD

  subroutine mhd_get_wCD(wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,cmin,cmax,&
                    ixI^L,ixO^L,idim,f)
  ! compute the intermediate state U*
  ! only needed where patchf=-1/1
  
  ! reference Li S., JCP, 203, 2005, 344-357
  ! reference T. Miyoski, Kusano JCP, 2008, 2005.
    use mod_global_parameters
    
    integer, intent(in)                                      :: ixI^L,ixO^L,idim
    double precision, dimension(ixI^S,1:nw), intent(in)      :: wRC,wLC
    double precision, dimension(ixI^S,1:nwflux), intent(in):: whll, Fhll
    double precision, dimension(ixI^S), intent(in)           :: lambdaCD
    double precision, dimension(ixI^S), intent(in)           :: cmax,cmin
    double precision, dimension(ixI^S,1:nwflux), intent(in):: fRC,fLC
    double precision, dimension(ixI^S,1:nwflux),intent(out):: f
    double precision, dimension(ixI^S,1:nw)        :: wCD,wSub
    double precision, dimension(ixI^S,1:nwflux)    :: fSub
    double precision, dimension(ixI^S)             :: vSub,cspeed,pCD,VdotBCD
    integer         , dimension(ixI^S), intent(in)           :: patchf

    integer                                        :: n, iw, idir,ix^D

    !-------------- auxiliary Speed and array-------------!
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
       if(patchf(ix^D)==1) then
         cspeed(ix^D)=cmax(ix^D)
         vSub(ix^D)=wRC(ix^D,mom(idim))/wRC(ix^D,rho_)
         wSub(ix^D,:)=wRC(ix^D,:)
         fSub(ix^D,:)=fRC(ix^D,:)
       else if(patchf(ix^D)==-1) then
         cspeed(ix^D)=cmin(ix^D)
         vSub(ix^D)=wLC(ix^D,mom(idim))/wLC(ix^D,rho_)
         wSub(ix^D,:)=wLC(ix^D,:)
         fSub(ix^D,:)=fLC(ix^D,:)
       end if
    {end do\}
    
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      if(abs(patchf(ix^D))==1) then
        wCD(ix^D,rho_) = wSub(ix^D,rho_)&
                         *(cspeed(ix^D)-vSub(ix^D))/(cspeed(ix^D)-lambdaCD(ix^D))
        do n=1,mhd_n_tracer
          iw = tracer(n)
          wCD(ix^D,iw) = wSub(ix^D,iw)*(cspeed(ix^D)-vSub(ix^D))&
             /(cspeed(ix^D)-lambdaCD(ix^D))
        end do
        !==== Magnetic field ====!
        do idir=1,ndir
          ! case from eq 31
          wCD(ix^D,mag(idir)) = whll(ix^D,mag(idir))
        end do
        !------- Momentum ------!
        do iw=1, ndir
          if(iw /= idim)then
            ! eq. 21 22
            wCD(ix^D,mom(iw))=(cspeed(ix^D)*wSub(ix^D,mom(iw))-fSub(ix^D,mom(iw))&
                -wCD(ix^D,mag(idim))*wCD(ix^D,mag(iw)))/&
                (cspeed(ix^D)-lambdaCD(ix^D))
          else
            ! eq. 20
            wCD(ix^D,mom(iw)) =  wCD(ix^D,rho_) * lambdaCD(ix^D)
          endif
        enddo
        if(mhd_energy) then
          VdotBCD(ix^D) = sum(whll(ix^D,mom(:))*whll(ix^D,mag(:)))/whll(ix^D,rho_)
          ! Eq 17
          pCD(ix^D)  = wsub(ix^D,rho_)*(cspeed(ix^D)-vSub(ix^D))&
                        *(lambdaCD(ix^D)-vSub(ix^D))&
                        +fSub(ix^D,mom(idim))-wsub(ix^D,mom(idim))*vSub(ix^D)&
                        + wCD(ix^D,mag(idim))**2
          ! Eq 31
          wCD(ix^D,e_) = (cspeed(ix^D)*wSub(ix^D,e_) &
                          -fSub(ix^D,e_)+lambdaCD(ix^D)*pCD(ix^D)&
                          -VdotBCD(ix^D)*wCD(ix^D,mag(idim)))&
                          /(cspeed(ix^D)-lambdaCD(ix^D))
        end if
      end if
    {end do\}

    do iw=1,nwflux
     if(iw == mag(idim)) then
       f(ixO^S,iw)=zero
     else if(mhd_glm .and. iw == psi_) then
       f(ixO^S,iw)=zero
     else
       where(abs(patchf(ixO^S))==1)
         ! f_i=fsub+lambda (wCD-wSub)
         f(ixO^S,iw)=fsub(ixO^S,iw)+cspeed(ixO^S)*(wCD(ixO^S,iw)-wsub(ixO^S,iw))
       endwhere
     end if
    end do

  end subroutine mhd_get_wCD

end module mod_mhd_hllc
