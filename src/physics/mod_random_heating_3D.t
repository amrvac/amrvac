module mod_random_heating_3D
  use mod_global_parameters
  use mod_comm_lib, only: mpistop

  double precision, allocatable :: va^D(:)
  double precision, allocatable :: dtimearr(:),timearray(:)
  ! This module is for imposing randomized heating in any-D setup
  ! Users can define the parameters "trelax,tramp,ntimes,vlim^D,periods,variation,nwaves" in the mod_usr.t file.
  !periods=300,variation=75  The time durations in seconds of each heating episode is typically 300s +/- 75s
  !ntimes = 100              100 heating episodes, e.g., 500 min. Enough for the simulation
  !nwaves = 1000             1000 waves for each episode
  !si=5./6.                  spectral index
  !vlim^D could be chosen as twice of the maximal spatial resolution in the simulation

contains

  subroutine getlQgrid(lQgrid,ixI^L,ixO^L,qt,trelax,tramp,x,ntimes,vlim^D)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L, ntimes, vlim^D
    double precision, intent(in)    :: qt, trelax, tramp,x(ixI^S,1:ndim)
    double precision, intent(inout) :: lQgrid(ixI^S)

    double precision                :: xcar(ixI^S,1:ndim)
    double precision                :: tshift,tr,stt^D,xcarmin^D
    double precision                :: tl1,tl2,tl3,tl4,tl5,tl6
    double precision                :: pulse1,pulse2,pulse3,pulse4,pulse5,pulse6
    integer                         :: iip^D,ix^D,i

    tshift = qt - trelax
    
    if(qt .lt. trelax) then
      tr=zero
    elseif(qt .lt. trelax+tramp) then
      tr=(qt-trelax)/tramp
    else
      tr=one
    endif
    

    do i = 4, ntimes
        if(qt .lt. trelax + timearray(i)) then
            tl1 = dexp(-(tshift - timearray(i-3))**2/(0.5d0*dtimearr(i-3))**2)
            tl2 = dexp(-(tshift - timearray(i-2))**2/(0.5d0*dtimearr(i-2))**2)
            tl3 = dexp(-(tshift - timearray(i-1))**2/(0.5d0*dtimearr(i-1))**2)
            tl4 = dexp(-(tshift - timearray(i)  )**2/(0.5d0*dtimearr(i)  )**2)
            tl5 = dexp(-(tshift - timearray(i+1))**2/(0.5d0*dtimearr(i+1))**2)
            tl6 = dexp(-(tshift - timearray(i+2))**2/(0.5d0*dtimearr(i+2))**2)
            {do ix^DB = ixOmin^DB,ixOmax^DB\}
              select case(coordinate)
              case (Cartesian)
                 ! remapped on the unit cube
                 ^D&stt^D=(xprobmax^D-xprobmin^D)/(dble(vlim^D)-1.d0);
                 ^D&iip^D=floor((x(ix^DD,^D)-xprobmin^D)/stt^D)+1;
                 ^D&iip^D=max(1,iip^D);
              case (cylindrical)
                  {^IFTWOD
                  call mpistop("random_heating_3D to be adjusted for 2D cylindrical")
                  }
                  {^IFTHREED
                  ! from R-Z-phi to x-y-z
                  xcar(ix1,ix2,ix3,1)=x(ix1,ix2,ix3,r_)*dcos(x(ix1,ix2,ix3,phi_))
                  xcar(ix1,ix2,ix3,2)=x(ix1,ix2,ix3,r_)*dsin(x(ix1,ix2,ix3,phi_))
                  xcar(ix1,ix2,ix3,3)=x(ix1,ix2,ix3,z_)
                  ! minimal coordinate x-y-z from Rmax and Z-range
                  xcarmin1=-xprobmax1
                  xcarmin2=-xprobmax1
                  xcarmin3=xprobmin2
                  ! x-coordinate between -Rmax,+Rmax
                  stt1=(2.0d0*xprobmax1)/(dble(vlim1)-1.d0)
                  ! y-coordinate between -Rmax,+Rmax
                  stt2=(2.0d0*xprobmax1)/(dble(vlim2)-1.d0)
                  ! z-coordinate between Zmin,Zmax
                  stt3=(xprobmax2-xprobmin2)/(dble(vlim3)-1.d0)
                  }
                 ^D&iip^D=floor((xcar(ix^DD,^D)-xcarmin^D)/stt^D)+1;
                 ^D&iip^D=max(1,iip^D);
              case(spherical)
                 call mpistop("random_heating_3D to be adjusted for this geometry")
              case default
                 call mpistop("random_heating_3D to be adjusted for this geometry")
              end select
              pulse1=(^D&va^D(iip^D+(i-4)*vlim^D)|*)
              pulse2=(^D&va^D(iip^D+(i-3)*vlim^D)|*)
              pulse3=(^D&va^D(iip^D+(i-2)*vlim^D)|*)
              pulse4=(^D&va^D(iip^D+(i-1)*vlim^D)|*)
              pulse5=(^D&va^D(iip^D+(i)*vlim^D)|*)
              pulse6=(^D&va^D(iip^D+(i+1)*vlim^D)|*)
              lQgrid(ix^D) = tl1 * pulse1 +  tl2 *pulse2  &
                            +tl3 * pulse3 +  tl4 *pulse4  &
                            +tl5 * pulse5 +  tl6 *pulse6
            {enddo\}
            exit 
        endif 
      enddo 

      lQgrid(ixO^S)=tr*lQgrid(ixO^S)
  
  end subroutine getlQgrid

  subroutine generateTV(ntimes,vlim^D,periods,variation,nwaves,si)
    integer, intent(in)          :: ntimes,vlim^D,nwaves
    double precision, intent(in) :: periods,variation,si
    integer                      :: i,vx^D,j

    ! generate random T, typically 300s +/- 75s   

    allocate(dtimearr(ntimes))  
    allocate(timearray(ntimes)) 

    if (mype==0) then 
      call randomT(timearray,ntimes,periods,variation)  
      dtimearr = timearray
      do i = 2, ntimes
        dtimearr(i) = timearray(i) - timearray(i-1)
      enddo
    endif   
    call MPI_BARRIER(icomm,ierrmpi)
    if (npe>1) then
      call MPI_BCAST(timearray,ntimes,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(dtimearr,ntimes,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    endif

    ! spatial distribution
  
    {^D& 
    vx^D=ntimes*vlim^D
    allocate(va^D(vx^D)) 
    if (mype==0) then 
       call randomV(^D,va^D,ntimes,vlim^D,nwaves,si)
    endif
    call MPI_BARRIER(icomm,ierrmpi)
    if (npe>1) then
      call MPI_BCAST(va^D,vx^D,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    endif 
    }

  end subroutine generateTV

  subroutine randomT(tb,ntimes,periods,variation)
    integer, intent(in)             :: ntimes
    double precision, intent(in)    :: periods,variation
    double precision, intent(inout) :: tb(:)

    character(len=100)  :: filename
    double precision                :: tt1,tt,mm
    integer                         :: i

    tt = 0.d0
    do i = 1, ntimes
      call random_number(mm)  
      tt1 = periods + variation * 2. * (mm - 0.5)
      tt = tt + tt1
      tb(i) = tt / unit_time
    end do
 
    write(filename,"(a)") "randomtimes.dat"
    open(unit=21,file=filename,form='formatted',status='replace')
    write(21,'(es12.4)') tb
    close(21)
     
  end subroutine randomT

  subroutine randomV(ival,va,ntimes,vlim,nwaves,si)
    integer, intent(in)             :: ival,ntimes,vlim,nwaves
    double precision, intent(in)    :: si
    double precision, intent(inout) :: va(:)

    character(len=100)  :: filename
    character(len=10)   :: xdirection

    double precision, allocatable   :: vn(:,:),vm(:)
    double precision                :: lambda,ampl,lambda_min,lambda_max,rms,phase
    integer                         :: i, j, kk
 
    allocate(vn(nwaves, vlim))
    allocate(vm(vlim))
    lambda_min = 1.d0/dble(vlim)
    lambda_max = 1.d0
    do kk = 1, ntimes
        do i = 1, nwaves
          lambda = (dble(i-1)/dble(nwaves-1))*(lambda_max-lambda_min)+lambda_min
          ampl = lambda**si 
          call random_number(phase)        
          ! phase to vary between -2pi-->+2pi
          phase=4.0d0*dpi*(phase-0.5d0)
          do j = 1, vlim  
            vn(i,j) = ampl*dsin(2.d0*dpi*(dble(j-1)/dble(vlim-1))/lambda+phase)
          end do
        end do     
        rms = 0.d0
        vm = 0.d0    
        do j = 1, vlim  
          ! sum all the waves in grid point j
          do i = 1, nwaves
            vm(j) = vm(j) + vn(i, j) 
          end do
          rms = rms + vm(j)**2
        end do         
        vm = vm / dsqrt(rms/vlim)    ! normalization        
        do j = 1, vlim
          va(vlim*(kk-1)+j) = vm(j)**2
        end do       
    end do

    write(xdirection,'(I2.2)') ival
    write(filename,"(a)") "randompositions"//TRIM(xdirection)//".dat"
    open(unit=22,file=filename,form='formatted',status='replace')
    write(22,'(es12.4)') va
    close(22)

    deallocate(vm)
    deallocate(vn)

  end subroutine randomV

end module mod_random_heating_3D
