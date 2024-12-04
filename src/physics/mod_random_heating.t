module mod_random_heating
  use mod_global_parameters
  use mod_comm_lib, only: mpistop
  double precision, allocatable :: va1(:),ta1(:),tb1(:),va2(:),ta2(:),tb2(:)
  !> This module is for imposing randomized heating at x-axis for 2.5D simulation, for details please see Zhou et al. (2020, Nature Astronomy,4,994) and
  !> Li et al. (2022,ApJ,926,216). 
  !> Users can define the parameters "trelax,vtim,vlim,periods,variation,num,mnx" in the mod_usr.t file.
  !> Example for using this module please see tests/mhd/coronal_rain_2.5D
  !trelax is the relaxization time for the simulation, if the relaxization is done seperately, then trelax=0.d0
  !periods=300,variation=75  !The time durations of each heating episode is typically 300s +/- 75s
  !vtim = 100                !100 heating episodes, e.g., 500 min. Enough for the simulation
  !num = 1000                !1000 waves for each episode
  !si=5./6.                  !spectral index
  !mnx could be chosen as twice of the spatial resolution in the simulation
  !vlim could be chosen as the same with mnx

contains

  subroutine getlQgrid(lQgrid,ixI^L,ixO^L,qt,trelax,x,mode,vtim,vlim)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, mode, vtim, vlim
    double precision, intent(in)    :: qt, trelax, x(ixI^S,1:ndim)
    double precision, intent(inout) :: lQgrid(ixI^S)

    double precision                :: tl1,tl2,tl3,tl4,tl5,tl6,t2,stt
    double precision                :: tr1,tr2,tr3,tr4,tr5,tr6
    double precision                :: va0
    integer                         :: Bp,Bt1,Bt2,ix^D,i

    t2 = qt - trelax
    va0=0.4d0
    stt = (xprobmax1- xprobmin1) / (dble(vlim) - 1.d0)

    select case(mode)
    ! mode 1 for gaussian, mode 2 for sinc
    case(1)
      if(qt .lt. trelax) then
        lQgrid(ixO^S)=0.d0
      else
        do i = 4, vtim
          if(qt .lt. trelax + tb1(i)) then
            tl1 = dexp(-(t2 - tb1(i-3))**2/(0.5d0*ta1(i-3))**2)
            tl2 = dexp(-(t2 - tb1(i-2))**2/(0.5d0*ta1(i-2))**2)
            tl3 = dexp(-(t2 - tb1(i-1))**2/(0.5d0*ta1(i-1))**2)
            tl4 = dexp(-(t2 - tb1(i))**2/(0.5d0*ta1(i))**2)
            tl5 = dexp(-(t2 - tb1(i+1))**2/(0.5d0*ta1(i+1))**2)
            tl6 = dexp(-(t2 - tb1(i+2))**2/(0.5d0*ta1(i+2))**2)
            {do ix^DB = ixOmin^DB,ixOmax^DB\}
              Bp = floor((x(ix^D,1) - xprobmin1) / stt) + 1
              if (Bp .lt. 1)  Bp=1
              if(x(ix^D,1) .lt. 0.d0) then
                lQgrid(ix^D) = tl1 * (va1(Bp + (i-4)*vlim)*(one-va0)+va0) &
                                + tl2 * (va1(Bp + (i-3)*vlim)*(one-va0)+va0) &
                                + tl3 * (va1(Bp + (i-2)*vlim)*(one-va0)+va0) &
                                + tl4 * (va1(Bp + (i-1)*vlim)*(one-va0)+va0) &
                                + tl5 * (va1(Bp + (i)*vlim)*(one-va0)+va0) &
                                + tl6 * (va1(Bp + (i+1)*vlim)*(one-va0)+va0)
              endif
            {enddo\}
            exit
          endif
        enddo

        do i = 4, vtim
          if(qt .lt. trelax + tb2(i)) then
            tr1 = dexp(-(t2 - tb2(i-3))**2/(0.5d0*ta2(i-3))**2)
            tr2 = dexp(-(t2 - tb2(i-2))**2/(0.5d0*ta2(i-2))**2)
            tr3 = dexp(-(t2 - tb2(i-1))**2/(0.5d0*ta2(i-1))**2)
            tr4 = dexp(-(t2 - tb2(i))**2/(0.5d0*ta2(i))**2)
            tr5 = dexp(-(t2 - tb2(i+1))**2/(0.5d0*ta2(i+1))**2)
            tr6 = dexp(-(t2 - tb2(i+2))**2/(0.5d0*ta2(i+2))**2)
            {do ix^DB = ixOmin^DB,ixOmax^DB\}
              Bp = floor((x(ix^D,1) - xprobmin1) / stt) + 1   
              if (Bp .lt. 1)  Bp=1
              if(x(ix^D,1) .ge. 0.d0) then
                lQgrid(ix^D) = tr1 * (va2(Bp + (i-4)*vlim)*(one-va0)+va0) &
                                + tr2 * (va2(Bp + (i-3)*vlim)*(one-va0)+va0) &
                                + tr3 * (va2(Bp + (i-2)*vlim)*(one-va0)+va0) &
                                + tr4 * (va2(Bp + (i-1)*vlim)*(one-va0)+va0) &
                                + tr5 * (va2(Bp + (i)*vlim)*(one-va0)+va0) &
                                + tr6 * (va2(Bp + (i+1)*vlim)*(one-va0)+va0)
              endif
            {enddo\}
            exit
          endif
        enddo
      endif
    case(2)
      if(qt .lt. trelax) then
        tr1 = 0.d0
        tr3 = 0.d0
        Bt1 = 0
      else if(qt .lt. trelax + tb1(1)) then
        tr1 = dsin(t2 * dpi / tb1(2))
        tr3 = 0.d0
        Bt1 = 0
      else if(qt .lt. trelax + tb1(2)) then
        tr1 = dsin(t2 * dpi / tb1(2))
        tr3 = dsin((t2 - tb1(1)) * dpi / (tb1(3) - tb1(1)))
        Bt1 = 0
      else
        do i = 2, vtim
          if(t2 .lt. tb1(i)) then
            Bt1 = i - 1
            exit
          endif
        enddo
        tr1 = dsin((t2 - tb1(Bt1-1)) * dpi / (tb1(Bt1+1) - tb1(Bt1-1)))
        tr3 = dsin((t2 - tb1(Bt1)) * dpi / (tb1(Bt1+2) - tb1(Bt1)))
      endif
      if(qt .lt. trelax) then
        tr2 = 0.d0
        tr4 = 0.d0
        Bt2 = 0
      else if(qt .lt. trelax + tb2(1)) then
        tr2 = dsin(t2 * dpi / tb2(2))
        tr4 = 0.d0
        Bt2 = 0
      else if(qt .lt. trelax + tb2(2)) then
        tr1 = dsin(t2 * dpi / tb2(2))
        tr3 = dsin((t2 - tb2(1)) * dpi / (tb2(3) - tb2(1)))
        Bt2 = 0
      else
        do i = 2, vtim
          if(t2 .lt. tb2(i)) then
            Bt2 = i - 1
            exit
          endif
        enddo
        tr2 = dsin((t2 - tb2(Bt2-1)) * dpi / (tb2(Bt2+1) - tb2(Bt2-1)))
        tr4 = dsin((t2 - tb2(Bt2)) * dpi / (tb2(Bt2+2) - tb2(Bt2)))
      endif
      {do ix^DB = ixOmin^DB,ixOmax^DB\}
          Bp = floor((x(ix^D,1) - xprobmin1) / stt) + 1
          if(x(ix^D,1) .lt. 0.d0) then
            lQgrid(ix^D) = tr1 * (va1(Bp + Bt1*vlim)*(one-va0)+va0) + tr3 * (va1(Bp+(Bt1+1)*vlim)*(one-va0)+va0)
          else
            lQgrid(ix^D) = tr2 * (va2(Bp + Bt2*vlim)*(one-va0)+va0) + tr4 * (va2(Bp+(Bt2+1)*vlim)*(one-va0)+va0)
          endif
      {enddo\}
    case default
      call mpistop('unknow mode')
    end select
  
  end subroutine getlQgrid

  subroutine generateTV(vtim,vlim,periods,variation,mnx,num,si)
    integer, intent(in)          :: vtim,vlim,mnx,num
    double precision, intent(in) :: periods,variation,si
    integer                      :: i,vx,j
    logical                      :: alive
    character(len=100)           :: filename

    ! generate random T, typically 300s +/- 75s   

    allocate(ta1(vtim))  ! left footpoint
    allocate(tb1(vtim)) 
    allocate(ta2(vtim))  ! right footpoint
    allocate(tb2(vtim))

    if (ndim/=2) call mpistop("Randomized heating for 2D, 2.5D simulations only!")
    if (mype==0) then 
      write(filename,"(a)") "LeftT.dat"
      inquire(file=filename,exist=alive)
      if(alive) then
        open(unit=21,file=filename,status='old')
        do j=1,vtim
          read(21,*) tb1(j)
        enddo
        close(21)
      else
        call randomT(tb1,vtim,periods,variation,filename)  
      endif 

      write(filename,"(a)") "RightT.dat"
      inquire(file=filename,exist=alive)
      if(alive) then
        open(unit=21,file=filename,status='old')
        do j=1,vtim
          read(21,*) tb2(j)
        enddo
        close(21)
      else
        call randomT(tb2,vtim,periods,variation,filename)  
      endif 

      ta1 = tb1
      ta2 = tb2
      do i = 2, vtim
        ta1(i) = tb1(i) - tb1(i-1)
        ta2(i) = tb2(i) - tb2(i-1)
      enddo
    endif   
    call MPI_BARRIER(icomm,ierrmpi)
    if (npe>1) then
      call MPI_BCAST(tb1,vtim,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(ta1,vtim,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(tb2,vtim,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(ta2,vtim,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    endif

    ! spatial distribution
  
    vx=vtim*vlim
    allocate(va1(vx)) 
    allocate(va2(vx))
    if (mype==0) then 
      write(filename,"(a)") "LeftV.dat"
      inquire(file=filename,exist=alive)
      if(alive) then
        open(unit=22,file=filename,status='old')
        do j=1,vx
          read(22,*) va1(j)
        enddo
        close(22)
      else     
        call randomV(va1,vtim,vlim,mnx,num,si,filename)
      endif
 
      write(filename,"(a)") "RightV.dat"
      inquire(file=filename,exist=alive)
      if(alive) then
        open(unit=22,file=filename,status='old')
        do j=1,vx
          read(22,*) va2(j)
        enddo
        close(22)
      else   
        call randomV(va2,vtim,vlim,mnx,num,si,filename)
      endif
    endif
    if (npe>1) then
      call MPI_BCAST(va1,vx,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(va2,vx,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    endif 

  end subroutine generateTV


  subroutine randomT(tb,vtim,periods,variation,filename)
    integer, intent(in)             :: vtim
    double precision, intent(in)    :: periods,variation
    double precision, intent(inout) :: tb(:)
    character(len=100), intent(in)  :: filename
    double precision                :: tt1,tt,mm
    integer                         :: i

    ! generate random T, typically 300s +/- 75s    
 
    tt = 0.d0
    do i = 1, vtim
      call random_number(mm)  
      tt1 = periods + variation * 2. * (mm - 0.5)
      tt = tt + tt1
      tb(i) = tt / unit_time
    end do
 
    open(unit=21,file=filename,form='formatted',status='new')
    write(21,'(es12.4)') tb
    close(21)
     
  end subroutine randomT

  subroutine randomV(va,vtim,vlim,mnx,num,si,filename)
    integer, intent(in)             :: vtim,vlim,mnx,num
    double precision, intent(in)    :: si
    double precision, intent(inout) :: va(:)
    character(len=100), intent(in)  :: filename
    double precision, allocatable   :: vn(:,:),vm(:)
    double precision                :: lambda,c,lambda_min,lambda_max,rms,ll
    integer                         :: i
    integer                         :: j, kk, lth = 1
 
    allocate(vn(num, mnx))
    allocate(vm(mnx))
    do kk = 1, vtim
        lambda_min = lth * 1.d0 / (mnx * 1.d0)
        lambda_max = lth * 2.d0
        do i = 1, num
          lambda = i * 1.d0 / num * (lambda_max - lambda_min) + lambda_min
          c = lambda ** si 
          call random_number(ll)        
          do j = 1, mnx  
            vn(i,j) = c*sin(2.*dpi*(j*1.d0/(mnx-1)*lth)/lambda+2.*dpi*2.*(ll-0.5))
          end do
        end do     
        rms = 0.d0
        vm = 0.d0    
        do j = 1, mnx  
          do i = 1, num
            vm(j) = vm(j) + vn(i, j) 
          end do
          rms = rms + vm(j)**2
        end do         
        vm = vm / dsqrt(rms/mnx)    ! normalization        
        do j = 1, vlim
          va(vlim*(kk-1)+j) = vm(j)**2
        end do       
    end do

    open(unit=22,file=filename,form='formatted',status='new')
    write(22,'(es12.4)') va
    close(22)

  end subroutine randomV

end module mod_random_heating
