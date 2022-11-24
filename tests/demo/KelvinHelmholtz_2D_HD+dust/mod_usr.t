!=============================================================================
! Simulate 2D KH evolution in 2-layer medium (Hillier 2019), use random vy
! setting for 2D hydro with dust species:
!=============================================================================
module mod_usr
  use mod_hd

  double precision:: min_ar,max_ar
  double precision:: delvy,delv
  double precision :: delrho=0.01d0
  logical:: reproduce
  double precision,allocatable, save:: xgrid(:),ygrid(:)
  integer, save:: Nx,Ny,ntmp
  double precision,allocatable, save:: rand(:,:)


contains

  subroutine usr_init()

    use mod_physics, only: phys_req_diagonal
    unit_length=1.0d18        ! in cm or 10^13 km (1/3 pc)
    unit_numberdensity= 1.0d3 ! number density per cubic cm
    unit_velocity= 1.0d6      ! in cm/s, or 10km/s

    call usr_params_read(par_files)

    usr_init_one_grid   => initonegrid_usr
    usr_set_parameters  => initglobaldata_usr
    usr_aux_output      => extra_var_output
    usr_add_aux_names   => extra_var_names_output
    usr_var_for_errest  => myvar_for_errest
    usr_refine_grid     => specialrefine_grid

    call set_coordinate_system("Cartesian_2D")
    call hd_activate()
    phys_req_diagonal = .true.
  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_dust

   integer,dimension(:),allocatable:: seed
   integer::  seed_size,ix
   double precision :: Lx, Ly

   double precision :: r(0:dust_n_species)
   double precision :: dsdust(dust_n_species)
   integer          :: i
   logical, save    :: first = .true.

   hd_gamma=5.0d0/3.0d0

   Ly=xprobmax2-xprobmin2
   Lx=xprobmax1-xprobmin1
   Nx=nint(Lx/dx(1,1))
   if (.not. allocated(xgrid)) allocate(xgrid(Nx))
   do ix=1,Nx
      xgrid(ix)=xprobmin1+ix*dx(1,1)
   enddo
   Ny=nint(Ly/dx(2,1))
   if (.not. allocated(ygrid)) allocate(ygrid(Ny))
   do ix=1,Ny
      ygrid(ix)=xprobmin2+ix*dx(2,1)
   enddo
   ntmp=Nx*Ny
   if (.not. allocated(rand)) allocate(rand(Nx,Ny))

   if(mype==0)then
      print *,'Number of modes (x,y)=',Nx,Ny
      call random_seed(SIZE=seed_size)
      allocate(seed(seed_size))
      if(reproduce)then
         seed = 20180815
        call random_seed(PUT=seed(1:seed_size))
      else
        call random_seed(GET=seed(1:seed_size))
      endif
      call random_number(rand) ! generate ntmp random numbers [0,1]
   endif
   call MPI_BARRIER(icomm,ierrmpi)
   if(npe>1) call MPI_BCAST(rand,ntmp,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)

   if (dust_n_species > 0) then
      ! specific dust grain density
      dust_density(1:dust_n_species) = 3.3d0  ! in g/cm   
      ! range in dust grain sizes
      min_ar              =   5.0d-7   ! minimum dust grain size hence 5 nm
      max_ar              = 250.0d-7   ! maximum dust grain size hence 250 nm
      ! === rescale dust quantities to dimensionless scale === !
      dust_density(1:dust_n_species) = dust_density(1:dust_n_species)/unit_density
      min_ar  = min_ar/unit_length
      max_ar  = max_ar/unit_length
      !-------------------------------

      ! here the dust sizes are defined. Ndust bins, with all bins having equal total mass.
      ! To do this, assume the particle distribution goes as r^-3.5
      r(0) = min_ar
      do i=1,dust_n_species
        r(i) = (sqrt(r(i-1)) +(sqrt(max_ar)- &
             sqrt(min_ar))/dust_n_species)**2.0d0
        dsdust(i) = r(i)-r(i-1)
      end do

      ! now calculate the weighted mean size of each bin, again assuming n goes as r^-3.5
      ! this assumes we weigh with r^2
      do i=1,dust_n_species
         dust_size(i) = dsqrt(5.0d0*(r(i)**(-0.5d0) - r(i-1)**(-0.5d0)) &
                                   /(r(i)**(-2.5d0) - r(i-1)**(-2.5d0)))
      end do

      if(first)then
         if(mype==0)then
            write(*,*) '*****************************************'
            if(SI_unit)then
               write(*,*) 'Units system in SI'
            else
               write(*,*) 'Units system in cgs'
            endif
            write(*,*) 'He_abundance is       =',He_abundance
            write(*,*) 'unit length is        =',unit_length
            write(*,*) 'unit number density is=',unit_numberdensity
            write(*,*) 'unit velocity is      =',unit_velocity
            write(*,*) 'unit time is          =',unit_time
            write(*,*) 'unit density is       =',unit_density
            write(*,*) 'unit pressure is      =',unit_pressure
            write(*,*) 'unit temperature is   =',unit_temperature
            write(*,*) 'specific heat ratio is=',hd_gamma
            write(*,*) '*****************************************'
            write(*,*) 'Dust included using ',dust_n_species,' dust species'
            write(*,*) 'Dust bins all have specific density rhop ',dust_density(1)
            write(*,*) '   in cgs units specific density is rhop ',dust_density(1)*unit_density
            write(*,*) 'Dust bins between min=',min_ar,' and max=',max_ar
            write(*,*) '   in cgs units from ',min_ar*unit_length,' to ',max_ar*unit_length
            do i=1,dust_n_species
              write(*,*) 'Dust type ',i,': grain radius r              =', dust_size(i)*unit_length
              write(*,*) 'Dust type ',i,': dimensionless grain radius r=', dust_size(i)
              write(*,*) 'Dust type ',i,': dimensionless rhop x r      =', dust_size(i)*dust_density(i)
            end do
            write(*,*) '*****************************************'
         endif
         first=.false.
      endif
   end if

 end subroutine initglobaldata_usr

 subroutine usr_params_read(files)
    use mod_global_parameters, only:unitpar
    use mod_constants
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ delvy,delv,delrho,reproduce
    do n=1, size(files)
      open(unitpar, file=trim(files(n)), status='old')
      read(unitpar, usr_list, end=111)
      111 close(unitpar)
    end do

 end subroutine usr_params_read

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    use mod_dust
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    integer:: i,ix^D,n
    double precision :: vely(ixI^S)
    double precision :: rand_vec(1:Ny)
    double precision :: tmp,rand_avg,rand_std
    logical, save :: first=.true.

    if(first)then
       if(mype==0) then
         write(*,*)'Simulate 2D KH with random vy with amplitudes=',delvy,delv
       endif
       first=.false.
    endif
    rand_avg=SUM(rand)/DBLE(ntmp)
    rand_std=SQRT(SUM((rand-rand_avg)**2)/DBLE(Ntmp-1))
    vely(ixO^S)=zero
    do ix1=ixOmin1,ixOmax1
       do ix2=ixOmin2,ixOmax2
           do i=1,Ny ! interpolate at x1 along xgrid for each y-position
              call interp1d(Nx,xgrid,rand(1:Nx,i),x(ix1,ix2,1),tmp)
              rand_vec(i)=tmp
           enddo
           call interp1d(Ny,ygrid,rand_vec(1:Ny),x(ix1,ix2,2),tmp)
           ! scale to delrho above average value
           vely(ix1,ix2)=vely(ix1,ix2)+delvy*(tmp-rand_avg)/rand_std
       enddo
    enddo

    ! reference density at inlet
    where (x(ixO^S,2)>0.0d0) 
        w(ixO^S,rho_)=1.5d0
        w(ixO^S,mom(1))=delv/2.5d0
    elsewhere
        w(ixO^S,rho_)=1.0d0
        w(ixO^S,mom(1))=-delv*1.5d0/2.5d0
    endwhere
    ! velocity perturbation
    w(ixO^S,mom(2))=vely(ixO^S)
    ! pressure to make bottom sound speed unity
    w(ixO^S,e_)=1.0d0/hd_gamma

   if (dust_n_species > 0) then
    do n = 1, dust_n_species
      w(ixO^S, dust_mom(:, n)) = 0.0d0

      where(x(ixO^S,2)>0.0d0)
        w(ixO^S,dust_rho(n))    = delrho/dust_n_species
        w(ixO^S,dust_mom(1, n)) = delv/2.5d0
      elsewhere
        !w(ixO^S,dust_rho(n)) = 0.0d0
        w(ixO^S,dust_rho(n)) =  0.1*delrho/dust_n_species
        w(ixO^S,dust_mom(1,n))=-delv*1.5d0/2.5d0
      endwhere
    end do
   endif

    call hd_to_conserved(ixI^L, ixO^L, w, x)

  end subroutine initonegrid_usr


  subroutine extra_var_output(ixI^L, ixO^L, w, x, normconv)
    use mod_physics
    integer, intent(in)              :: ixI^L,ixO^L
    double precision, intent(in)     :: x(ixI^S,1:ndim)
    double precision                 :: w(ixI^S,nw+nwauxio)
    double precision                 :: normconv(0:nw+nwauxio)
    double precision                 :: pp(ixI^S),vv(ixI^S),cs2(ixI^S)
    double precision                 :: wlocal(ixI^S,1:nw)
    double precision :: rho(ixI^S),gradrho(ixI^S),drho(ixI^S)
    double precision :: kk,kk0,grhomax,kk1
    integer :: idims

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    call hd_get_pthermal(wlocal,x,ixI^L,ixO^L,pp)
    w(ixO^S,nw+1) = pp(ixO^S)/w(ixO^S,rho_)
    vv(ixO^S)=w(ixO^S,mom(1))/w(ixO^S,rho_)
    call hd_get_csound2(w,x,ixI^L,ixO^L,cs2)
    w(ixO^S,nw+2) =vv(ixO^S)/dsqrt(cs2(ixO^S))
    ! then compute schlieren density plot
    rho(ixI^S)=w(ixI^S,rho_)
    gradrho(ixO^S)=zero
    do idims=1,ndim
      select case(typegrad)
        case("central")
          call gradient(rho,ixI^L,ixO^L,idims,drho)
        case("limited")
          call gradientS(rho,ixI^L,ixO^L,idims,drho)
      end select
      gradrho(ixO^S)=gradrho(ixO^S)+drho(ixO^S)**2.0d0
    enddo
    gradrho(ixO^S)=dsqrt(gradrho(ixO^S))
    kk=5.0d0
    kk0=0.001d0
    kk1=1.0d0
    grhomax=2000.0d0
    w(ixO^S,nw+3)=dexp(-kk*(gradrho(ixO^S)-kk0*grhomax)/(kk1*grhomax-kk0*grhomax))
    w(ixO^S,nw+4)=dlog10(w(ixO^S,rho_))

  end subroutine extra_var_output

  subroutine extra_var_names_output(varnames)
    character(len=*)  :: varnames

    varnames = "T Mx schlier logrho"

  end subroutine extra_var_names_output

subroutine interp1d(nin,xin,yin,xout1,yout1)
! subroutine to interpolate values from a 1D function
  integer,intent(in) :: nin
  DOUBLE PRECISION,intent(in) :: xin(1:nin),yin(1:nin)
  DOUBLE PRECISION,intent(in) :: xout1
  DOUBLE PRECISION,intent(out) :: yout1
  integer*4 :: ixa
  DOUBLE PRECISION :: xa1,xb1,ya1,yb1,xtmp

  ixa=minloc(abs(xin-xout1),1)
  xtmp=xin(ixa)
  ! not in range of xin[0] and xin[nin]
  ! use extraplotion
  if (ixa == 1) then
     xa1=xin(1)
     xb1=xin(2)
     ya1=yin(1)
     yb1=yin(2)
     yout1=ya1+(yb1-ya1)*(xout1-xa1)/(xb1-xa1)
  endif
  if (ixa == nin) then
     xa1=xin(ixa-1)
     xb1=xin(ixa)
     ya1=yin(ixa-1)
     yb1=yin(ixa)
     yout1=ya1+(yb1-ya1)*(xout1-xa1)/(xb1-xa1)
  endif

  if (ixa < nin .and. ixa > 1) then
    if (xtmp > xout1) then
      ixa=ixa-1
    endif
    xa1=xin(ixa)
    xb1=xin(ixa+1)
    ya1=yin(ixa)
    yb1=yin(ixa+1)
    yout1=ya1+(yb1-ya1)*(xout1-xa1)/(xb1-xa1)
  endif

end subroutine interp1d

  subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    ! you must set consistent values for integers refine/coarsen:
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen
    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    if(qt<10.d0.and.any(dabs(x(ix^S,2))<= 0.02d0) ) then
      coarsen = -1
      refine  = 1
    endif

  end subroutine specialrefine_grid

  subroutine myvar_for_errest(ixI^L,ixO^L,iflag,w,x,var)
      use mod_global_parameters
      integer, intent(in)           :: ixI^L,ixO^L,iflag
      double precision, intent(in)  :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
      double precision, intent(out) :: var(ixI^S)

      if (iflag >nw+1)call mpistop(' iflag error')
      call hd_get_pthermal(w,x,ixI^L,ixO^L,var)
      var(ixO^S)=var(ixO^S)/w(ixO^S,rho_)

  end subroutine myvar_for_errest

end module mod_usr
