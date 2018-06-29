
!! Copyright (C) 2002-2007 BigDFT group 
!! This file is distributed under the terms of the
!! GNU General Public License, see ~/COPYING file
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the list of contributors, see ~/AUTHORS 


!!****p* BigDFT/Poisson_Solver
!! NAME
!!   Poisson_Solver
!!
!! FUNCTION
!!    Program test for Poisson
!!    Laplacian V = 4pi rho
!!    May work either in parallel or in serial case
!!    And for different geometries
!!
!! AUTHOR
!!    Luigi Genovese
!!
!! COPYRIGHT
!!    Copyright (C) 2007 CEA
!! CREATION DATE
!!    February 2007
!!
!! SOURCE
!!
program PoissonSolver

  use poisson_solver

  use mpi
  implicit none
  !Order of interpolating scaling function
  !integer, parameter :: itype_scf=8
  real(kind=8), parameter :: a_gauss = 1.0d0,a2 = a_gauss**2
  !Length of the box
  real(kind=8), parameter :: acell = 10.d0
  character(len=50) :: chain
  character(len=2) :: geocode
  character(len=1) :: datacode
  real(kind=8), dimension(:,:,:), allocatable :: density, rhopot,potential,pot_ion
  real(kind=8), pointer :: karray(:)
  real(kind=8) :: hx,hy,hz,max_diff,eh,exc,vxc,hgrid,diff_parser,offset
  real(kind=8) :: ehartree,eexcu,vexcu,diff_par,diff_ser
  integer :: n01,n02,n03,n1,n2,n3,itype_scf,i_all,i_stat
  integer :: i1,i2,i3,j1,j2,j3,i1_max,i2_max,i3_max,iproc,nproc,ierr,i3sd,ncomp
  integer :: n_cell,ixc,n3d,n3p,n3pi,i3xcsh,i3s
  logical :: alsoserial,onlykernel
  integer :: n_args

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  n_args = command_argument_count()

  if (n_args < 3 .or. n_args > 6) then
     if (iproc == 0) then
        print *, "Arguments: nx ny nz [ixc geocode datacode]"
        print *, "nx, ny, nz: number of grid points"
        print *, "ixc: control exchange correlation functional"
        print *, "     default: ixc = 0 (disabled)"
        print *, "geocode:  P - periodic, F - free, S - surface)"
        print *, "     default: F (free)"
        print *, "datacode:  G - global, D - distributed"
        print *, "     default: G (global)"
     end if
     call MPI_FINALIZE(ierr)
     stop
  end if

  !Use arguments
  call getarg(1,chain)
  read(unit=chain,fmt=*) n01
  call getarg(2,chain)
  read(unit=chain,fmt=*) n02
  call getarg(3,chain)
  read(unit=chain,fmt=*) n03

  if (n_args >= 4) then
     call getarg(4,chain)
     read(unit=chain,fmt=*) ixc
  else
     ixc = 0
  end if

  if (n_args >= 5) then
     call getarg(5,chain)
     read(unit=chain,fmt=*) geocode
  else
     geocode = 'F'
  end if

  if (n_args >= 6) then
     call getarg(6,chain)
     read(unit=chain,fmt=*) datacode
  else
     datacode = 'G'
  end if

  !perform also the comparison wit the serial case
  alsoserial=.false.
  onlykernel=.false.
  !code for the Poisson Solver in the parallel case
  !datacode='G'

  if (geocode == 'P') then

     if (iproc==0) print *,"PSolver, periodic BC: ",n01,n02,n03,'processes',nproc

  else if (geocode == 'S') then

     if (iproc==0) print *,"PSolver for surfaces: ",n01,n02,n03,'processes',nproc

  else if (geocode == 'F') then

     if (iproc==0) print *,"PSolver, free BC: ",n01,n02,n03,'processes',nproc

  end if

  !initialize memory counting
  call memocc(0,iproc,'count','start')

  !Step size
  n_cell = max(n01,n02,n03)
  hx=acell/real(n01,kind=8)
  hy=acell/real(n02,kind=8)
  hz=acell/real(n03,kind=8)

  !grid for the free BC case
  hgrid=max(hx,hy,hz)
  !hgrid=hx

  !we must choose properly a test case with a positive density
  itype_scf=14

  call timing(iproc,'parallel      ','IN')

  call createKernel(geocode,n01,n02,n03,hx,hy,hz,itype_scf,iproc,nproc,karray)

  if (.not. onlykernel) then
     !Allocations
     !Density
     allocate(density(n01,n02,n03),stat=i_stat)
     call memocc(i_stat,product(shape(density))*kind(density),'density','poisson_solver')
     !Density then potential
     allocate(rhopot(n01,n02,n03),stat=i_stat)
     call memocc(i_stat,product(shape(rhopot))*kind(rhopot),'rhopot','poisson_solver')
     allocate(potential(n01,n02,n03),stat=i_stat)
     call memocc(i_stat,product(shape(potential))*kind(potential),'potential','poisson_solver')
     !ionic potential
     allocate(pot_ion(n01,n02,n03),stat=i_stat)
     call memocc(i_stat,product(shape(pot_ion))*kind(pot_ion),'pot_ion','poisson_solver')

     call test_functions(geocode,ixc,n01,n02,n03,acell,a_gauss,hx,hy,hz,&
          density,potential,rhopot,pot_ion)


     !offset, used only for the periodic solver case
     if (ixc==0) offset=potential(1,1,1)!-pot_ion(1,1,1)

     !dimension needed for allocations
     call PS_dim4allocation(geocode,datacode,iproc,nproc,n01,n02,n03,ixc,n3d,n3p,n3pi,i3xcsh,i3s)

     !dimension for comparison in the global or distributed poisson solver
     if (datacode == 'G') then
        i3sd=1
        ncomp=n03
     else if (datacode == 'D') then
        i3sd=i3s
        ncomp=n3p
     end if

!!$  print *,'iproc,i3xcsh,i3s',iproc,i3xcsh,i3s

     !apply the Poisson Solver (case with distributed potential
     call PSolver(geocode,datacode,iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
          density(1,1,i3sd),karray,pot_ion(1,1,i3s+i3xcsh),ehartree,eexcu,vexcu,offset,.true.,1)

  end if

  i_all=-product(shape(karray))*kind(karray)
  deallocate(karray,stat=i_stat)
  call memocc(i_stat,i_all,'karray','poisson_solver')

  call timing(iproc,'              ','RE')

  if (.not. onlykernel) then

     !comparison (each process compare its own part)
     call compare(n01,n02,ncomp,potential(1,1,i3sd+i3xcsh),density(1,1,i3sd+i3xcsh),&
          i1_max,i2_max,i3_max,max_diff)

!!$  print *,'iproc,i3xcsh,i3s,max_diff',iproc,i3xcsh,i3s,max_diff

     !extract the max
     if (nproc > 1) then
        call MPI_ALLREDUCE(max_diff,diff_par,1,MPI_double_precision,  &
             MPI_MAX,MPI_COMM_WORLD,ierr)
     else
        diff_par=max_diff
     end if


     if (iproc == 0) then

        write(*,*) '--------------------'
        write(*,*) 'Parallel calculation '
        write(unit=*,fmt="(1x,a,3(1pe20.12))") "eht, exc, vxc:",ehartree,eexcu,vexcu
        write(*,'(a,3(i0,1x))') '  Max diff at: ',i1_max,i2_max,i3_max
        write(unit=*,fmt="(1x,a,1pe20.12)") '    Max diff:',diff_par,&
             '      result:',density(i1_max,i2_max,i3_max),&
             '    original:',potential(i1_max,i2_max,i3_max)
     end if

  end if

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !Serial case
  if (alsoserial) then
     call timing(0,'             ','IN')

     call createKernel(geocode,n01,n02,n03,hx,hy,hz,itype_scf,0,1,karray)

     if (.not. onlykernel) then
        !offset, used only for the periodic solver case
        if(ixc==0) offset=potential(1,1,1)!-pot_ion(1,1,1)

        !apply the Poisson Solver (case with distributed potential
        call PSolver(geocode,'G',0,1,n01,n02,n03,ixc,hx,hy,hz,&
             rhopot,karray,pot_ion,eh,exc,vxc,offset,.true.,1)

     end if

     i_all=-product(shape(karray))*kind(karray)
     deallocate(karray,stat=i_stat)
     call memocc(i_stat,i_all,'karray','poisson_solver')

     call timing(iproc,'              ','RE')

     if (.not. onlykernel) then
        !Maximum difference
        call compare(n01,n02,n03,potential,rhopot,i1_max,i2_max,i3_max,diff_ser)

!!$     print *,'iproc,diff_ser',iproc,diff_ser

        if (iproc==0) then
           write(*,*) '------------------'
           write(*,*) 'Serial Calculation'
           write(*,"(1x,a,3(1pe20.12))") "eht, exc, vxc:",eh,exc,vxc
           write(*,'(a,3(i0,1x))') '  Max diff at: ',i1_max,i2_max,i3_max
           write(*,"(1x,a,1pe20.12)") '     Max diff:',diff_ser,&
                '       result:',rhopot(i1_max,i2_max,i3_max),&
                '     original:',potential(i1_max,i2_max,i3_max)
        end if

        !Maximum difference, parallel-serial
        call compare(n01,n02,ncomp,rhopot(1,1,i3sd+i3xcsh),density(1,1,i3sd+i3xcsh),&
             i1_max,i2_max,i3_max,max_diff)

!!$     print *,'max_diff,i1_max,i2_max,i3_max,i3s,i3xcsh,n3p',max_diff,i1_max,i2_max,i3_max,&
!!$          i3s,i3xcsh,n3p

        if (nproc > 1) then
           !extract the max
           call MPI_ALLREDUCE(max_diff,diff_parser,1,MPI_double_precision,  &
                MPI_MAX,MPI_COMM_WORLD,ierr)
        else
           diff_parser=max_diff
        end if

        if (iproc==0) then
           write(*,*) '------------------'
           write(*,'(a,3(i0,1x))')&
                'difference parallel-serial, at',i1_max,i2_max,i3_max
           write(*,"(1x,a,1pe12.4)")&
                '    Max diff:',diff_parser,&
                '    parallel:',density(i1_max,i2_max,i3_max),&
                '      serial:',rhopot(i1_max,i2_max,i3_max)
           write(*,"(1x,a,3(1pe12.4))")&
                "energy_diffs:",ehartree-eh,eexcu-exc,vexcu-vxc
        end if
     end if
  end if
  if (iproc==0 .and. .not. onlykernel) then

     call regroup_data(geocode,n01,n02,n03,hx,hy,hz,diff_par,diff_parser)

     i2=i2_max
     do i3=1,n03
        do i1=1,n01
           j1=n1/2+1-abs(n1/2+1-i1)
           j2=n2/2+1-abs(n2/2+1-i2)
           j3=n3/2+1-abs(n3/2+1-i3)
           write(11,*)i1,i3,rhopot(i1,i2,i3),potential(i1,i2,i3),&
                density(i1,i2,i3)
        end do
     end do
     i3=i3_max
     do i2=1,n02
        do i1=1,n01
           j1=n1/2+1-abs(n1/2+1-i1)
           j2=n2/2+1-abs(n2/2+1-i2)
           j3=n3/2+1-abs(n3/2+1-i3)
           write(12,*)i1,i2,rhopot(i1,i2,i3),potential(i1,i2,i3),&
                density(i1,i2,i3)
        end do
     end do

  end if

  if (.not. onlykernel) then
     i_all=-product(shape(density))*kind(density)
     deallocate(density,stat=i_stat)
     call memocc(i_stat,i_all,'density','poisson_solver')
     i_all=-product(shape(rhopot))*kind(rhopot)
     deallocate(rhopot,stat=i_stat)
     call memocc(i_stat,i_all,'rhopot','poisson_solver')
     i_all=-product(shape(potential))*kind(potential)
     deallocate(potential,stat=i_stat)
     call memocc(i_stat,i_all,'potential','poisson_solver')
     i_all=-product(shape(pot_ion))*kind(pot_ion)
     deallocate(pot_ion,stat=i_stat)
     call memocc(i_stat,i_all,'pot_ion','poisson_solver')
  end if

  !finalize memory counting
  call memocc(0,0,'count','stop')

  call MPI_FINALIZE(ierr)  

end program PoissonSolver

subroutine regroup_data(geocode,n01,n02,n03,hx,hy,hz,max_diff,diff_parser)
  implicit none
  character(len=2), intent(in) :: geocode
  integer, intent(in) ::n01,n02,n03
  real(kind=8), intent(in) :: hx,hy,hz,max_diff,diff_parser
  !local variables
  character(len=14) :: string
  real(kind=8) :: tcp1,tcp2,tcm1,tcm2,tk1,tk2,txc1,txc2,pcp,pcm,pk,pxc
  real(kind=8) :: tcp,tcm,tk,txc,hgrid
  hgrid=max(hx,hy,hz)
  if (geocode == 'S') hgrid=hy

  open(unit=60,file='time.par',status='unknown')
  read(60,*)
  read(60,*)string,tcp1,tcp2,pcp
  read(60,*)string,tcm1,tcm2,pcm
  read(60,*)string,tk1,tk2,pk
  read(60,*)string,txc1,txc2,pxc
  close(60)
  tcp=tcp2
  tcm=tcm2
  tk=tk2
  txc=txc2
  
  write(99,'(a2,3(i4),1pe9.2,1pe10.3,4(1pe9.2),4(0pf5.1),1pe9.2)')&
       geocode,n01,n02,n03,hgrid,max_diff,tcp,tcm,tk,txc,pcp,pcm,pk,pxc,diff_parser
  
end subroutine regroup_data

subroutine compare(n01,n02,n03,potential,density,i1_max,i2_max,i3_max,max_diff)
  implicit none
  integer, intent(in) :: n01,n02,n03
  real(kind=8), dimension(n01,n02,n03), intent(in) :: potential,density
  integer, intent(out) :: i1_max,i2_max,i3_max
  real(kind=8), intent(out) :: max_diff

  !local variables
  integer :: i1,i2,i3
  real(kind=8) :: factor
  max_diff = 0.d0
  i1_max = 1
  i2_max = 1
  i3_max = 1
  do i3=1,n03
     do i2=1,n02 
        do i1=1,n01
           factor=abs(potential(i1,i2,i3)-density(i1,i2,i3))
           if (max_diff < factor) then
              max_diff = factor
              i1_max = i1
              i2_max = i2
              i3_max = i3
           end if
        end do
     end do
  end do
end subroutine compare


! this subroutine builds some analytic functions that can be used for 
! testing the poisson solver.
! The default choice is already well-tuned for comparison.
! WARNING: not all the test functions can be used for all the boundary conditions of
! the poisson solver, in order to have a reliable analytic comparison.
! The parameters of the functions must be adjusted in order to have a sufficiently localized
! function in the isolated direction and an explicitly periodic function in the periodic ones.
! Beware of the high-frequency components that may falsify the results when hgrid is too high.
subroutine test_functions(geocode,ixc,n01,n02,n03,acell,a_gauss,hx,hy,hz,&
     density,potential,rhopot,pot_ion)
  implicit none
  character(len=2), intent(in) :: geocode
  integer, intent(in) :: n01,n02,n03,ixc
  real(kind=8), intent(in) :: acell,a_gauss,hx,hy,hz
  real(kind=8), dimension(n01,n02,n03), intent(out) :: density,potential,rhopot,pot_ion

  !local variables
  integer :: i1,i2,i3,nu,ifx,ify,ifz
  real(kind=8) :: x,x1,x2,x3,y,length,denval,pi,a2,derf,factor,r,r2
  real(kind=8) :: fx,fx2,fy,fy2,fz,fz2,a,ax,ay,az,bx,by,bz,tt,potion_fac

  if (trim(geocode) == 'P') then
     !parameters for the test functions
     length=acell
     a=0.5d0/a_gauss**2
     !test functions in the three directions
     ifx=5
     ify=5
     ifz=5
     !parameters of the test functions
     ax=length
     ay=length
     az=length
     bx=2.d0!real(nu,kind=8)
     by=2.d0!real(nu,kind=8)
     bz=2.d0

!!$     !plot of the functions used
!!$     do i1=1,n03
!!$        x = hx*real(i1,kind=8)!valid if hy=hz
!!$        y = hz*real(i1,kind=8) 
!!$        call functions(x,ax,bx,fx,fx2,ifx)
!!$        call functions(y,az,bz,fz,fz2,ifz)
!!$        write(20,*)i1,fx,fx2,fz,fz2
!!$     end do

     !Initialization of density and potential
     denval=0.d0 !value for keeping the density positive
     do i3=1,n03
        x3 = hz*real(i3-n03/2-1,kind=8)
        call functions(x3,az,bz,fz,fz2,ifz)
        do i2=1,n02
           x2 = hy*real(i2-n02/2-1,kind=8)
           call functions(x2,ay,by,fy,fy2,ify)
           do i1=1,n01
              x1 = hx*real(i1-n01/2-1,kind=8)
              call functions(x1,ax,bx,fx,fx2,ifx)
              density(i1,i2,i3) = fx2*fy*fz+fx*fy2*fz+fx*fy*fz2
              denval=max(denval,-density(i1,i2,i3))
              potential(i1,i2,i3) = fx*fy*fz!density(i1,i2,i3)
           end do
        end do
     end do


!plane capacitor oriented along the y direction
!!$     do i2=1,n02
!!$        if (i2==n02/4) then
!!$           do i3=1,n03
!!$              do i1=1,n01
!!$                 density(i1,i2,i3)=1.d0!real(i2,kind=8)
!!$              end do
!!$           end do
!!$        else if (i2==3*n02/4) then
!!$           do i3=1,n03
!!$              do i1=1,n01
!!$                 density(i1,i2,i3)=-1.d0!real(i2,kind=8)
!!$              end do
!!$           end do
!!$        else
!!$           do i3=1,n03
!!$              do i1=1,n01
!!$                 density(i1,i2,i3)=0.d0
!!$              end do
!!$           end do
!!$        end if
!!$     end do
!!$     denval=0.d0

     if (ixc==0) denval=0.d0

  else if (trim(geocode) == 'S') then
     !parameters for the test functions
     length=acell
     a=0.5d0/a_gauss**2
     !test functions in the three directions
     ifx=5
     ifz=5
     !non-periodic dimension
     ify=6
     !parameters of the test functions
     ax=length
     az=length
     bx=real(nu,kind=8)
     bz=real(nu,kind=8)
     !non-periodic dimension
     ay=length!4.d0*a
     by=a
     density(:,:,:) = 0.d0!1d-20 !added

     !plot of the functions used
     do i1=1,n02
        x = hx*real(i1-n02/2-1,kind=8)!valid if hy=hz
        y = hy*real(i1-n02/2-1,kind=8) 
        call functions(x,ax,bx,fx,fx2,ifx)
        call functions(y,ay,by,fy,fy2,ify)
        write(20,*)i1,fx,fx2,fy,fy2
     end do

     !Initialisation of density and potential
     !Normalisation
     do i3=1,n03
        x3 = hz*real(i3-n03/2-1,kind=8)
        call functions(x3,az,bz,fz,fz2,ifz)
        do i2=1,n02
           x2 = hy*real(i2-n02/2-1,kind=8)
           call functions(x2,ay,by,fy,fy2,ify)
           do i1=1,n01
              x1 = hx*real(i1-n02/2-1,kind=8)
              call functions(x1,ax,bx,fx,fx2,ifx)
              density(i1,i2,i3) = fx2*fy*fz+fx*fy2*fz+fx*fy*fz2
              denval=max(denval,-density(i1,i2,i3))
              potential(i1,i2,i3) = fx*fy*fz
           end do
        end do
     end do

!plane capacitor oriented along the y direction
!!$     do i2=1,n02
!!$        if (i2==1) then
!!$           do i3=1,n03
!!$              do i1=1,n01
!!$                 density(i1,i2,i3)=1.d0!real(i2,kind=8)
!!$              end do
!!$           end do
!!$        else if (i2==n02) then
!!$           do i3=1,n03
!!$              do i1=1,n01
!!$                 density(i1,i2,i3)=-1.d0!real(i2,kind=8)
!!$              end do
!!$           end do
!!$        else
!!$           do i3=1,n03
!!$              do i1=1,n01
!!$                 density(i1,i2,i3)=0.d0
!!$              end do
!!$           end do
!!$        end if
!!$     end do

     if (ixc==0) denval=0.d0

  else if (trim(geocode) == 'F') then

     !grid for the free BC case
     !hgrid=max(hx,hy,hz)

     pi = 4.d0*atan(1.d0)
     a2 = a_gauss**2

     !Normalization
     factor = 1.d0/(a_gauss*a2*pi*sqrt(pi))
     !gaussian function
     do i3=1,n03
        x3 = hz*real(i3-n03/2,kind=8)
        do i2=1,n02
           x2 = hy*real(i2-n02/2,kind=8)
           do i1=1,n01
              x1 = hx*real(i1-n01/2,kind=8)
              r2 = x1*x1+x2*x2+x3*x3
              density(i1,i2,i3) = factor*exp(-r2/a2)
              r = sqrt(r2)
              !Potential from a gaussian
              if (r == 0.d0) then
                 potential(i1,i2,i3) = 2.d0/(sqrt(pi)*a_gauss)
              else
                 potential(i1,i2,i3) = derf(r/a_gauss)/r
              end if
           end do
        end do
     end do

!plane capacitor oriented along the y direction
!!$     do i2=1,n02
!!$        if (i2==n02/4) then
!!$           do i3=1,n03
!!$              do i1=1,n01
!!$                 density(i1,i2,i3)=1.d0!real(i2,kind=8)
!!$              end do
!!$           end do
!!$        else if (i2==3*n02/4) then
!!$           do i3=1,n03
!!$              do i1=1,n01
!!$                 density(i1,i2,i3)=-1.d0!real(i2,kind=8)
!!$              end do
!!$           end do
!!$        else
!!$           do i3=1,n03
!!$              do i1=1,n01
!!$                 density(i1,i2,i3)=0.d0
!!$              end do
!!$           end do
!!$        end if
!!$     end do
     
     denval=0.d0

  else

     print *,'geometry code not admitted',geocode
     stop

  end if

! For ixc/=0 the XC potential is added to the solution, and an analytic comparison is no more
! possible. In that case the only possible comparison is between the serial and the parallel case
! To ease the comparison between the serial and the parallel case we add a random pot_ion
! to the potential.

  if (ixc==0) then
     potion_fac=0.d0
  else
     potion_fac=1.d0
  end if

     rhopot(:,:,:) = density(:,:,:) + denval
     do i3=1,n03
        do i2=1,n02
           do i1=1,n01
              call random_number(tt)
              !tt=0.d0!1.d0
              pot_ion(i1,i2,i3)=tt
              potential(i1,i2,i3)=potential(i1,i2,i3)+potion_fac*tt
!!$              !for the ixc/=0 case
!!$              call random_number(tt)
!!$              rhopot(i1,i2,i3)=abs(tt)
           end do
        end do
     end do
     if (denval /= 0.d0) density=rhopot

end subroutine test_functions

subroutine functions(x,a,b,f,f2,whichone)
  implicit none
  integer, intent(in) :: whichone
  real(kind=8), intent(in) :: x,a,b
  real(kind=8), intent(out) :: f,f2
  !local variables
  real(kind=8) :: r,r2,y,yp,ys,factor,pi,g,h,g1,g2,h1,h2
  real(kind=8) :: length,frequency,nu,sigma,agauss

  pi = 4.d0*datan(1.d0)
  select case(whichone)
  case(1)
     !constant
     f=1.d0
     f2=0.d0
  case(2)
     !gaussian of sigma s.t. a=1/(2*sigma^2)
     r2=a*x**2
     f=dexp(-r2)
     f2=(-2.d0*a+4.d0*a*r2)*dexp(-r2)
  case(3)
     !gaussian "shrinked" with a=length of the system
     length=a
     r=pi*x/length
     y=dtan(r)
     yp=pi/length*1.d0/(dcos(r))**2
     ys=2.d0*pi/length*y*yp
     factor=-2.d0*ys*y-2.d0*yp**2+4.d0*yp**2*y**2
     f2=factor*dexp(-y**2)
     f=dexp(-y**2)
  case(4)
     !cosine with a=length, b=frequency
     length=a
     frequency=b
     r=frequency*pi*x/length
     f=dcos(r)
     f2=-(frequency*pi/length)**2*dcos(r)
  case(5)
     !exp of a cosine, a=length
     nu=2.d0
     r=pi*nu/a*x
     y=dcos(r)
     yp=dsin(r)
     f=dexp(y)
     factor=(pi*nu/a)**2*(-y+yp**2)
     f2=factor*f
  case(6)
     !gaussian times "shrinked" gaussian, sigma=length/10
     length=a
     r=pi*x/length
     y=dtan(r)
     yp=pi/length*1.d0/(dcos(r))**2
     ys=2.d0*pi/length*y*yp
     factor=-2.d0*ys*y-2.d0*yp**2+4.d0*yp**2*y**2
     g=dexp(-y**2)
     g1=-2.d0*y*yp*g
     g2=factor*dexp(-y**2)
     
     sigma=length/10
     agauss=0.5d0/sigma**2
     r2=agauss*x**2
     h=dexp(-r2)
     h1=-2.d0*agauss*x*h
     h2=(-2.d0*agauss+4.d0*agauss*r2)*dexp(-r2)
     f=g*h
     f2=g2*h+g*h2+2.d0*g1*h1
  case(7)
     !sine with a=length, b=frequency
     length=a
     frequency=b
     r=frequency*pi*x/length
     f=dsin(r)
     f2=-(frequency*pi/length)**2*dsin(r)
  end select

end subroutine functions

!!***

!!$!fake ABINIT subroutines
!!$subroutine wrtout(unit,message,mode_paral)
!!$  implicit none
!!$
!!$  !Arguments ------------------------------------
!!$  integer,intent(in) :: unit
!!$  character(len=4),intent(in) :: mode_paral
!!$  character(len=500),intent(inout) :: message
!!$
!!$  print *,message
!!$end subroutine wrtout
!!$
!!$subroutine leave_new(mode_paral)
!!$  implicit none
!!$
!!$  !Arguments ------------------------------------
!!$  character(len=4),intent(in) :: mode_paral
!!$
!!$  print *,'exiting...'
!!$  stop
!!$end subroutine leave_new


