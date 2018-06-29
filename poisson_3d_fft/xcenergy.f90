
!! Copyright (C) 2002-2007 BigDFT group 
!! This file is distributed under the terms of the
!! GNU General Public License, see ~/COPYING file
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the list of contributors, see ~/AUTHORS 


!!****h* BigDFT/xcenergy
!! NAME
!!    xcenergy
!!
!! FUNCTION
!!    Calculate the XC terms from the given density in a distributed way.
!!    it assign also the proper part of the density to the zf array 
!!    which will be used for the core of the FFT procedure.
!!    Following the values of ixc and of sumpion, the array pot_ion is either summed or assigned
!!    to the XC potential, or even ignored.
!!
!! SYNOPSIS
!!    geocode  Indicates the boundary conditions (BC) of the problem:
!!            'F' free BC, isolated systems.
!!                The program calculates the solution as if the given density is
!!                "alone" in R^3 space.
!!            'S' surface BC, isolated in y direction, periodic in xz plane                
!!                The given density is supposed to be periodic in the xz plane,
!!                so the dimensions in these direction mus be compatible with the FFT
!!                Beware of the fact that the isolated direction is y!
!!            'P' periodic BC.
!!                The density is supposed to be periodic in all the three directions,
!!                then all the dimensions must be compatible with the FFT.
!!                No need for setting up the kernel.
!!    m1,m2,m3    global dimensions in the three directions.
!!    md1,md2,md3 dimensions of the arrays compatible with the FFT in the three directions.
!!    nproc       number of processors
!!    iproc       label of the process,from 0 to nproc-1
!!    ixc         eXchange-Correlation code. Indicates the XC functional to be used 
!!                for calculating XC energies and potential. 
!!                ixc=0 indicates that no XC terms are computed. The XC functional codes follow
!!                the ABINIT convention.
!!    hx,hy,hz    grid spacings. 
!!    rhopot      density in the distributed format.
!!    karray      kernel of the poisson equation. It is provided in distributed case, with
!!                dimensions that are related to the output of the PS_dim4allocation routine
!!                it MUST be created by following the same geocode as the Poisson Solver.
!!    pot_ion     additional external potential that is added to the output, 
!!                when the XC parameter ixc/=0. It is always provided in the distributed form,
!!                clearly without the overlapping terms which are needed only for the XC part
!!    exc,vxc     XC energy and integral of $\rho V_{xc}$ respectively
!!    nxc         value of the effective distributed dimension in the third direction
!!    nwb         enlarged dimension for calculating the WB correction
!!    nxt         enlarged dimension for calculating the GGA case 
!!                (further enlarged for compatibility with WB correction if it is the case)
!!    nwbl,nwbr
!!    nxcl,nxcr   shifts in the three directions to be compatible with the relation
!!                nxc+nxcl+nxcr-2=nwb, nwb+nwbl+nwbr=nxt.
!!    sumpion     logical value which states whether to sum pot_ion to the final result or not
!!                if sumpion==.true. zfionxc will be pot_ion+vxci
!!                if sumpion==.false. zfionxc will be vxci
!!                this value is ignored when ixc=0. In that case zfionxc is untouched
!!    zf          output array corresponding to the density which can be passed to FFT part
!!    zfionxc     output array which will contain pot_ion+vxci or vxci, following sumpion
!! WARNING
!!    The dimensions of pot_ion must be compatible with geocode, datacode, nproc, 
!!    ixc and iproc. Since the arguments of these routines are indicated with the *, it
!!    is IMPERATIVE to refer to PSolver routine for the correct allocation sizes.
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    February 2007
!!
!! SOURCE
!!
subroutine xc_energy(geocode,m1,m2,m3,md1,md2,md3,nxc,nwb,nxt,nwbl,nwbr,&
     nxcl,nxcr,ixc,hx,hy,hz,rhopot,pot_ion,sumpion,zf,zfionxc,exc,vxc,iproc,nproc,nspden)

  implicit none

  !Arguments----------------------
  character(len=1), intent(in) :: geocode
  logical, intent(in) :: sumpion
  integer, intent(in) :: m1,m2,m3,nxc,nwb,nxcl,nxcr,nxt,md1,md2,md3,ixc,iproc,nproc,nspden
  integer, intent(in) :: nwbl,nwbr
  real(kind=8), intent(in) :: hx,hy,hz
  real(kind=8), dimension(m1,m3,nxt,nspden), intent(inout) :: rhopot
  real(kind=8), dimension(*), intent(in) :: pot_ion
  real(kind=8), dimension(md1,md3,md2/nproc), intent(out) :: zf
  real(kind=8), dimension(md1,md3,md2/nproc,nspden), intent(out) :: zfionxc
  real(kind=8), intent(out) :: exc,vxc

  !Local variables----------------
  real(kind=8), dimension(:,:,:), allocatable :: exci,d2vxci
  real(kind=8), dimension(:,:,:,:), allocatable :: vxci,dvxci,dvxcdgr
  real(kind=8), dimension(:,:,:,:,:), allocatable :: gradient
  real(kind=8) :: elocal,vlocal,rho,potion,hgrid,sfactor
  integer :: npts,i_all,order,offset,i_stat,ispden
  integer :: i1,i2,i3,j1,j2,j3,jp2,jpp2,jppp2
  integer :: ndvxc,nvxcdgr,ngr2

  !interface with drivexc
  interface
     subroutine drivexc(exc,ixc,npts,nspden,order,rho_updn,vxc,ndvxc,ngr2,nvxcdgr,&
          dvxc,d2vxc,grho2_updn,vxcgr,exexch)    !Optional arguments 
       implicit none
       !Arguments ------------------------------------
       !scalars
       integer,intent(in) :: ixc,ndvxc,ngr2,npts,nspden,nvxcdgr,order
       integer,intent(in),optional :: exexch
       !arrays
       real(kind=8),intent(in) :: rho_updn(npts,nspden)
       real(kind=8),intent(in),optional :: grho2_updn(npts,ngr2)
       real(kind=8),intent(out) :: exc(npts),vxc(npts,nspden)
       real(kind=8),intent(out),optional :: d2vxc(npts),dvxc(npts,ndvxc)
       real(kind=8),intent(out),optional :: vxcgr(npts,nvxcdgr)
     end subroutine drivexc
  end interface

  !Body

  !check for the dimensions
  if (  nwb/=nxcl+nxc+nxcr-2 .or. nxt/=nwbr+nwb+nwbl) then
     print *,'the XC dimensions are not correct'
     print *,'nxc,nwb,nxt,nxcl,nxcr,nwbl,nwbr',nxc,nwb,nxt,nxcl,nxcr,nwbl,nwbr
     stop
  end if

  !these are always the same
!  nspden=1
  order=1
  
  !useful for the freeBC case
  hgrid=max(hx,hy,hz)

  !starting point of the density array for the GGA cases in parallel
  offset=nwbl+1
  if (ixc/=0) then
     !divide by two the density to applicate it in the ABINIT xc routines
     if(nspden==1) then
        do i3=1,nxt
           do i2=1,m3
              do i1=1,m1
                 rhopot(i1,i2,i3,nspden)=.5d0*rhopot(i1,i2,i3,nspden)
              end do
           end do
        end do
     end if
!     rewind(301)
!     do ispden=1,nspden
!        do i3=1,nxt
!           do i2=1,m3
!              do i1=1,m1
!                 write(301,'(f18.12)') rhopot(i1,i2,i3,ispden)
!              end do
!           end do
!        end do
!     end do

     !Allocations of the exchange-correlation terms, depending on the ixc value
     call size_dvxc(ixc,ndvxc,ngr2,nspden,nvxcdgr,order)

     if (ixc >= 11 .and. ixc <= 16) then
        !computation of the gradient
        allocate(gradient(m1,m3,nwb,2*nspden-1,0:3),stat=i_stat)
        call memocc(i_stat,product(shape(gradient))*kind(gradient),'gradient','xc_energy')

        !the calculation of the gradient will depend on the geometry code
        if (geocode=='F') then
           call calc_gradient(m1,m3,nxt,nwb,nwbl,nwbr,rhopot,nspden,&
                hgrid,hgrid,hgrid,gradient)
        else
        print *,'geocode=',geocode,&
             ':the calculation of the gradient is still to be performed in this case'
        stop
        end if

     end if

     !Allocations
     allocate(exci(m1,m3,nwb),stat=i_stat)
     call memocc(i_stat,product(shape(exci))*kind(exci),'exci','xc_energy')
     allocate(vxci(m1,m3,nwb,nspden),stat=i_stat)
     call memocc(i_stat,product(shape(vxci))*kind(vxci),'vxci','xc_energy')

     if (ndvxc/=0) then
        allocate(dvxci(m1,m3,nwb,ndvxc),stat=i_stat)
        call memocc(i_stat,product(shape(dvxci))*kind(dvxci),'dvxci','xc_energy')
     end if
     if (nvxcdgr/=0) then
        allocate(dvxcdgr(m1,m3,nwb,nvxcdgr),stat=i_stat)
        call memocc(i_stat,product(shape(dvxcdgr))*kind(dvxcdgr),'dvxcdgr','xc_energy')
     end if
     if ((ixc==3 .or. (ixc>=7 .and. ixc<=15)) .and. order==3) then
        allocate(d2vxci(m1,m3,nwb),stat=i_stat)
        call memocc(i_stat,product(shape(d2vxci))*kind(d2vxci),'d2vxci','xc_energy')
     end if

     if (.not.allocated(gradient) .and. nxc/=nxt ) then
        print *,'xc_energy: if nxt/=nxc the gradient must be allocated'
        stop
     end if

     !this part can be commented out if you don't want to use ABINIT modules
     !of course it must be substituted with an alternative XC calculation
     npts=m1*m3*nwb
     !let us apply ABINIT routines
     !case with gradient
     if (ixc >= 11 .and. ixc <= 16) then
      if (order**2 <= 1 .or. ixc == 16) then
         if (ixc /= 13) then             
           call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset,1),vxci,ndvxc,ngr2,nvxcdgr,&
                 &grho2_updn=gradient,vxcgr=dvxcdgr) 
         else
           call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset,1),vxci,ndvxc,ngr2,nvxcdgr,&
                 &grho2_updn=gradient) 
         end if
      else if (order /= 3) then
         if (ixc /= 13) then             
           call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset,1),vxci,ndvxc,ngr2,nvxcdgr,&
                 &dvxc=dvxci,grho2_updn=gradient,vxcgr=dvxcdgr) 
         else
           call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset,1),vxci,ndvxc,ngr2,nvxcdgr,&
                 &dvxc=dvxci,grho2_updn=gradient) 
         end if
      else if (order == 3) then
         if (ixc /= 13) then             
           call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset,1),vxci,ndvxc,ngr2,nvxcdgr,&
                 &dvxc=dvxci,d2vxc=d2vxci,grho2_updn=gradient,vxcgr=dvxcdgr) 
         else
           call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset,1),vxci,ndvxc,ngr2,nvxcdgr,&
                 &dvxc=dvxci,d2vxc=d2vxci,grho2_updn=gradient) 
         end if
      end if

     !do not calculate the White-Bird term in the Leeuwen Baerends XC case
        if (ixc/=13 .and. geocode == 'F') then
           call vxcpostprocessing(m1,m3,nwb,nxc,nxcl,nxcr,nspden,nvxcdgr,gradient,&
                hgrid,hgrid,hgrid,dvxcdgr,vxci)
        end if

        !restore the density array in the good position if it was shifted for the parallel GGA
        if (nspden==2 .and. nxt /= nwb) then
           j3=nwb+1
           do i3=nwb-nwbr,1,-1
              j3=j3-1
              do i2=1,m3
                 do i1=1,m1
                    rhopot(i1,i2,nwbl+j3,2)=rhopot(i1,i2,i3,2)
                 end do
              end do
           end do
           do i3=nxt,nwb+nwbl+1,-1 !we have nwbr points
              j3=j3-1
              do i2=1,m3
                 do i1=1,m1
                    rhopot(i1,i2,nwbl+j3,2)=rhopot(i1,i2,i3,1)
                 end do
              end do
           end do
        end if

        !cases without gradient
     else
        if (order**2 <=1 .or. ixc >= 31 .and. ixc<=34) then
           call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset,1),vxci,ndvxc,ngr2,nvxcdgr)
        else if (order==3 .and. (ixc==3 .or. ixc>=7 .and. ixc<=10)) then
           call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset,1),vxci,ndvxc,ngr2,nvxcdgr,&
                &dvxc=dvxci,d2vxc=d2vxci)
        else
           call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset,1),vxci,ndvxc,ngr2,nvxcdgr,&
                &dvxc=dvxci)
        end if
     end if
     !end of the part that can be commented out

     if (allocated(dvxci)) then
        i_all=-product(shape(dvxci))*kind(dvxci)
        deallocate(dvxci,stat=i_stat)
        call memocc(i_stat,i_all,'dvxci','xc_energy')
     end if
     if (allocated(dvxcdgr)) then
        i_all=-product(shape(dvxcdgr))*kind(dvxcdgr)
        deallocate(dvxcdgr,stat=i_stat)
        call memocc(i_stat,i_all,'dvxcdgr','xc_energy')
     end if
     if (allocated(d2vxci)) then
        i_all=-product(shape(d2vxci))*kind(d2vxci)
        deallocate(d2vxci,stat=i_stat)
        call memocc(i_stat,i_all,'d2vxci','xc_energy')
     end if
     if (allocated(gradient)) then
        i_all=-product(shape(gradient))*kind(gradient)
        deallocate(gradient,stat=i_stat)
        call memocc(i_stat,i_all,'gradient','xc_energy')
     end if

!     rewind(300)
!     do ispden=1,nspden
!        do i3=1,nxt
!           do i2=1,m3
!              do i1=1,m1
!                 write(300,'(f18.12)') rhopot(i1,i2,i3,ispden)
!              end do
!           end do
!        end do
!     end do


     exc=0.d0
     vxc=0.d0
     sfactor=1.0d0
     if(nspden==1) sfactor=2.0d0
     if (sumpion) then
        !summing the xc potential into the zfionxc array with pot_ion
        ispden=1
        do jp2=1,nxc
           j2=offset+jp2+nxcl-2
           jppp2=jp2+nxcl-1
           do j3=1,m3
              do j1=1,m1
                 jpp2=j1+(j3-1)*m1+(jp2-1)*m1*m3
                 rho=rhopot(j1,j3,j2,ispden)
                 potion=pot_ion(jpp2)
                 elocal=exci(j1,j3,jppp2)
                 vlocal=vxci(j1,j3,jppp2,ispden)
                 exc=exc+elocal*rho
                 vxc=vxc+vlocal*rho
                 zf(j1,j3,jp2)=sfactor*rho !restore the original normalization
                 zfionxc(j1,j3,jp2,ispden)=potion+vlocal
              end do
              do j1=m1+1,md1
                 zf(j1,j3,jp2)=0.d0
                 zfionxc(j1,j3,jp2,ispden)=0.d0
              end do
           end do
           do j3=m3+1,md3
              do j1=1,md1
                 zf(j1,j3,jp2)=0.d0
                 zfionxc(j1,j3,jp2,ispden)=0.d0
              end do
           end do
        end do
        do jp2=nxc+1,md2/nproc
           do j3=1,md3
              do j1=1,md1
                 zf(j1,j3,jp2)=0.d0
                 zfionxc(j1,j3,jp2,ispden)=0.d0
              end do
           end do
        end do
        !spin-polarised case
        if (nspden==2) then
           ispden=2
           do jp2=1,nxc
              j2=offset+jp2+nxcl-2
              jppp2=jp2+nxcl-1
              do j3=1,m3
                 do j1=1,m1
                    jpp2=j1+(j3-1)*m1+(jp2-1)*m1*m3
                    rho=rhopot(j1,j3,j2,ispden)
                    potion=pot_ion(jpp2)
                    elocal=exci(j1,j3,jppp2)
                    vlocal=vxci(j1,j3,jppp2,ispden)
                    exc=exc+elocal*rho
                    vxc=vxc+vlocal*rho
                    zf(j1,j3,jp2)=zf(j1,j3,jp2)+sfactor*rho !restore original normalization
                    zfionxc(j1,j3,jp2,ispden)=potion+vlocal
                 end do
                 do j1=m1+1,md1
                    zfionxc(j1,j3,jp2,ispden)=0.d0
                 end do
              end do
              do j3=m3+1,md3
                 do j1=1,md1
                    zfionxc(j1,j3,jp2,ispden)=0.d0
                 end do
              end do
           end do
           do jp2=nxc+1,md2/nproc
              do j3=1,md3
                 do j1=1,md1
                    zfionxc(j1,j3,jp2,ispden)=0.d0
                 end do
              end do
           end do
        end if
     else
        !the zfionxc aray contain only the XC potential. pot_ion is ignored
        ispden=1
        do jp2=1,nxc
           j2=offset+jp2+nxcl-2
           jppp2=jp2+nxcl-1
           do j3=1,m3
              do j1=1,m1
                 rho=rhopot(j1,j3,j2,ispden)
                 elocal=exci(j1,j3,jppp2)
                 vlocal=vxci(j1,j3,jppp2,ispden)
                 exc=exc+elocal*rho
                 vxc=vxc+vlocal*rho
                 zf(j1,j3,jp2)=sfactor*rho!restore the original normalization
                 zfionxc(j1,j3,jp2,ispden)=vlocal
              end do
              do j1=m1+1,md1
                 zf(j1,j3,jp2)=0.d0
                 zfionxc(j1,j3,jp2,ispden)=0.d0
              end do
           end do
           do j3=m3+1,md3
              do j1=1,md1
                 zf(j1,j3,jp2)=0.d0
                 zfionxc(j1,j3,jp2,ispden)=0.d0
              end do
           end do
        end do
        do jp2=nxc+1,md2/nproc
           do j3=1,md3
              do j1=1,md1
                 zf(j1,j3,jp2)=0.d0
                 zfionxc(j1,j3,jp2,ispden)=0.d0
              end do
           end do
        end do
        !spin-polarised case
        if (nspden==2) then
           ispden=2
           do jp2=1,nxc
              j2=offset+jp2+nxcl-2
              jppp2=jp2+nxcl-1
              do j3=1,m3
                 do j1=1,m1
                    rho=rhopot(j1,j3,j2,ispden)
                    elocal=exci(j1,j3,jppp2)
                    vlocal=vxci(j1,j3,jppp2,ispden)
                    exc=exc+elocal*rho
                    vxc=vxc+vlocal*rho
                    zf(j1,j3,jp2)=zf(j1,j3,jp2)+sfactor*rho!restore the original normalization
                    zfionxc(j1,j3,jp2,ispden)=vlocal
                 end do
                 do j1=m1+1,md1
                    zfionxc(j1,j3,jp2,ispden)=0.d0
                 end do
              end do
              do j3=m3+1,md3
                 do j1=1,md1
                    zfionxc(j1,j3,jp2,ispden)=0.d0
                 end do
              end do
           end do
           do jp2=nxc+1,md2/nproc
              do j3=1,md3
                 do j1=1,md1
                    zfionxc(j1,j3,jp2,ispden)=0.d0
                 end do
              end do
           end do
        end if
     end if
     !the two factor is due to the 
     !need of using the density of states in abinit routines
     exc=sfactor*hx*hy*hz*exc
     vxc=sfactor*hx*hy*hz*vxc


     !De-allocations
     i_all=-product(shape(exci))*kind(exci)
     deallocate(exci,stat=i_stat)
     call memocc(i_stat,i_all,'exci','xc_energy')
     i_all=-product(shape(vxci))*kind(vxci)
     deallocate(vxci,stat=i_stat)
     call memocc(i_stat,i_all,'vxci','xc_energy')


  else

     !case without XC terms
     !distributing the density in the zf array
     exc=0.d0
     vxc=0.d0
     do jp2=1,nxc
        j2=offset+jp2+nxcl-2
        jpp2=jp2     
        do j3=1,m3
           do j1=1,m1
              zf(j1,j3,jp2)=rhopot(j1,j3,j2,1)
           end do
           do j1=m1+1,md1
              zf(j1,j3,jp2)=0.d0
           end do
        end do
        do j3=m3+1,md3
           do j1=1,md1
              zf(j1,j3,jp2)=0.d0
           end do
        end do
     end do
     do jp2=nxc+1,md2/nproc
        do j3=1,md3
           do j1=1,md1
              zf(j1,j3,jp2)=0.d0
           end do
        end do
     end do

  end if

end subroutine xc_energy

!!****f* BigDFT/vxcpostprocessing
!! NAME
!! vxcpostprocessing
!!
!! FUNCTION
!! Correct the XC potential with the White-Bird formula, to be used for the 
!! GGA case. Works either in parallel of in serial, by proper change of the 
!! arguments.
!!
!! SOURCE
subroutine vxcpostprocessing(n01,n02,n03,n3eff,wbl,wbr,nspden,nvxcdgr,gradient,hx,hy,hz,dvxcdgr,wb_vxc)
  implicit none
  integer, intent(in) :: n01,n02,n03,n3eff,wbl,wbr,nspden,nvxcdgr
  real(kind=8), intent(in) :: hx,hy,hz
  real(kind=8), dimension(n01,n02,n03,2*nspden-1,0:3), intent(in) :: gradient
  real(kind=8), dimension(n01,n02,n03,nvxcdgr), intent(in) :: dvxcdgr
  real(kind=8), dimension(n01,n02,n03,nspden), intent(inout) :: wb_vxc
  !Local variables
  integer :: i1,i2,i3,dir_i,i_all,i_stat
  real(kind=8) :: dnexcdgog,grad_i,rho_up,rho_down,rho_tot
  real(kind=8), dimension(:,:,:,:,:), allocatable :: f_i

  !Body

  allocate(f_i(n01,n02,n03,3,nspden),stat=i_stat)
  call memocc(i_stat,product(shape(f_i))*kind(f_i),'f_i','vxcpostprocessing')
  
  !let us first treat the case nspden=1
  if (nspden == 1) then
     !Let us construct the object we have to manipulate with another gradient
     if (nvxcdgr == 3) then
        do dir_i=1,3
           !Let us construct the object we have to manipulate with another gradient
           do i3=1,n03
              do i2=1,n02
                 do i1=1,n01
                    dnexcdgog=0.5d0*dvxcdgr(i1,i2,i3,1) + dvxcdgr(i1,i2,i3,3)
                    grad_i=2.d0*gradient(i1,i2,i3,1,dir_i)
                    f_i(i1,i2,i3,dir_i,1)=dnexcdgog*grad_i
                 end do
              end do
           end do
        end do
     else
        do dir_i=1,3
           !Let us construct the object we have to manipulate with another gradient
           do i3=1,n03
              do i2=1,n02
                 do i1=1,n01
                    dnexcdgog=0.5d0*dvxcdgr(i1,i2,i3,1)
                    grad_i=2.d0*gradient(i1,i2,i3,1,dir_i)
                    f_i(i1,i2,i3,dir_i,1)=dnexcdgog*grad_i
                 end do
              end do
           end do
        end do
     end if

  !then the spin-polarized case
  else

     if (nvxcdgr == 3) then
        do dir_i=1,3
           do i3=1,n03
              do i2=1,n02
                 do i1=1,n01
                    rho_up=gradient(i1,i2,i3,1,dir_i)  !rho_ instead of grad_ for ABINIT comp.
                    rho_down=gradient(i1,i2,i3,2,dir_i)
                    rho_tot=gradient(i1,i2,i3,3,dir_i)
                    f_i(i1,i2,i3,dir_i,1)=rho_up*dvxcdgr(i1,i2,i3,1)+&
                         rho_tot*dvxcdgr(i1,i2,i3,3)
                    f_i(i1,i2,i3,dir_i,2)=rho_down*dvxcdgr(i1,i2,i3,2)+&
                         rho_tot*dvxcdgr(i1,i2,i3,3)
                 end do
              end do
           end do
        end do
     else
        do dir_i=1,3
           do i3=1,n03
              do i2=1,n02
                 do i1=1,n01
                    rho_up=gradient(i1,i2,i3,1,dir_i)
                    rho_down=gradient(i1,i2,i3,2,dir_i)
                    f_i(i1,i2,i3,dir_i,1)=rho_up*dvxcdgr(i1,i2,i3,1)
                    f_i(i1,i2,i3,dir_i,2)=rho_down*dvxcdgr(i1,i2,i3,2)
                 end do
              end do
           end do
        end do
     end if

  !end of spin-polarized if statement
  end if

  !let us now calculate the gradient and correct the result
  call wb_correction(n01,n02,n03,n3eff,wbl,wbr,f_i,&
       hx,hy,hz,nspden,wb_vxc)


  i_all=-product(shape(f_i))*kind(f_i)
  deallocate(f_i,stat=i_stat)
  call memocc(i_stat,i_all,'f_i','vxcpostprocessing')

end subroutine vxcpostprocessing
