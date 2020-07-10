! Raytracing module for MPI-AMRVAC
! Has been tested in 1D and 2D.
!
! Allard Jan van Marle
!
!
! Last update 07/01/2013
!


!=============================================================================
module mod_raytracing
!
! Contains global variables for raytracing
!
  implicit none


 integer             :: nrays         ! nr of rays
 integer             :: npoints       ! nr of points per ray
 integer,allocatable :: npointsmax(:) ! highest point in use along each ray

 integer,allocatable :: nrayp(:)              ! number of rays per processor
 integer,allocatable :: minrayp(:),maxrayp(:) ! min and max number of rays per processor 


 integer          :: nu        ! nr of variables that are unique to the rays
 integer          :: cdens_    ! index oc column density
 integer          :: ion_      ! index of ionization state

 character*79     :: raysource,raytype,gridtoray,raytogrid &
                     ,raycomm,quicksearch,rayphysics
		     
 logical          :: writerays,checkintersect,regular_rays

 double precision,allocatable :: time_ray_interpol(:)  ! time spent interpolating between grid and rays
 double precision,allocatable :: time_ray_comm(:)      ! time spent communicating rays between processors
 double precision,allocatable :: time_ray_phys(:)      ! time spent doing physics on rays
 double precision,allocatable :: time_ray_init(:)      ! time spent initializing rays
 

 double precision :: xsource^D  ! location of source for pointsources

 double precision,allocatable :: xraystart(:,:),xrayend(:,:) ! first and last point of each ray
 double precision,allocatable :: xraypoint(:,:,:)            ! position of each point
 double precision,allocatable :: xraymin(:,:),xraymax(:,:)   ! min and max position of each ray per dimension
 double precision,allocatable :: wraypoint(:,:,:)            ! grid data on point

 double precision,allocatable :: uraypoint(:,:,:)            ! variables that exist only on the ray
 double precision,allocatable :: rayvar(:,:)                 ! variable that is defined per ray

 logical, allocatable :: ptfound(:,:)                        ! flag for having found a point in a grid
 logical, allocatable :: ptgathered(:,:)                     ! flag for having gathered a point onto the central processor


 character*79 :: rayfilename


end module mod_raytracing

!=============================================================================

subroutine init_rays
!
! Initializes rays. Has to be called from 'init_global data'
!
use mod_raytracing

use mod_global_parameters


 double precision :: drayx^D
 double precision, allocatable :: raylength(:)

 integer :: iray, ipoint, ip

 integer :: fileid
 character(79) :: filename
 
 double precision :: time1,time2

!----------------------------------------------------------------------------
time1 = MPI_WTIME()

if (.not. (fixprocess) ) call MPISTOP( 'getting the grid data to the rays would be a good idea!' )


raysource = TRIM(raysource)
raytype   = TRIM(raytype)
raytogrid = TRIM(raytogrid)
gridtoray = TRIM(gridtoray)
raycomm   = TRIM(raycomm)

allocate( xraystart(1:nrays,1:ndim) )
allocate( xrayend(1:nrays,1:ndim) )
allocate( xraypoint(1:nrays,1:npoints,1:ndim) )
allocate( xraymin(1:nrays,1:ndim) )
allocate( xraymax(1:nrays,1:ndim) )
allocate( wraypoint(1:nrays,1:npoints,1:nw) )
allocate( uraypoint(1:nrays,1:npoints,1:nu) )

allocate( rayvar(1:nrays,1:nu) )

allocate( ptfound(1:nrays,1:npoints) )
allocate( ptgathered(1:nrays,1:npoints) )
allocate( npointsmax(1:nrays) )
allocate( raylength(1:nrays ) )

ptfound(1:nrays,1:npoints)    = .false.
ptgathered(1:nrays,1:npoints) = .false.

wraypoint(1:nrays,1:npoints,1:nw) = zero
uraypoint(1:nrays,1:npoints,1:nu) = zero

npointsmax(1:nrays) = npoints


allocate( nrayp(0:npe-1) )
allocate( minrayp(0:npe-1) )
allocate( maxrayp(0:npe-1) )

allocate( time_ray_interpol(0:npe-1) )
allocate( time_ray_comm(0:npe-1) )
allocate( time_ray_phys(0:npe-1) )
allocate( time_ray_init(0:npe-1) )

time_ray_interpol(0:npe-1) = zero
time_ray_comm(0:npe-1)     = zero
time_ray_phys(0:npe-1)     = zero

!
! Distribute rays over processors
!
if( nrays>1) then
  nrayp(0:npe-2) = (nrays/(npe))
  nrayp(npe-1) = nrays-(npe-1)*nrayp(0)

  minrayp(0)=1
  maxrayp(0)=nrayp(0)
  do ip=1,npe-1
      minrayp(ip) = maxrayp(ip-1)+1
      maxrayp(ip) = minrayp(ip)+nrayp(ip)-1
  enddo
else
  nrayp(0)   = 1
  minrayp(0) = 1
  maxrayp(0) = 1
  if( npe>1 ) then
      nrayp(1:npe-1)   = 0
      minrayp(1:npe-1) = 0
      maxrayp(1:npe-1) = 0
  endif
endif

 call MPI_BARRIER(icomm,ierrmpi)
!
!write(*,*) 'I am proc. ', mype,' and I deal with rays ', minrayp(mype), ' through ', maxrayp(mype)
!

!
! Define the end poins for each ray
!

select case( raysource )
 case( 'point' )
 {^IFTWOD

      {^D&xraystart(1:nrays,^D)=xsource^D;} 
      
      xraystart(nrays/2+1,1) = xsource1    
      xraystart(nrays/2+1,2) = xsource2 
      xrayend(nrays/2+1,1)   = xprobmax1 - smalldouble 
      xrayend(nrays/2+1,2)   = xprobmax2 - smalldouble 
      raylength(nrays/2+1)   = dsqrt((xrayend(nrays/2+1,1)-xraystart(nrays/2+1,1))**2 &
                            + (xrayend(nrays/2+1,2)-xraystart(nrays/2+1,2))**2 )
               
      npointsmax(nrays/2+1)  = npoints

      drayx2 = (xprobmax2-xprobmin2)/(nrays/2)
      drayx1 = (xprobmax1-xprobmin1)/(nrays/2)
            
      
      do iray=nrays/2,1,-1
         xraystart(iray,1) = xsource1 
	 xraystart(iray,2) = xsource2       
         xrayend(iray,1)   = max(xrayend(iray+1,1) - drayx1,xprobmin1+smalldouble)
	 xrayend(iray,2)   = xprobmax2 - smalldouble      
         raylength(iray)      = max( smalldouble,dsqrt((xrayend(iray,1)-xraystart(iray,1))**2 &
                           + (xrayend(iray,2)-xraystart(iray,2))**2 ) )
      
         npointsmax(iray) = max( 2,INT(npoints *raylength(iray)/raylength(nrays/2+1)) )	     
      enddo

      do iray=nrays/2+2,nrays
         xraystart(iray,1) = xsource1
	 xraystart(iray,2) = xsource2 
         xrayend(iray,1)   = xprobmax1 - smalldouble	 
	 xrayend(iray,2)   = max(xrayend(iray-1,2) - drayx2,xprobmin2+smalldouble)
         raylength(iray)   = max( smalldouble,dsqrt((xrayend(iray,1)-xraystart(iray,1))**2 &
                           + (xrayend(iray,2)-xraystart(iray,2))**2 ) )
      
         npointsmax(iray) = max( 2,INT(npoints *raylength(iray)/raylength(nrays/2+1)) )	     
      enddo      
}     
           
 case( 'x1dir' )
{^IFONED
      xraystart(1,1) = xprobmin1 + smalldouble
      xrayend(1,1)   = xprobmax1 - smalldouble   
      raylength(1)   = dabs( xrayend(1,1) - xraystart(1,1) )

      npointsmax(1)  = npoints
}
{^IFTWOD
      xraystart(1,1)     = xprobmin1 + smalldouble   
      xraystart(1,2)     = xprobmin2 + smalldouble   
      xraystart(nrays,1) = xprobmin1 + smalldouble
      xraystart(nrays,2) = xprobmax2 - smalldouble
      xrayend(1,1)       = xprobmax1 - smalldouble   
      xrayend(1,2)       = xprobmin2 + smalldouble   
      xrayend(nrays,1)   = xprobmax1 - smalldouble
      xrayend(nrays,2)   = xprobmax2 - smalldouble

      drayx2 = (xraystart(nrays,2) - xraystart(1,2))/(nrays-1)
      
      do iray = 2,nrays-1
      	 xraystart(iray,1) = xprobmin1 + smalldouble
         xraystart(iray,2) = xraystart(iray-1,2) + drayx2
         xrayend(iray,1)   = xprobmax1 - smalldouble	 
	 xrayend(iray,2)   = xrayend(iray-1,2) + drayx2	
      enddo
}
{^IFTHREED
   call mpistop( 'Under construction' )
}
      npointsmax(1:nrays) = npoints

 case( 'x2dir' )
{^IFTWOD
      xraystart(1,1)     = xprobmin1 + smalldouble  
      xraystart(1,2)     = xprobmin2 + smalldouble 
      xraystart(nrays,1) = xprobmax1 - smalldouble
      xraystart(nrays,2) = xprobmin2 + smalldouble
      xrayend(1,1)       = xprobmin1 + smalldouble  
      xrayend(1,2)       = xprobmax2 - smalldouble   
      xrayend(nrays,1)   = xprobmax1 - smalldouble
      xrayend(nrays,2)   = xprobmax2 - smalldouble     
     
      drayx1 = (xraystart(nrays,1) - xraystart(1,1))/(nrays-1)

      do iray = 2,nrays-1
         xraystart(iray,1) = xraystart(iray-1,1) + drayx1
	 xraystart(iray,2) = xprobmin2 + smalldouble
         xrayend(iray,1)   = xrayend(iray-1,1) + drayx1	 
	 xrayend(iray,2)   = xprobmax2 - smalldouble
      enddo     
}
{^IFTHREED
   call mpistop( 'Under construction' )
}

      npointsmax(1:nrays) = npoints 
 case( 'angle' )
{^IFTWOD
      if( mod( nrays,2 )==0 )then 
         call MPISTOP( 'Number of rays must be uneven for this distribution' )
      endif
       
      xraystart(nrays/2+1,1) = xprobmin1 + smalldouble   
      xraystart(nrays/2+1,2) = xprobmin2 + smalldouble   
      xrayend(nrays/2+1,1)   = xprobmax1 - smalldouble 
      xrayend(nrays/2+1,2)   = xprobmax2 - smalldouble 
      raylength(nrays/2+1)   = dsqrt((xrayend(nrays/2+1,1)-xraystart(nrays/2+1,1))**2 &
                            + (xrayend(nrays/2+1,2)-xraystart(nrays/2+1,2))**2 )
               
      npointsmax(nrays/2+1)  = npoints

      drayx2 = (xprobmax2-xprobmin2)/(nrays/2)
      drayx1 = (xprobmax1-xprobmin1)/(nrays/2)            
      
      do iray=nrays/2,1,-1
         xraystart(iray,1) = xprobmin1 + smalldouble
	 xraystart(iray,2) = min( xraystart(iray+1,2) + drayx2,xprobmax2-smalldouble )
         xrayend(iray,1)   = max( xrayend(iray+1,1) - drayx1,xprobmin1+smalldouble )
	 xrayend(iray,2)   = xprobmax2 - smalldouble      
         raylength(iray)      = max( smalldouble,dsqrt((xrayend(iray,1)-xraystart(iray,1))**2 &
                           + (xrayend(iray,2)-xraystart(iray,2))**2 ) )
      
         npointsmax(iray) = max( 1,INT(npoints *raylength(iray)/raylength(nrays/2+1)) )	     
      enddo

      do iray=nrays/2+2,nrays
         xraystart(iray,1) = min( xraystart(iray-1,1) + drayx1,xprobmax1-smalldouble )
	 xraystart(iray,2) = xprobmin2 + smalldouble
         xrayend(iray,1)   = xprobmax1 - smalldouble	 
	 xrayend(iray,2)   = max( xrayend(iray-1,2) - drayx2,xprobmin2+smalldouble )      
         raylength(iray)   = max( smalldouble,dsqrt((xrayend(iray,1)-xraystart(iray,1))**2 &
                           + (xrayend(iray,2)-xraystart(iray,2))**2 ) )
      
         npointsmax(iray) = max( 1,INT(npoints *raylength(iray)/raylength(nrays/2+1)) )	     
      enddo
}       
 case( 'usrraysource' )
    call mpistop( 'This is where you have to specify your own source' )
 case default
    call mpistop( 'This type of ray has not been defined' )
end select    

!
! Determine the points along each ray
!
select case( raytype )
 case( 'fixed' )
 
 do iray = 1,nrays
   {^D&xraypoint(iray,1,^D) = xraystart(iray,^D);}
   {^D&xraypoint(iray,npoints,^D) = xrayend(iray,^D);}
   {^D&drayx^D = (xrayend(iray,^D)-xraystart(iray,^D))/(npointsmax(iray)-1);}
    do ipoint = 2,npoints-1
   {^D&xraypoint(iray,ipoint,^D) =  xraypoint(iray,ipoint-1,^D) + drayx^D;}
    enddo
   {^D&xraymin(iray,^D) = minval(xraypoint(iray,1:npointsmax(iray),^D));}
   {^D&xraymax(iray,^D) = maxval(xraypoint(iray,1:npointsmax(iray),^D));}

enddo
  
 case default
    call mpistop( 'This type of ray has not been defined' )
end select


if(mype == 0) then
   fileid = 10
   write(filename,"(a,a)") TRIM(base_filename), ".rays"
   open(fileid,file=filename, status='unknown')
   write(fileid,*) 'TITLE="RAYS"'
{^IFONED      write(10,*) 'VARIABLES = X'}
{^IFTWOD      write(10,*) 'VARIABLES = X,Y'}
{^IFTHREED    write(10,*) 'VARIABLES = X,Y,Z'}
   do iray = 1,nrays
      write(fileid,*) 'ZONE T= "',iray,'",F=POINT' 
      do ipoint = 1,npointsmax(iray)
{^IFONED    write(fileid,1001) xraypoint(iray,ipoint,1)}
{^IFTWOD    write(fileid,1001) xraypoint(iray,ipoint,1), xraypoint(iray,ipoint,2)}
{^IFTHREED  write(fileid,1001) xraypoint(iray,ipoint,1), xraypoint(iray,ipoint,2), xraypoint(iray,ipoint,3)}
      enddo
      write(fileid,*)
   enddo
   close(fileid)
endif

deallocate( raylength )


time2 = MPI_WTIME()
time_ray_init(mype) = time2-time1

return
1001 format(4(1x,1pe12.5))
end subroutine init_rays

!==============================================================================

subroutine getgriddata(ixI^L,ixO^L,x,w)
!
! Interpolates grid data onto the ray points
!
!
! If we call this at the start of the new timestep, 
! we could take one gridcell extra on each side.
!
use mod_raytracing

use mod_global_parameters

integer, intent(in)          :: ixI^L, ixO^L
double precision, intent(in) :: x(ixI^S,1:ndim),w(ixI^S,1:nw)

double precision :: xpoint^D, dx^D,xmin^D,xmax^D

logical          :: ingrid,ingridlast
integer          :: iray,ipoint,igridpoint^D,istartlastgrid,ipointstart

double precision :: time1,time2


!----------------------------------------------------------------------------

 time1 = MPI_WTIME()

!
! Determine if the point is inside the grid
!

{^D&xmin^D=minval(x(ixO^S,^D));}
{^D&xmax^D=maxval(x(ixO^S,^D));}
{^D&dx^D=(xmax^D-xmin^D)/(ixOmax^D-ixOmin^D);}
{^D&xmin^D=xmin^D-dx^D;}
{^D&xmax^D=xmax^D+dx^D;}

 do iray = 1,nrays
   select case( TRIM(quicksearch) )
   case( 'basic' )
{^IFTWOD
       if( xraymin(iray,1) > xmax1 .or. xraymax(iray,1)<xmin1 ) cycle 
       if( xraymin(iray,2) > xmax2 .or. xraymax(iray,2)<xmin2 ) cycle
}
{^IFTHREED
           xraymin(iray,3) > xmax3 .or. xraymax(iray,3)<xmin3 ) cycle
}
   case( 'none' )     
   case default
      call MPISTOP( 'This search method has not been implemented' )
   end select

    ingridlast = .false.

    do ipoint = 1,npointsmax(iray)
       ingrid = .false.      
       {^D&xpoint^D= xraypoint(iray,ipoint,^D);}       

        ingrid = ( xpoint1>=xmin1 .and. xpoint1<=xmax1 )
        if( .not.(ingrid) ) cycle
{^IFTWOD
        ingrid = ( ingrid &	 
	        .and. xpoint2>=xmin2 .and. xpoint2<=xmax2 )
}
        if( .not.(ingrid )) cycle
{^IFTHREED
        ingrid = ( ingrid  .and. xpoint3>=xmin3 .and. xpoint3<=xmax3 )
}		
        if( .not.(ingrid ) )then 
            cycle
        else
	select case( gridtoray )
	case( 'nearest' )
	   call findnearestgridpoint( ixI^L,ixI^L,xpoint^D,x,igridpoint^D) ! call with ixI to include ghostcells
	   wraypoint(iray,ipoint,1:nw) = w(igridpoint^D,1:nw)
	case('linear' )
	   call interpolategrid( ixI^L,ixO^L,iray,ipoint,x,w )
        case default
	   call mpistop( 'This type of interpolation has not been implemented' )
	end select
	ptfound(iray,ipoint) = ingrid 
	endif    
      ! 
      ! Check if the ray has passed through the grid
      !
        if( (.not.(ingrid)) .and. ingridlast ) exit
        ingridlast=ingrid
    enddo
 enddo

time2 = MPI_WTIME()
time_ray_interpol(mype) = time_ray_interpol(mype) + (time2-time1)


return

end subroutine getgriddata

!==============================================================================

subroutine interpolategrid(ixI^L,ixO^L,iray,ipoint,x,w)
!
! Interpolation of grid values on a raypoint
!
use mod_raytracing

use mod_global_parameters

 integer, intent(in)           :: ixI^L,ixO^L, iray,ipoint
 double precision, intent(in)  :: x(ixI^S,1:ndim),w(ixI^S,1:nw)
 
 double precision :: xpoint^D
 integer :: ix^D, il^D,ih^D

 double precision :: dxgrid^D
{^IFONED double precision :: dwgrid(1:nw)}
{^IFTWOD double precision :: wtemp1(1:nw), wtemp2(1:nw), rtogrid(1:4)}
!---------------------------------------------------------------------------------

 {^D&xpoint^D = xraypoint(iray,ipoint,^D);}
!
! Find surrounding gridpoints
!
{^IFONED
do ix1=ixOmin1,ixOmax1
   if( x(ix1,1) >= xpoint1 ) then 
      il1 =ix1-1;ih1=ix1
      exit
   endif
enddo
}

{^IFTWOD
do ix1=ixOmin1,ixOmax1+2
   if( x(ix1,3,1) >= xpoint1 ) then 
      il1 =ix1-1;ih1=ix1
      exit
   endif
enddo
do ix2=ixOmin2,ixOmax2+2
   if( x(3,ix2,2) >= xpoint2 ) then 
      il2 =ix2-1;ih2=ix2
      exit
   endif
enddo
}

{^IFTHREED
do ix1=ixOmin1,ixOmax1+2
   if( x(ix1,3,3,1) >= xpoint1 ) then 
      il1 =ix1-1;ih1=ix1
      exit
   endif
enddo
do ix2=ixOmin2,ixOmax2+2
   if( x(3,ix2,3,2) >= xpoint2 ) then 
      il2 =ix2-1;ih2=ix2
      exit
   endif
enddo
do ix3=ixOmin3,ixOmax3+2
   if( x(3,3,ix3,3) >= xpoint3 ) then 
      il3 =ix3-1;ih3=ix3
      exit
   endif
enddo
}

select case( TRIM(gridtoray) )
 case( 'linear' )
{^IFONED
  dxgrid1      =  x(ih1,1) -x(il1,1)
  dwgrid(1:nw) =  w(ih1,1:nw)-w(il1,1:nw)
  wraypoint(iray,ipoint,1:nw) =  w(il1,1:nw) + (xpoint1-x(il1,1)) * dwgrid(1:nw)/dxgrid1
}
{^IFTWOD
  dxgrid1      =  x(ih1,il2,1) - x(il1,il2,1)
  dxgrid2      =  x(3,ih2,2) - x(3,il2,2)
  wtemp1(1:nw) =  (x(ih1,il2,1)-xpoint1)*w(il1,il2,1:nw) /dxgrid1  & 
               +  (xpoint1-x(il1,il2,1))*w(ih1,il2,1:nw) /dxgrid1
  wtemp2(1:nw) =  (x(ih1,il2,1)-xpoint1)*w(il1,ih2,1:nw) /dxgrid1  & 
               +  (xpoint1-x(il1,il2,1))*w(ih1,ih2,1:nw) /dxgrid1
  wraypoint(iray,ipoint,1:nw) =  (x(il1,ih2,2)-xpoint2)*wtemp1 /dxgrid2 &
               + (xpoint2-x(il1,il2,2))*wtemp1 /dxgrid2
}
{^IFTHREED
    call MPISTOP( ' under construction')
}
 case default 
    call MPISTOP( 'This interpolation method has not been implemented' )
end select

return
end subroutine interpolategrid

!==============================================================================

subroutine findnearestgridpoint(ixI^L,ixO^L,xpoint^D,x,igridpoint^D)
!
! Finds which gridpoint lies nearest to the ray point
!
use mod_raytracing

use mod_global_parameters

integer, intent(in)          :: ixI^L,ixO^L
double precision, intent(in) :: xpoint^D,x(ixI^S,1:ndim)

integer, intent(out) :: igridpoint^D

!double precision :: rtopoint(ixG^T)
{^IFONED   double precision             :: rtopoint(1:2)}
{^IFTWOD   double precision,allocatable :: rtopoint(:,:)}
{^IFTHREED double precision,allocatable :: rtopoint(:,:,:)}


integer :: ix^D, i ,ilo^D,ihi^D     !,i^T    !, il^D,ih^D
double precision :: rmin
!---------------------------------------------------------------------------------
 

{^IFONED
do ix1=ixOmin1,ixOmax1
   if( x(ix1,1) >= xpoint1 ) then 
      ilo1 =ix1-1;ihi1=ix1
      exit
   endif
enddo

rtopoint(1) = (x(ilo1,^D)-xpoint^D)**2
rtopoint(2) = (x(ihi1,^D)-xpoint^D)**2
if( rtopoint(1) <= rtopoint(2) ) then
   igridpoint1 = ilo1
else
   igridpoint1 = ihi1
endif
return
}

{^IFTWOD
do ix1=ixOmin1,ixOmax1+2
   if( x(ix1,3,1) >= xpoint1 ) then 
      ilo1 =ix1-1;ihi1=ix1
      exit
   endif
enddo
do ix2=ixOmin2,ixOmax2+2
   if( x(3,ix2,2) >= xpoint2 ) then 
      ilo2 =ix2-1;ihi2=ix2
      exit
   endif
enddo
}

{^IFTHREED
do ix1=ixOmin1,ixOmax1+2
   if( x(ix1,3,3,1) >= xpoint1 ) then 
      ilo1 =ix1-1;ihi1=ix1
      exit
   endif
enddo
do ix2=ixOmin2,ixOmax2+2
   if( x(3,ix2,3,2) >= xpoint2 ) then 
      ilo2 =ix2-1;ihi2=ix2
      exit
   endif
enddo
do ix3=ixOmin3,ixOmax3+2
   if( x(3,3,ix3,3) >= xpoint3 ) then 
      ilo3 =ix3-1;ihi3=ix3
      exit
   endif
enddo
}

{^NOONED
allocate( rtopoint({^D&ilo^D:ihi^D}) )

rtopoint(i^T) =({^D&((x(i^T,^D)-xpoint^D)**2)+})
rmin = 1.0d32
{do ix^D = ilo^D,ihi^D\}
   if( rtopoint(ix^D) <=smalldouble ) then
      {^D&igridpoint^D = ix^D;}
      exit
   else if( rtopoint(ix^D) < rmin ) then
   rmin = rtopoint(ix^D)
   {^D&igridpoint^D = ix^D;}
   endif
{enddo^D&\}


deallocate( rtopoint)
}


return
end subroutine findnearestgridpoint
!==============================================================================
!
subroutine ray_columndens
!
! Calculate columndensity along rays, using 2 point trapezoid intgeral
!
use mod_raytracing

use mod_global_parameters

double precision :: dray(1:npoints)  
integer          :: iray,ipoint
!
!---------------------------------------------------------------------------
! 
 uraypoint(minrayp(mype):maxrayp(mype),1:npoints,cdens_) = zero
 
 do iray=minrayp(mype),maxrayp(mype)
    dray(1) = 0.5d0*dsqrt({^D&((xraypoint(iray,2,^D)-xraypoint(iray,1,^D))**2)+})
    dray(npointsmax(iray)) = 0.5d0*dsqrt({^D&((xraypoint(iray,npointsmax(iray),^D) &
                           - xraypoint(iray,npointsmax(iray)-1,^D))**2)+})

    uraypoint(iray,1,cdens_) = zero
    do ipoint=2,npointsmax(iray)-1
      dray(ipoint) = dsqrt({^D&((xraypoint(iray,ipoint,^D) &
                   - xraypoint(iray,ipoint-1,^D))**2)+})
      uraypoint(iray,ipoint,cdens_ ) = uraypoint(iray,ipoint-1,cdens_ )                            &
                + 0.5d0*(wraypoint(iray,ipoint-1,rho_)+wraypoint(iray,ipoint,rho_))*dray(ipoint) 
    enddo			  
			  
    uraypoint(iray,npointsmax(iray),cdens_ ) = uraypoint(iray,npointsmax(iray)-1,cdens_ )          &
            + 0.5d0*(wraypoint(iray,npointsmax(iray),rho_)*dray(npointsmax(iray)))
 enddo
 
 
 return 
 end subroutine ray_columndens
!==============================================================================
!
subroutine ray_stromgren
!
! Calculate Stromgren radius along ray, using 2 point trapezoid intgeral
! From Garcia-Segura & Franco 1996
!
! N.B. Strictly speaking, for this routine, no interpolation back to the hydrogrid is required. 
! The code only needs to knopw the Stromgren radius as a function of latitude.
! This would greatly speed up the simulation.
!
! Also, we could make alpha_h a variable and draw it inside the integral, based on local temperature
! To do this, we have to apply 'ray_primitive'
!
use mod_raytracing

use mod_global_parameters

double precision :: rray(1:npoints),drray, F, dF, Tray
integer          :: iray,ipoint
!
!---------------------------------------------------------------------------
! 
 uraypoint(minrayp(mype):maxrayp(mype),1:npoints,ion_) = zero
 
 call primitive_rays
 
 do iray=minrayp(mype),maxrayp(mype)
    select case( coordinate )
    case( spherical )
      select case( TRIM(raysource) )
      case( 'x1dir' )
          rray(1:npoints) = dabs(xraypoint(iray,1:npoints,1))
      case default 
          call MPISTOP( ' Stromgren radii are always defined in the radial direction. ' )      
      end select
    case default
       call MPISTOP( ' The Stromgren radius concept is based on spherical symmetry. ' )
    end select

    drray = rray(1)

    F    = eqpar(Fstar_)
    Tray = (wraypoint( iray,1,p_ )/wraypoint( iray,1,rho_ )) /eqpar(Tscale_)
    dF   = ((wraypoint(iray,1,rho_)*w_convert_factor(rho_))/mhydro)**2
    dF   = (dF * (4.0d0*(drray**3)* eqpar(alphah_))/3.0D0)   
    F    = F-dF

    if( F < zero) cycle
    
    uraypoint(iray,1,ion_ ) =  one
    do ipoint=2,npointsmax(iray)
     
      Tray  = (wraypoint( iray,ipoint,p_ ) /wraypoint( iray,ipoint,rho_ )) /eqpar(Tscale_)
      drray = rray(ipoint) - rray(ipoint-1)
    
      dF    = half*((rray(ipoint-1)*(wraypoint(iray,ipoint-1,rho_)*w_convert_factor(rho_)/mhydro))**2 &
            + (rray(ipoint)* (wraypoint(iray,ipoint,rho_)*w_convert_factor(rho_)/mhydro))**2)
      dF    = dF * drray * eqpar(alphah_)  
      F     = F-dF
      if( F> zero ) then 
          uraypoint(iray,ipoint,ion_ ) =  one
      else
          rayvar( iray,ion_ ) = rray(ipoint-1)
          exit
      endif
    enddo	    
 enddo
 
 
 return 
 end subroutine ray_stromgren
 
!==============================================================================

subroutine update_rays
!
! Get the ray-data to the relevant processors
!
use mod_raytracing

use mod_global_parameters


double precision :: time2, time1

!---------------------------------------------------------------------------

time1 = MPI_WTIME()

select case(npe) 
 case(1)
!
! No point in calling a comm. routine if there is only one processor
!
 case default
   select case( TRIM(raycomm) )
    case( 'central_single' )
       call  allgather_griddata
!    case( 'central_fast' )
!       call  allgather_griddata_fast
    case( 'central_block' )
       call allgather_griddata_block
    case( 'distribute_single' )
       call get_mygriddata
    case( 'distribute_block' )
       call get_mygriddata_block
    case default
       call MPISTOP( 'This communication method has not been implemented' )
   end select
end select

time2 = MPI_WTIME()
time_ray_comm(mype) = time_ray_comm(mype) + (time2-time1)


time1 = MPI_WTIME()

select case( TRIM(rayphysics) )
 case( 'none' )
!
!    In case you want to trace rays without actually doing anything with them
!
 case(  'columndens' )
    call ray_columndens
 case(  'stromgren' )
!
! Technically, for the Stromgren radius you don't need to return the whole ray-grid
! Just  R_s as a function of Theta should be enough
!
!
    call ray_stromgren
 case( 'starflux' )
    call MPISTOP( 'This will be iomplemented soon' )   
 case default
    call MPISTOP( 'This type of physics has not been implemented' )    
end select

time2 = MPI_WTIME()
time_ray_phys(mype) = time_ray_phys(mype) + (time2-time1)

ptfound( 1:nrays,1:npoints )    = .false.
ptgathered( 1:nrays,1:npoints ) = .false.

time1 = MPI_WTIME()
 
select case(npe) 
 case(1)
!
! No point in calling a comm. routine if there is only one processor
!
 case default
   select case( TRIM(raycomm) )
    case( 'central_single' )
       call allput_raydata
    case( 'central_block' )
       call allput_raydata_block
    case( 'distribute_single' )
       call allput_raydata
    case( 'distribute_block' )
       call allput_raydata_block
    case default
       call MPISTOP( 'This communication method has not been implemented' )
   end select
end select

time2 = MPI_WTIME()
time_ray_comm(mype) = time_ray_comm(mype) + (time2-time1)


return
end subroutine update_rays

!==========================================================================

subroutine allgather_griddata
!
! gather all teh ray-data on teh central node, then send to all processors
!
use mod_raytracing

use mod_global_parameters

integer :: iray,ipoint,ipsend,iprecv
integer :: status(MPI_STATUS_SIZE)

!----------------------------------------------------
 call MPI_BARRIER(icomm,ierrmpi)

 if(mype==0) ptgathered(1:nrays,1:npoints) = ptfound(1:nrays,1:npoints)

 do ipsend = 1,npe-1
   if(mype==0 ) then
      ptfound(1:nrays,1:npoints) = .false.
      call MPI_RECV(ptfound(1:nrays,1:npoints),nrays*npoints,mpi_logical,ipsend,1,icomm,status,ierrmpi )
   else if ( mype==ipsend) then
      call MPI_SEND(ptfound(1:nrays,1:npoints),nrays*npoints,mpi_logical,0,1,icomm,ierrmpi )
   endif   

   if( mype==0 .or. mype==ipsend) then
      do iray=1,nrays
      if( .not.(any(ptfound(iray,1:npoints))) ) cycle
      do ipoint=1,npointsmax(iray)
         if( .not. ptfound(iray,ipoint) )cycle
         if(mype==0 ) then
            call MPI_RECV(wraypoint(iray,ipoint,1:nw),nw,mpi_double_precision,ipsend,2,icomm,status,ierrmpi )
            ptgathered(iray,ipoint) = .true.
         else if ( mype==ipsend) then
            call MPI_SEND(wraypoint(iray,ipoint,1:nw),nw,mpi_double_precision,0,2,icomm,ierrmpi )
         endif
      enddo
      enddo
   endif
 enddo

do iprecv = 1,npe-1
   if(mype==iprecv) then
      call MPI_RECV(wraypoint(minrayp(iprecv):maxrayp(iprecv),1:npoints,1:nw),nrayp(iprecv)*npoints*nw &
                   ,mpi_double_precision,0,4,icomm,status,ierrmpi )
   else if( mype == 0 ) then
      call MPI_SEND(wraypoint(minrayp(iprecv):maxrayp(iprecv),1:npoints,1:nw),nrayp(iprecv)*npoints*nw &
                   ,mpi_double_precision,iprecv,4,icomm,ierrmpi )
   endif
enddo

 call MPI_BARRIER(icomm,ierrmpi)

return
end subroutine allgather_griddata
!==============================================================================

subroutine allgather_griddata_fast
!
! Gather all ray-data on central node, than redistribute, 
! but send per block, rather than as individual points
!
use mod_raytracing

use mod_global_parameters

integer :: iray,ipoint,ipsend,iprecv
integer :: status(MPI_STATUS_SIZE)

integer :: isendmin,isendmax

!----------------------------------------------------

 call MPI_BARRIER(icomm,ierrmpi)

 if(mype==0) ptgathered(1:nrays,1:npoints) = ptfound(1:nrays,1:npoints)

 do ipsend = 1,npe-1
   if(mype==0 ) then
      ptfound(1:nrays,1:npoints) = .false.
      call MPI_RECV(ptfound(1:nrays,1:npoints),nrays*npoints,mpi_logical,ipsend,1,icomm,status,ierrmpi )
   else if ( mype==ipsend) then
      call MPI_SEND(ptfound(1:nrays,1:npoints),nrays*npoints,mpi_logical,0,1,icomm,ierrmpi )
   endif   

   if( mype==0 .or. mype==ipsend) then
      do iray=1,nrays
      isendmin=0; isendmax=0
      if( .not.(any(ptfound(iray,1:npoints))) ) cycle
      do ipoint=1,npointsmax(iray)
          if( isendmin==0 .and. .not.(ptfound(iray,ipoint)) ) cycle
!
!  Establish begin and end of available data, then send it as a single block
!
         if( ptfound(iray,ipoint) .and. (.not.(ptfound(iray,ipoint-1)) .or. ipoint==1) ) then 
	    isendmin = ipoint
 	 endif
         if( isendmin> 0 .and. ((.not.(ptfound(iray,ipoint))) .or. ipoint==npointsmax(iray))  ) then   
	   
	   if(ipoint ==npointsmax(iray) ) then
	       isendmax = ipoint
	   else
               isendmax = ipoint-1 
	   endif
           if(mype==0 ) then
               call MPI_RECV(wraypoint(iray,isendmin:isendmax,1:nw),nw*(1+isendmax-isendmin),mpi_double_precision,ipsend,2,icomm,status,ierrmpi )
               ptgathered(iray,isendmin:isendmax) = .true.
	       isendmin = 0; isendmax=0
            else if ( mype==ipsend) then
               call MPI_SEND(wraypoint(iray,isendmin:isendmax,1:nw),nw*(1+isendmax-isendmin),mpi_double_precision,0,2,icomm,ierrmpi )
	       isendmin=0; isendmax=0
	    endif
	 endif     
      enddo
      enddo

   endif
 enddo

do iprecv = 1,npe-1
   if(mype==iprecv) then
      call MPI_RECV(wraypoint(minrayp(iprecv):maxrayp(iprecv),1:npoints,1:nw),nrayp(iprecv)*npoints*nw &
                   ,mpi_double_precision,0,4,icomm,status,ierrmpi )
   else if( mype == 0 ) then
      call MPI_SEND(wraypoint(minrayp(iprecv):maxrayp(iprecv),1:npoints,1:nw),nrayp(iprecv)*npoints*nw &
                   ,mpi_double_precision,iprecv,4,icomm,ierrmpi )
   endif
enddo

 call MPI_BARRIER(icomm,ierrmpi)


return
end subroutine allgather_griddata_fast
!==============================================================================

subroutine allgather_griddata_block
!
! Collects the data from the processors and 
! redistributes
!
! This routine transfers all available data as a single buffer
!
! Keep in mind that at the end each processor only has data for its own rays. 
! Only the central node has all data


use mod_raytracing

use mod_global_parameters

integer :: iray,ipoint,ipsend,iprecv
integer :: status(MPI_STATUS_SIZE)


double precision,allocatable :: wraysend(:,:)
integer,allocatable          :: pointraysend(:,:)
integer :: ipointsend,npointsend

!----------------------------------------------------



 allocate( wraysend(1:npoints*nrays,1:nw) )
 allocate( pointraysend(1:nrays*npoints,1:2) )

 call MPI_BARRIER(icomm,ierrmpi)


 if(mype==0) ptgathered(1:nrays,1:npoints) = ptfound(1:nrays,1:npoints)


 do ipsend = 1,npe-1
!
! Determine how many points the cpu has
!
   npointsend = 0
   if(mype==0 ) then
         call MPI_RECV( npointsend,1,mpi_integer,ipsend,1,icomm,status,ierrmpi )  
   else if( mype==ipsend) then
      do iray=1,nrays
      if( .not.(any(ptfound(iray,1:npoints))) ) cycle
      do ipoint=1,npointsmax(iray)
         if( ptfound(iray,ipoint)  ) then 
	    npointsend = npointsend+1
!
! Fill buffer with data
!
	    wraysend( npointsend,1:nw )  = wraypoint( iray,ipoint,1:nw )
	    pointraysend( npointsend,1 ) = iray
	    pointraysend( npointsend,2 ) = ipoint
 	 endif
      enddo
      enddo
      call MPI_SEND( npointsend,1,mpi_integer,0,1,icomm,ierrmpi )   
   endif

   if( npointsend == 0 ) cycle  
!
!  Send entire buffer
!	    
   if(mype==0 ) then
        call MPI_RECV( wraysend(1:npointsend,1:nw),nw*npointsend,mpi_double_precision,ipsend,2,icomm,status,ierrmpi )
   else if ( mype==ipsend) then
        call MPI_SEND( wraysend(1:npointsend,1:nw),nw*npointsend,mpi_double_precision,0,2,icomm,ierrmpi )
   endif

   if(mype==0 ) then
        call MPI_RECV( pointraysend(1:npointsend,1:2),2*npointsend,mpi_integer,ipsend,3,icomm,status,ierrmpi )
   else if ( mype==ipsend) then
        call MPI_SEND( pointraysend(1:npointsend,1:2),2*npointsend,mpi_integer,0,3,icomm,ierrmpi )
   endif
!
! Distribute buffer over the grid
!
   if( mype==0 ) then
      do ipointsend=1,npointsend
         iray   = pointraysend(ipointsend,1)
         ipoint = pointraysend(ipointsend,2)
         wraypoint(iray,ipoint,1:nw ) = wraysend(ipointsend,1:nw)
	 ptgathered(iray,ipoint) = .true.
      enddo
   endif
 enddo

do iprecv = 1,npe-1
   if(mype==iprecv) then
      call MPI_RECV(wraypoint(minrayp(iprecv):maxrayp(iprecv),1:npoints,1:nw),nrayp(iprecv)*npoints*nw &
                   ,mpi_double_precision,0,4,icomm,status,ierrmpi )
   else if( mype == 0 ) then
      call MPI_SEND(wraypoint(minrayp(iprecv):maxrayp(iprecv),1:npoints,1:nw),nrayp(iprecv)*npoints*nw &
                   ,mpi_double_precision,iprecv,4,icomm,ierrmpi )
   endif
enddo
 call MPI_BCAST(ptgathered(1:nrays,1:npoints),nrays*npoints,mpi_logical,0,icomm,ierrmpi )

 call MPI_BARRIER(icomm,ierrmpi)

 deallocate( wraysend )
 deallocate( pointraysend )

return
end subroutine allgather_griddata_block
!==============================================================================

subroutine get_mygriddata
!
! Gathers the relevant data for each processor on that processor. 
! Keep in mind that at the end each processor only has data for its own rays. 
!
use mod_raytracing

use mod_global_parameters

integer :: iray,ipoint,ipsend,iprecv
logical, allocatable :: ptfound_send(:,:)

integer :: status(MPI_STATUS_SIZE)

!----------------------------------------------------

allocate( ptfound_send(1:nrays,1:npoints) )

 call MPI_BARRIER(icomm,ierrmpi)

 do iprecv = 0,npe-1
 do ipsend = 0,npe-1
   if(iprecv == ipsend ) cycle
   
   if(mype==iprecv) then
      ptfound_send(1:nrays,1:npoints) = .false.
      call MPI_RECV(ptfound_send(minrayp(iprecv):maxrayp(iprecv),1:npoints),nrayp(iprecv)*npoints,mpi_logical,ipsend,1,icomm,status,ierrmpi )
   else if ( mype==ipsend) then
      ptfound_send(minrayp(iprecv):maxrayp(iprecv),1:npoints) = ptfound(minrayp(iprecv):maxrayp(iprecv),1:npoints)
      call MPI_SEND(ptfound_send(minrayp(iprecv):maxrayp(iprecv),1:npoints),nrayp(iprecv)*npoints,mpi_logical,iprecv,1,icomm,ierrmpi )
   endif   

   do iray=minrayp(iprecv),maxrayp(iprecv)
   if( .not.(any(ptfound_send(iray,1:npoints))) ) cycle
   do ipoint=1,npoints
      if( .not. ptfound_send(iray,ipoint) )cycle      
      if(mype==iprecv ) then
         call MPI_RECV(wraypoint(iray,ipoint,1:nw),nw,mpi_double_precision,ipsend,1,icomm,status,ierrmpi )
         ptgathered(iray,ipoint) = .true.
      else if ( mype==ipsend) then
         call MPI_SEND(wraypoint(iray,ipoint,1:nw),nw,mpi_double_precision,iprecv,1,icomm,ierrmpi )
      endif
   enddo
   enddo

 enddo
 enddo

 call MPI_BARRIER(icomm,ierrmpi)


deallocate( ptfound_send )

return
end subroutine get_mygriddata
!==============================================================================

subroutine get_mygriddata_block
!
! Gets the relevant data to each processor
! This routine transfers all available data as a single buffer
!
! Keep in mind that at the end each processor only has data for its own rays. 
!

use mod_raytracing

use mod_global_parameters

integer :: iray,ipoint
integer :: status(MPI_STATUS_SIZE)


double precision,allocatable :: wraysend(:,:)
integer,allocatable          :: pointraysend(:,:)
integer                      :: ipointsend,npointsend

integer                      :: ipsend,iprecv

!----------------------------------------------------

 allocate( wraysend(1:npoints*nrays,1:nw) )
 allocate( pointraysend(1:nrays*npoints,1:2) )

 call MPI_BARRIER(icomm,ierrmpi)


 do iprecv = 0,npe-1
 do ipsend = 0,npe-1
   if( iprecv==ipsend ) cycle
!
! Determine how relevantmany points the sending cpu has
!
   npointsend = 0
   if(mype==iprecv ) then
         call MPI_RECV( npointsend,1,mpi_integer,ipsend,1,icomm,status,ierrmpi )  
   else if( mype==ipsend) then
      do iray=minrayp(iprecv),maxrayp(iprecv)
      if( .not.(any(ptfound(iray,1:npoints))) ) cycle
      do ipoint=1,npointsmax(iray)
         if( ptfound(iray,ipoint)  ) then 
	    npointsend = npointsend+1
!
! Fill buffer with data
!
	    wraysend( npointsend,1:nw )  = wraypoint( iray,ipoint,1:nw )
	    pointraysend( npointsend,1 ) = iray
	    pointraysend( npointsend,2 ) = ipoint
 	 endif
      enddo
      enddo
      call MPI_SEND( npointsend,1,mpi_integer,iprecv,1,icomm,ierrmpi )   
   endif

   if( npointsend == 0 ) cycle  
!
!  Send entire buffer
!	    

   if(mype==iprecv ) then
        call MPI_RECV( wraysend(1:npointsend,1:nw),nw*npointsend,mpi_double_precision,ipsend,2,icomm,status,ierrmpi )
   else if ( mype==ipsend) then
        call MPI_SEND( wraysend(1:npointsend,1:nw),nw*npointsend,mpi_double_precision,iprecv,2,icomm,ierrmpi )
   endif

   if(mype==iprecv ) then
        call MPI_RECV( pointraysend(1:npointsend,1:2),2*npointsend,mpi_integer,ipsend,3,icomm,status,ierrmpi )
   else if ( mype==ipsend) then
        call MPI_SEND( pointraysend(1:npointsend,1:2),2*npointsend,mpi_integer,iprecv,3,icomm,ierrmpi )
   endif
!
! Distribute buffer over the grid
!
   if( mype==iprecv ) then
      do ipointsend=1,npointsend
         iray   = pointraysend(ipointsend,1)
         ipoint = pointraysend(ipointsend,2)
         wraypoint(iray,ipoint,1:nw ) = wraysend(ipointsend,1:nw)
	 ptgathered(iray,ipoint) = .true.
      enddo
   endif
 enddo
 enddo

 call MPI_BARRIER(icomm,ierrmpi)

 deallocate( wraysend )
 deallocate( pointraysend )

return
end subroutine get_mygriddata_block
!==========================================================================

subroutine allput_raydata
!
! Send all ray-specific variables to central node, then redistribute
!
use mod_raytracing

use mod_global_parameters

integer :: iray,ipoint,ipsend
integer :: status(MPI_STATUS_SIZE)

!----------------------------------------------------
 call MPI_BARRIER(icomm,ierrmpi)

 if(mype==0) ptgathered(1:nrays,1:npoints) = ptfound(1:nrays,1:npoints)
!
! Send all data to central node
!
 do ipsend = 1,npe-1
   if( mype==0 .or. mype==ipsend) then
      do iray=minrayp(ipsend),maxrayp(ipsend)
      do ipoint=1,npointsmax(iray)
         if(mype==0 ) then
            call MPI_RECV(uraypoint(iray,ipoint,1:nu),nu,mpi_double_precision,ipsend,1,icomm,status,ierrmpi )
            ptgathered(iray,ipoint) = .true.
         else if ( mype==ipsend) then
            call MPI_SEND(uraypoint(iray,ipoint,1:nu),nu,mpi_double_precision,0,1,icomm,ierrmpi )
         endif
      enddo
         if(mype==0 ) then
            call MPI_RECV(rayvar(iray,1:nu)    &
	                 ,1,mpi_double_precision,ipsend,1,icomm,status,ierrmpi )
         else if ( mype==ipsend) then
            call MPI_SEND(rayvar(iray,1:nu)    &
	                 ,1,mpi_double_precision,0,1,icomm,ierrmpi )
         endif

      enddo
   endif
 enddo

!
! Broadcast to all processors
!

 call MPI_BCAST(uraypoint(1:nrays,1:npoints,1:nu),nrays*npoints*nu,mpi_double_precision,0,icomm,ierrmpi )

 call MPI_BARRIER(icomm,ierrmpi)

return
end subroutine allput_raydata
!==============================================================================
subroutine allput_raydata_block
!
! Send all ray-specific variables to central node, then redistribute
!
! moves data per block
!
use mod_raytracing

use mod_global_parameters

integer :: iray,ipoint,ipsend
integer :: status(MPI_STATUS_SIZE)

!----------------------------------------------------
 call MPI_BARRIER(icomm,ierrmpi)

 if(mype==0) ptgathered(1:nrays,1:npoints) = ptfound(1:nrays,1:npoints)

 do ipsend = 1,npe-1
   if( mype==0 .or. mype==ipsend) then
         if(mype==0 ) then
            call MPI_RECV(uraypoint(minrayp(ipsend):maxrayp(ipsend),1:npoints,1:nu)    &
	                 ,nu*nrayp(ipsend)*npoints  &
	                 ,mpi_double_precision,ipsend,1,icomm,status,ierrmpi )
         else if ( mype==ipsend) then
            call MPI_SEND(uraypoint(minrayp(ipsend):maxrayp(ipsend),1:npoints,1:nu)    &
	                 ,nu*nrayp(ipsend)*npoints  &
	                 ,mpi_double_precision,0,1,icomm,ierrmpi )
         endif

         if(mype==0 ) then
            call MPI_RECV(rayvar(minrayp(ipsend):maxrayp(ipsend),1:nu)    &
	                 ,nu*nrayp(ipsend),mpi_double_precision,ipsend,1,icomm,status,ierrmpi )
         else if ( mype==ipsend) then
            call MPI_SEND(rayvar(minrayp(ipsend):maxrayp(ipsend),1:nu)    &
	                 ,nu*nrayp(ipsend),mpi_double_precision,0,1,icomm,ierrmpi )
         endif

   endif
 enddo

 call MPI_BCAST(uraypoint(1:nrays,1:npoints,1:nu),nrays*npoints*nu,mpi_double_precision,0,icomm,ierrmpi )


 call MPI_BARRIER(icomm,ierrmpi)

return
end subroutine allput_raydata_block
!==============================================================================
subroutine write_ray_data( outname )
!
!  Get all data from the rays and write them to a file
!
 use mod_raytracing

 use mod_global_parameters

 character*79, intent(in) :: outname

 integer :: iray,ipoint

!----------------------------------------

 if( .not. (writerays) ) return

select case( TRIM(raycomm) )
 case( 'central_block' )
    call allgather_griddata_block
    call allput_raydata_block
 case( 'cenral_single' )
    call allgather_griddata_block
    call allput_raydata
 case( 'distribute_block' )
    call allgather_griddata_block
    call allput_raydata_block
 case( 'distribute_single' )
    call allgather_griddata_block
    call allput_raydata
 case default
    call MPISTOP( 'This communication method has not been implemented' )
end select 
 
 
 call primitive_rays
 
 if( mype==0) then
    open(11,file=TRIM(outname), status='unknown')
    write(11,*) 'TITLE="RAYS"'
{^IFONED   write(11,*) 'VARIABLES = X,RHO,V1,V2,P,U1' }
{^IFTWOD    write(11,*) 'VARIABLES = X,Y,RHO,V1,V2,P,U1' }
!!!{^IFTWOD    write(11,*) 'VARIABLES = X,Y,Z,RHO,V1,V2,V3,P,U1' }
   
    do iray = 1,nrays
       write(11, '(a,i4,a,i4a,1pe12.5,a)' ) 'ZONE T= "',iray,'",I=',npointsmax(iray) ,',SOLUTIONTIME=',global_time*time_convert_factor,',F=POINT'
       do ipoint = 1,npointsmax(iray)
          write(11,1001) {^D&xraypoint(iray,ipoint,^D)},wraypoint(iray,ipoint,1:nw)*w_convert_factor(1:nw) &
	                 ,uraypoint(iray,ipoint,1:nu)
       enddo
       write(11,*)
    enddo
    close(11)
 endif
 
 call conservative_rays
 
return
1001 format(20(1x,1pe12.5))
end subroutine write_ray_data
!==============================================================================

subroutine primitive_rays
!
! Transform from conservative to primitive
! For the moment only for hd and mhd. (mhd version not yet tested)

 use mod_raytracing

 use mod_global_parameters 
 
 integer :: iray
!-------------------------------------------------------------------------------

do iray = 1,nrays
{^C&wraypoint(iray,1:npointsmax(iray),v^C_)=wraypoint(iray,1:npointsmax(iray),m^C_)/wraypoint(iray,1:npointsmax(iray),rho_);}

wraypoint(iray,1:npointsmax(iray),p_)=(eqpar(gamma_)-one) * (wraypoint(iray,1:npointsmax(iray),e_) &
                -half*wraypoint(iray,1:npointsmax(iray),rho_)  &
		*{^C&(wraypoint(iray,1:npointsmax(iray),v^C_)**2)+} )
{^IFMHD
wraypoint(iray,1:npointsmax(iray),p_)=wraypoint(iray,1:npointsmax(iray),p_)  &
                -(eqpar(gamma_)-one) &
	        *half*({^C&wraypoint(iray,1:npointsmax(iray),b^C_)**2.0d0+})
}
enddo


return
end subroutine primitive_rays
!==============================================================================

subroutine conservative_rays
!
! Transform from  primitive to conservative
! For the moment only for hd and mhd.  (mhd version not yet tested)

 use mod_raytracing

 use mod_global_parameters
  
 integer :: iray
!-------------------------------------------------------------------------------


do iray = 1,nrays
{^C&wraypoint(iray,1:npointsmax(iray),m^C_)=wraypoint(iray,1:npointsmax(iray),m^C_)*wraypoint(iray,1:npointsmax(iray),rho_);}

wraypoint(iray,1:npointsmax(iray),e_)=wraypoint(iray,1:npointsmax(iray),p_)/(eqpar(gamma_)-one)     &
                + half *({^C&(wraypoint(iray,1:npointsmax(iray),v^C_)**2)+}) &
		/wraypoint(iray,1:npointsmax(iray),rho_) 

{^IFMHD	
wraypoint(iray,1:npointsmax(iray),e_)=wraypoint(iray,1:npointsmax(iray),e_)  &
                + half*({^C&wraypoint(iray,1:npointsmax(iray),b^C_)**2.0d0+})
}
enddo

return
end subroutine conservative_rays

!==============================================================================

subroutine getraydata( ixI^L,ixO^L,x,u )
!
!  Get ray data to the grid
!
use mod_raytracing

use mod_global_parameters

integer, intent(in)           :: ixI^L, ixO^L
double precision, intent(in)  :: x(ixI^S,1:ndim)
double precision, intent(out) :: u(ixG^T,1:nu)


logical           :: patchray(1:nrays)
double precision  :: xpoint^D
double precision  :: rtoray(1:nrays),rmin          
integer           :: ix^D, iray,ipoint,iraynear,ipointnear
integer           :: ir1,ip1,ir2,ip2

double precision :: time2, time1

!---------------------------------------------------

time1 = MPI_WTIME()

patchray(1:nrays) = .true.

if ( checkintersect ) then 
   call intersect_ray_grid( ixI^L,ixO^L,x,patchray )
   if( .not.(any(patchray(1:nrays))) ) then
      call MPISTOP( 'Your radiative grid is too course' )
   endif
endif


{do ix^D = ixO^LIM^D\}
            xpoint1 = x(ix^D,1)
{^IFTWOD    xpoint2 = x(ix^D,2)}
{^IFTHREED  xpoint3 = x(ix^D,3)}

   if( regular_rays ) then 
      call findnearestraypoint_regular( xpoint^D, patchray, iray,ipoint )
      u(ix^D,1:nu ) = uraypoint(iray,ipoint,1:nu)
   else
      select case( raytogrid )
       case( 'nearest' )
          call findnearestraypoint( xpoint^D, patchray, iray,ipoint )
          u(ix^D,1:nu ) = uraypoint(iray,ipoint,1:nu)
       case(  'nearestray' )
          call point_to_ray( xpoint^D, patchray,rtoray )
          rmin = 1.0d32
          do iray=1,nrays
             if( .not. patchray(iray) ) cycle
	     if( rtoray(iray) < smalldouble ) then
	        iraynear = iray   
	        exit
	     else if( rtoray(iray ) < rmin ) then
               rmin     = rtoray(iray)
	       iraynear = iray   
            endif  
          enddo
          call findnearestonray( xpoint^D, iraynear,ipointnear )
          u(ix^D,1:nu ) = uraypoint(iraynear,ipointnear,1:nu)
       case( 'raytogrid_usr' )
           call MPISTOP( 'You have to define this for yourself' )    
       case default
          call MPISTOP( 'This method has not been implemented' )
       end select
   endif

{enddo^D&\}	    

time2 = MPI_WTIME()

time_ray_interpol(mype) = time_ray_interpol(mype) + (time2-time1)

return
end subroutine getraydata
!============================================================================

subroutine intersect_ray_grid( ixI^L,ixO^L,x,patchray )
!
! Find which rays come close o intersecting with the grid 
!
use mod_raytracing

use mod_global_parameters

integer, intent(in)          :: ixI^L, ixO^L
double precision, intent(in) :: x(ixI^S,1:ndim)

logical, intent(inout) :: patchray( 1:nrays )

double precision :: xmin^D,xmax^D,dx^D
integer :: iray

!---------------------------------------------------------------------------

!
! Exclude those rays that are more than a full grid-width away from the grid border
!
{^D&xmin^D=minval(x(ixO^S,^D));}
{^D&xmax^D=maxval(x(ixO^S,^D));}
{^D&dx^D=0.5d0*dabs(xmax^D-xmin^D);}
{^D&xmin^D=xmin^D-dx^D;}
{^D&xmax^D=xmax^D+dx^D;}

 do iray = 1,nrays

{^IFTWOD
       if( xraymin(iray,1) > xmax1 .or. xraymax(iray,1)<xmin1 ) patchray(iray) = .false.
       if( xraymin(iray,2) > xmax2 .or. xraymax(iray,2)<xmin2 ) patchray(iray) = .false.
}
{^IFTHREED
           xraymin(iray,3) > xmax3 .or. xraymax(iray,3)<xmin3 ) patchray(iray) = .false.
}
 enddo

return
end subroutine intersect_ray_grid

!============================================================================
subroutine findnearestraypoint(  xpoint^D,patchray,irayout,ipointout )
!
! Find the indices of the nearest ray-point
!
use mod_raytracing

use mod_global_parameters

double precision, intent(in) :: xpoint^D
logical, intent(in)          :: patchray(1:nrays)
integer, intent(out)         :: irayout,ipointout

double precision             :: r2min,dist2(1:nrays,1:npoints)  !,dist(1:nrays,1:npoints)
{^IFONED   double precision  :: r2topoint(1:2)}
{^IFTWOD   double precision,allocatable  :: r2topoint(:,:)}

integer                      :: iray,ipoint 
integer                      :: ilpoint,ihpoint,ilray,ihray

!--------------------------------------------------------------------------

irayout   = 1
ipointout = 1
r2min=1.0d32   


{^IFONED
do ipoint=1,npointsmax(1)
   if( xraypoint(1,ipoint,1) >= xpoint1 ) then 
      ilpoint =ipoint-1;ihpoint=ipoint
      exit
   endif
enddo
r2topoint(1) = (xraypoint(1,ilpoint,1)-xpoint1)**2
r2topoint(2) = (xraypoint(1,ihpoint,1)-xpoint1)**2
if( r2topoint(1)< r2topoint(2) ) then
  irayout   = 1
  ipointout = ilpoint
else
  irayout   = 1
  ipointout = ihpoint   
endif
}


{^IFTWOD
dist2(1:nrays,1:npoints) = {^D&((xraypoint(1:nrays,1:npoints,^D)-xpoint^D)**2)+}
   
do iray=1,nrays
   if( .not.(patchray(iray)) ) cycle
   do ipoint = 1,npointsmax(iray)
      if( dist2(iray,ipoint) <= smalldouble) then
         ipointout = ipoint
         irayout   = iray
	 return
      else if( dist2(iray,ipoint) < r2min ) then
         r2min     = dist2(iray,ipoint)
	 irayout   = iray
	 ipointout = ipoint
      endif
   enddo
enddo
}


{^IFTHREED
dist2(1:nrays,1:npoints) = {^D&((xraypoint(1:nrays,1:npoints,^D)-xpoint^D)**2)+} 

do iray=1,nrays
   if( .not.(patchray(iray)) ) cycle
   do ipoint = 1,npointsmax(iray)
       if( dist2(iray,ipoint) <= smalldouble) then
         ipointout = ipoint
         irayout  = iray
	 return
       else if( dist2(iray,ipoint) < rmin ) then
         rmin      = dist2(iray,ipoint)
	 irayout   = iray
	 ipointout = ipoint
       endif
   enddo
enddo

}


return
end subroutine findnearestraypoint 
!================================================================================
subroutine findnearestraypoint_regular( xpoint^D,patchray,irayout,ipointout )
!
! interpolation scheme for 'regular' ray grid
!
use mod_raytracing
use mod_global_parameters

double precision, intent(in) :: xpoint^D
logical, intent(in)          :: patchray(1:nrays)
integer, intent(out)         :: irayout,ipointout

double precision             :: r2min,dist2(1:nrays,1:npoints)  !,dist(1:nrays,1:npoints)
{^IFONED   double precision  :: r2topoint(1:2)}
{^IFTWOD   double precision,allocatable  :: r2topoint(:,:)}

integer                      :: iray,ipoint 
integer                      :: ilpoint,ihpoint,ilray,ihray

!-----------------------------------------------------------------------------------
irayout   = 1
ipointout = 1
r2min=1.0d32   


{^IFONED
do ipoint=1,npointsmax(1)
   if( xraypoint(1,ipoint,1) >= xpoint1 ) then 
      ilpoint =ipoint-1;ihpoint=ipoint
      exit
   endif
enddo
r2topoint(1) = (xraypoint(1,ilpoint,1)-xpoint1)**2
r2topoint(2) = (xraypoint(1,ihpoint,1)-xpoint1)**2
if( r2topoint(1)< r2topoint(2) ) then
  irayout   = 1
  ipointout = ilpoint
else
  irayout   = 1
  ipointout = ihpoint   
endif
}

{^IFTWOD
if( TRIM(raysource)== 'x1dir' ) then 
   do iray=1,nrays
      if( .not. (patchray(iray)) ) cycle
      if( xraypoint(iray,1,2) >= xpoint2 ) then 
         ilray =iray-1;ihray=iray
         exit
      endif
   enddo
   do ipoint=1,npointsmax(1)
      if( xraypoint(1,ipoint,1) >= xpoint1 ) then 
         ilpoint =ipoint-1;ihpoint=ipoint
         exit
      endif
   enddo
else if( TRIM(raysource) =='x2dir' ) then
   do iray=1,nrays
      if( .not. (patchray(iray)) ) cycle
      if( xraypoint(iray,1,1) >= xpoint1 ) then 
         ilray =iray-1;ihray=iray
         exit
      endif
   enddo
   do ipoint=1,npointsmax(1)
      if( xraypoint(1,ipoint,2) >= xpoint2 ) then 
         ilpoint =ipoint-1;ihpoint=ipoint
         exit
      endif
   enddo
else
   call MPISTOP( 'This is not a regular ray_grid' ) 
endif
   
allocate( r2topoint(ilray:ihray,ilpoint:ihpoint) )
r2topoint(ilray,ilpoint ) = (xraypoint(ilray,ilpoint,1)-xpoint1)**2 + (xraypoint(ilray,ilpoint,2)-xpoint2)**2
r2topoint(ihray,ilpoint ) = (xraypoint(ihray,ilpoint,1)-xpoint1)**2 + (xraypoint(ihray,ilpoint,2)-xpoint2)**2
r2topoint(ilray,ihpoint ) = (xraypoint(ilray,ihpoint,1)-xpoint1)**2 + (xraypoint(ilray,ihpoint,2)-xpoint2)**2
r2topoint(ihray,ihpoint ) = (xraypoint(ihray,ihpoint,1)-xpoint1)**2 + (xraypoint(ihray,ihpoint,2)-xpoint2)**2
do iray = ilray,ihray
do ipoint = ilpoint,ihpoint
   if( r2topoint(iray,ipoint) <=smalldouble ) then
      ipointout = ipoint
      irayout   = iray
      exit
   else if( r2topoint(iray,ipoint) < r2min ) then
      r2min = r2topoint(iray,ipoint)
      ipointout = ipoint
      irayout   = iray
   endif
enddo
enddo
deallocate( r2topoint )
}


{^IFTHREED
  call MPISTOP( 'This is not implemented in 3-D' ) 
}




return
end subroutine findnearestraypoint_regular
!============================================================================
subroutine findnearestonray(  xpoint^D,irayin,ipointout   )
!
! Find the indices of the nearest point on a given ray
!
use mod_raytracing

use mod_global_parameters

double precision, intent(in) :: xpoint^D
integer, intent(in)          :: irayin
integer, intent(out)         :: ipointout

double precision  :: r2min,dist2(1:npoints)
integer           :: ipoint 

!--------------------------------------------------------------------------

ipointout = 1
r2min     = 1.0d32   

dist2(1:npoints) = {^D&((xraypoint(irayin,1:npoints,^D)-xpoint^D)**2)+} 

do ipoint = 1,npointsmax(irayin)
   if( dist2(ipoint) < r2min ) then
      r2min = dist2(ipoint)
      ipointout = ipoint
   else
      exit
   endif
enddo
 
return
end subroutine findnearestonray

!=============================================================================

subroutine point_to_ray(xpoint^D, patchray,rtoray )
!
! returns for a given point the shortest distance to the rays
!
use mod_raytracing
use mod_global_parameters

double precision,intent(in)  :: xpoint^D
logical,intent(in)           :: patchray(1:nrays)


double precision,intent(out) :: rtoray(1:nrays)


integer :: iray
!-----------------------------------------------------------------

rtoray(1:nrays) = 1.0d32
{^IFONED
   rtoray(1) = zero
}

{^IFTWOD
do iray = 1,nrays
   if( .not.(patchray(iray)) ) cycle
   rtoray(iray) = dabs(((xrayend(iray,1)-xraystart(iray,1))*(xraystart(iray,2)-xpoint2) &
                - (xraystart(iray,1)-xpoint1)*(xrayend(iray,2)-xraystart(iray,2))) &
	        / dsqrt((xrayend(iray,1)-xraystart(iray,1))**2 + (xrayend(iray,2)-xraystart(iray,2))**2))
enddo
}

return
end subroutine point_to_ray

!=============================================================================

subroutine time_spent_on_rays
!
! Output of time spent on raytracing
!
use mod_raytracing
use mod_global_parameters

double precision :: trinit,trcomm,trinterpol,trphys

!---------------------------------------------------------------------------

 call MPI_REDUCE(time_ray_init(mype),trinit,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
 call MPI_REDUCE(time_ray_comm(mype),trcomm,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
 call MPI_REDUCE(time_ray_interpol(mype),trinterpol,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
 call MPI_REDUCE(time_ray_phys(mype),trphys,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)

if( mype ==0 ) then
   write(*,'(a,f12.3,a)')'Ray initialization took    : ',trinit/npe,' sec'
   write(*,'(a,f12.3,a)')'Ray communication took     : ',trcomm/npe,' sec'
   write(*,'(a,f12.3,a)')'Ray interpolation took     : ',trinterpol/npe,' sec'
   write(*,'(a,f12.3,a)')'Ray physics took           : ',trphys/npe,' sec'
end if

end subroutine time_spent_on_rays



!=============================================================================
!
! To implement the raytracing, you have toput a line in process_grid_usr:
!
! call  getgriddata(ixI^L,ixO^L,x,w)
!
! and in special_source (or wherever you need the ray data:
!
!   call getraydata( ixI^L,ixO^L,x,u) 
!
!  with u(ixG^T) a locally defined variable.
!
! You also have to replace the amrvac.t, advance.t and amrio.t 
! files with the new files that are part of this add-on.
!

