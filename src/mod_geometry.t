!> Module with geometry-related routines (e.g., divergence, curl)
module mod_geometry
  use mod_comm_lib, only: mpistop
  implicit none
  public

  integer :: coordinate=-1
  !$acc declare copyin(coordinate)
  integer, parameter :: Cartesian          = 0
  integer, parameter :: Cartesian_stretched= 1
  integer, parameter :: cylindrical        = 2
  integer, parameter :: spherical          = 3
  integer, parameter :: Cartesian_expansion= 4

  integer :: type_curl=0
  !$acc declare copyin(type_curl)
  integer, parameter :: central=1
  integer, parameter :: Gaussbased=2
  integer, parameter :: Stokesbased=3

contains

  !> Set the coordinate system to be used
  subroutine set_coordinate_system(geom)
    use mod_global_parameters
    character(len=*), intent(in) :: geom !< Name of the coordinate system

    ! Store the geometry name
    geometry_name = geom

    select case (geom)
    case ("Cartesian","Cartesian_1D","Cartesian_2D","Cartesian_3D")
      ndir = ndim
      coordinate=Cartesian
    case ("Cartesian_1D_expansion")
      if (ndim /= 1) call mpistop("Geometry Cartesian_1D_expansion but ndim /= 1")
      ndir = ndim
      coordinate=Cartesian_expansion
    case ("Cartesian_1.5D")
      if (ndim /= 1) call mpistop("Geometry Cartesian_1.5D but ndim /= 1")
      coordinate=Cartesian
      ndir = 2
    case ("Cartesian_1.75D")
      if (ndim /= 1) call mpistop("Geometry Cartesian_1.75D but ndim /= 1")
      coordinate=Cartesian
      ndir = 3
    case ("Cartesian_2.5D")
      if (ndim /= 2) call mpistop("Geometry Cartesian_2.5D but ndim /= 2")
      coordinate=Cartesian
      ndir = 3
    case ("cylindrical","cylindrical_2D","cylindrical_3D")
      ndir = ndim
      r_   = 1
      z_   = 2
      if(ndir==3) phi_ = 3
      coordinate=cylindrical
    case ("cylindrical_2.5D")
      if (ndim /= 2) call mpistop("Geometry cylindrical_2.5D but ndim /= 2")
      ndir = 3
      r_   = 1
      z_   = 2
      phi_ = 3
      coordinate=cylindrical
    case ("polar","polar_2D","polar_3D")
      ndir = ndim
      r_   = 1
      phi_ = 2
      if(ndir==3) z_ = 3
      coordinate=cylindrical
    case ("polar_1.5D")
       if (ndim /= 1) call mpistop("Geometry polar_1.5D but ndim /= 1")
       ndir = 2
       r_   = 1
       phi_ = 2
       coordinate=cylindrical
    case ("polar_2.5D")
      if (ndim /= 2) call mpistop("Geometry polar_2.5D but ndim /= 2")
      ndir = 3
      r_   = 1
      phi_ = 2
      z_   = 3
      coordinate=cylindrical
    case ("spherical","spherical_2D","spherical_3D")
      ndir = ndim
      r_   = 1
      if(ndir==3) phi_ = 3
      z_   = -1
      coordinate=spherical
    case ("spherical_2.5D")
      if (ndim /= 2) &
           call mpistop("Geometry spherical_2.5D requires ndim == 2")
      ndir = 3
      r_   = 1
      phi_ = 3
      z_   = -1
      coordinate=spherical
    case default
      call mpistop("Unknown geometry specified")
    end select
  end subroutine set_coordinate_system

  subroutine set_pole
    use mod_global_parameters

    select case (coordinate)
    case (spherical) {^IFTHREED
      ! For spherical grid, check whether phi-direction is periodic
      if(periodB(ndim)) then
        if(phi_/=3) call mpistop("phi_ should be 3 in 3D spherical coord!")
        if(mod(ng3(1),2)/=0) &
             call mpistop("Number of meshes in phi-direction should be even!")
        if(abs(xprobmin2)<smalldouble) then
          if(mype==0) write(unitterm,*) &
               "Will apply pi-periodic conditions at northpole!"
          poleB(1,2)=.true.
        else
          if(mype==0) write(unitterm,*) "There is no northpole!"
        end if
        if(abs(xprobmax2-dpi)<smalldouble) then
          if(mype==0) write(unitterm,*) &
               "Will apply pi-periodic conditions at southpole!"
          poleB(2,2)=.true.
        else
          if(mype==0) write(unitterm,*) "There is no southpole!"
        end if
      end if}
    case (cylindrical)
      {
      if (^D == phi_ .and. periodB(^D)) then
        if(mod(ng^D(1),2)/=0) then
          call mpistop("Number of meshes in phi-direction should be even!")
        end if

        if(abs(xprobmin1)<smalldouble) then
          if (mype==0) then
            write(unitterm,*) "Will apply pi-periodic conditions at r=0"
          end if
          poleB(1,1)=.true.
        else
          if (mype==0) then
            write(unitterm,*) "There is no cylindrical axis!"
          end if
        end if
      end if\}
    end select

  end subroutine set_pole

  !> Deallocate geometry-related variables
  subroutine putgridgeo(igrid)
    use mod_global_parameters
    integer, intent(in) :: igrid

    deallocate(ps(igrid)%surfaceC,ps(igrid)%surface,ps(igrid)%dvolume,ps(igrid)%dx,&
         psc(igrid)%dx,ps(igrid)%ds,psc(igrid)%dvolume)

  end subroutine putgridgeo

  !> calculate area of surfaces of cells
  subroutine get_surface_area(s,ixG^L)
    use mod_global_parameters
    use mod_usr_methods, only: usr_set_surface

    type(state) :: s
    integer, intent(in) :: ixG^L

    double precision :: x(ixG^S,ndim), xext(ixGmin1-1:ixGmax1,1), drs(ixG^S), drs_ext(ixGmin1-1:ixGmax1), dx2(ixG^S), dx3(ixG^S)
    double precision :: exp_factor_ext(ixGmin1-1:ixGmax1),del_exp_factor_ext(ixGmin1-1:ixGmax1),exp_factor_primitive_ext(ixGmin1-1:ixGmax1)
    double precision :: exp_factor(ixGmin1:ixGmax1),del_exp_factor(ixGmin1:ixGmax1),exp_factor_primitive(ixGmin1:ixGmax1)

    select case (coordinate)

    case (Cartesian_expansion)
      drs(ixG^S)=s%dx(ixG^S,1)
      x(ixG^S,1)=s%x(ixG^S,1)
      {^IFONED
      if(associated(usr_set_surface))then
           call usr_set_surface(ixGmin1,ixGmax1,x,drs,exp_factor,del_exp_factor,exp_factor_primitive)
           if (any(exp_factor <= zero)) call mpistop("The area must always be strictly positive!")
      endif
      s%surface(ixG^S,1)=exp_factor(ixG^S)
      xext(ixGmin1-1,1)=x(1,1)-half*drs(1)
      xext(ixG^S,1)=x(ixG^S,1)+half*drs(ixG^S)
      drs_ext(ixGmin1-1)=drs(1)
      drs_ext(ixG^S)=drs(ixG^S)
      if(associated(usr_set_surface)) call usr_set_surface(ixGmin1-1,ixGmax1,xext,drs_ext,exp_factor_ext,del_exp_factor_ext,exp_factor_primitive_ext)
      s%surfaceC(ixGmin1-1:ixGmax1,1)=exp_factor_ext(ixGmin1-1:ixGmax1)
      }

    case (Cartesian,Cartesian_stretched)
      drs(ixG^S)=s%dx(ixG^S,1)
      {^NOONED
      dx2(ixG^S)=s%dx(ixG^S,2)}
      {^IFTHREED
      dx3(ixG^S)=s%dx(ixG^S,3)}

      {^IFONED
      s%surfaceC(ixG^S,1)=1.d0
      s%surface(ixG^S,1) =1.d0
      }
      {^IFTWOD
      s%surfaceC(ixG^S,1)=dx2(ixG^S)
      s%surfaceC(ixG^S,2)=drs(ixG^S)
      s%surface(ixG^S,1) =dx2(ixG^S)
      s%surface(ixG^S,2)=drs(ixG^S)
      }
      {^IFTHREED
      s%surfaceC(ixG^S,1)= dx2(ixG^S)*dx3(ixG^S)
      s%surfaceC(ixG^S,2)= drs(ixG^S)*dx3(ixG^S)
      s%surfaceC(ixG^S,3)= drs(ixG^S)*dx2(ixG^S)
      s%surface(ixG^S,1)=s%surfaceC(ixG^S,1)
      s%surface(ixG^S,2)=s%surfaceC(ixG^S,2)
      s%surface(ixG^S,3)=s%surfaceC(ixG^S,3)
      }
      {s%surfaceC(0^D%ixG^S,^D)=s%surfaceC(1^D%ixG^S,^D);\}
    case (spherical)
      x(ixG^S,1)=s%x(ixG^S,1)
      {^NOONED
      x(ixG^S,2)=s%x(ixG^S,2)
      }
      drs(ixG^S)=s%dx(ixG^S,1)
      {^NOONED
      dx2(ixG^S)=s%dx(ixG^S,2)
      }
      {^IFTHREED
      dx3(ixG^S)=s%dx(ixG^S,3)
      }

      s%surfaceC(ixG^S,1)=(x(ixG^S,1)+half*drs(ixG^S))**2 {^NOONED &
           *two*dsin(x(ixG^S,2))*dsin(half*dx2(ixG^S))}{^IFTHREED*dx3(ixG^S)}

      {^NOONED
      s%surfaceC(ixG^S,2)=x(ixG^S,1)*drs(ixG^S)&
           *dsin(x(ixG^S,2)+half*dx2(ixG^S))}{^IFTHREED*dx3(ixG^S)}

      {^IFTHREED
      s%surfaceC(ixG^S,3)=x(ixG^S,1)*drs(ixG^S)*dx2(ixG^S)
      }

      {^IFONED
      s%surfaceC(0,1)=dabs(x(1,1)-half*drs(1))**2
      }
      {^IFTWOD
      s%surfaceC(0,ixGmin2:ixGmax2,1)=(x(1,ixGmin2:ixGmax2,1)-half*drs(1,&
         ixGmin2:ixGmax2))**2*two*dsin(x(1,ixGmin2:ixGmax2,2))*dsin(half*dx2(1,&
         ixGmin2:ixGmax2))
      s%surfaceC(ixGmin1:ixGmax1,0,2)=x(ixGmin1:ixGmax1,1,&
         1)*drs(ixGmin1:ixGmax1,1)*dsin(x(ixGmin1:ixGmax1,1,&
         2)-half*dx2(ixGmin1:ixGmax1,1))
      }
      {^IFTHREED
      s%surfaceC(0,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1)=(x(1,ixGmin2:ixGmax2,&
         ixGmin3:ixGmax3,1)-half*drs(1,ixGmin2:ixGmax2,&
         ixGmin3:ixGmax3))**2*two*dsin(x(1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
         2))*dsin(half*dx2(1,ixGmin2:ixGmax2,ixGmin3:ixGmax3))*dx3(1,&
         ixGmin2:ixGmax2,ixGmin3:ixGmax3)
      s%surfaceC(ixGmin1:ixGmax1,0,ixGmin3:ixGmax3,2)=x(ixGmin1:ixGmax1,1,&
         ixGmin3:ixGmax3,1)*drs(ixGmin1:ixGmax1,1,&
         ixGmin3:ixGmax3)*dsin(x(ixGmin1:ixGmax1,1,ixGmin3:ixGmax3,&
         2)-half*dx2(ixGmin1:ixGmax1,1,ixGmin3:ixGmax3))*dx3(ixGmin1:ixGmax1,1,&
         ixGmin3:ixGmax3)
      s%surfaceC(ixGmin1:ixGmax1,ixGmin2:ixGmax2,0,3)=&
         s%surfaceC(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1,3)
      }

      s%surface(ixG^S,1)=x(ixG^S,1)**2 {^NOONED &
           *two*dsin(x(ixG^S,2))*dsin(half*dx2(ixG^S))}{^IFTHREED*dx3(ixG^S)}
      {^NOONED
      s%surface(ixG^S,2)=x(ixG^S,1)*drs(ixG^S)&
           *dsin(x(ixG^S,2))}{^IFTHREED*dx3(ixG^S)}

      {^IFTHREED
      s%surface(ixG^S,3)=x(ixG^S,1)*drs(ixG^S)*dx2(ixG^S)}

    case (cylindrical)
      x(ixG^S,1)=s%x(ixG^S,1)
      drs(ixG^S)=s%dx(ixG^S,1)
      {^NOONED
      dx2(ixG^S)=s%dx(ixG^S,2)}
      {^IFTHREED
      dx3(ixG^S)=s%dx(ixG^S,3)}

      s%surfaceC(ixG^S,1)=dabs(x(ixG^S,1)+half*drs(ixG^S)){^DE&*dx^DE(ixG^S) }
      {^NOONED
      if (z_==2) s%surfaceC(ixG^S,2)=x(ixG^S,1)*drs(ixG^S){^IFTHREED*dx3(ixG^S)}
      if (phi_==2) s%surfaceC(ixG^S,2)=drs(ixG^S){^IFTHREED*dx3(ixG^S)}
      }
      {^IFTHREED
      if (z_==3) s%surfaceC(ixG^S,3)=x(ixG^S,1)*drs(ixG^S)*dx2(ixG^S)
      if (phi_==3) s%surfaceC(ixG^S,3)=drs(ixG^S)*dx2(ixG^S)
      }
      {^IFONED
      s%surfaceC(0,1)=dabs(x(1,1)-half*drs(1))
      }
      {^IFTWOD
      s%surfaceC(0,ixGmin2:ixGmax2,1)=dabs(x(1,ixGmin2:ixGmax2,1)-half*drs(1,&
         ixGmin2:ixGmax2))*dx2(1,ixGmin2:ixGmax2)
      }
      {^IFTHREED
      s%surfaceC(0,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1)=dabs(x(1,ixGmin2:ixGmax2,&
         ixGmin3:ixGmax3,1)-half*drs(1,ixGmin2:ixGmax2,ixGmin3:ixGmax3))*dx2(1,&
         ixGmin2:ixGmax2,ixGmin3:ixGmax3)*dx3(1,ixGmin2:ixGmax2,&
         ixGmin3:ixGmax3)
      }
      {s%surfaceC(0^DE%ixG^S,^DE)=s%surfaceC(1^DE%ixG^S,^DE);\}

      s%surface(ixG^S,1)=dabs(x(ixG^S,1)){^DE&*dx^DE(ixG^S) }
      {^NOONED
      if (z_==2) s%surface(ixG^S,2)=x(ixG^S,1)*drs(ixG^S){^IFTHREED*dx3(ixG^S)}
      if (phi_==2) s%surface(ixG^S,2)=drs(ixG^S){^IFTHREED*dx3(ixG^S)}}
      {^IFTHREED
      if (z_==3) s%surface(ixG^S,3)=x(ixG^S,1)*drs(ixG^S)*dx2(ixG^S)
      if (phi_==3) s%surface(ixG^S,3)=drs(ixG^S)*dx2(ixG^S)}

    case default
      call mpistop("Sorry, coordinate unknown")
    end select

  end subroutine get_surface_area

end module mod_geometry
