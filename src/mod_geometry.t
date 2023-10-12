!> Module with geometry-related routines (e.g., divergence, curl)
module mod_geometry
  implicit none
  public

  !> Indices for cylindrical and spherical coordinates, negative value when not used:
  integer :: r_     = -1
  integer :: theta_ = -1
  integer :: phi_   = -1
  integer :: z_     = -1

  integer :: coordinate=-1
  integer, parameter :: Cartesian          = 0
  integer, parameter :: Cartesian_stretched= 1
  integer, parameter :: cylindrical        = 2
  integer, parameter :: spherical          = 3

  integer :: type_curl=0
  integer, parameter :: central     = 1
  integer, parameter :: Gaussbased  = 2
  integer, parameter :: Stokesbased = 3

  ! --------- Metric variables --------- !
  !> Index of the lapse function
  integer, protected :: alp_ = -1
  !> Index of the conformal factor
  integer, protected :: psi_ = -1
  !> Indices of the shift vector 
  integer, protected, allocatable :: beta(:)
  !> Indices of the vector potential
  integer, protected, allocatable :: vecX(:)

  procedure(chris_indices), pointer :: i_chris => null()

  abstract interface
     integer function chris_indices(i,j,k)
       integer, intent(in) :: i,j,k
     end function chris_indices
  end interface

contains

  !> Set the coordinate system to be used
  subroutine set_coordinate_system(geom)
    use mod_global_parameters
    use mod_variables
    character(len=*), intent(in) :: geom !< Name of the coordinate system
    integer :: idir

    ! Store the geometry name
    geometry_name = geom

    select case (geom)
    case ("Cartesian","Cartesian_1D","Cartesian_2D","Cartesian_3D")
      ndir = ndim
      coordinate=Cartesian
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
      phi_ = 3
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
      z_   = 3
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
      ! Note that spherical 1.5D is not supported
      ndir   = ndim
      r_     = 1
      theta_ = 2
      phi_   = 3
      z_     = -1
      coordinate=spherical
    case ("spherical_2.5D")
      if (ndim /= 2) &
           call mpistop("Geometry spherical_2.5D requires ndim == 2")
      ndir    = 3
      r_      = 1
      theta_  = 2
      phi_    = 3
      z_      = -1
      coordinate=spherical
    case default
      call mpistop("Unknown geometry specified")
    end select

    call set_metric()
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

  !> calculate barycenter
  subroutine get_xbar(s,ixG^L)
    use mod_global_parameters

    type(mesh_t) :: s
    integer, intent(in) :: ixG^L

    double precision :: x(ixG^S,ndim), drs(ixG^S), dx2(ixG^S), dx3(ixG^S)
    integer          :: h1x^L{^NOONED, h2x^L}

    select case (coordinate)
    case (Cartesian_stretched)
      call mpistop("Cartesian_stretched is not supported yet")
    case (Cartesian)
      s%xbar=s%x
    case (cylindrical)
      !if ( z_ == 3 ) iphi = 2 !fixme: maybe put this in chris
      x(ixG^S,1)=s%x(ixG^S,1)
      drs(ixG^S)=s%dx(ixG^S,1)
      s%xbar(ixG^S,1)= s%x(ixG^S,1)
      !s%xbar(ixG^S,1)= (x(ixG^S,1)**2 + drs(ixG^S)**2/12.0d0) / x(ixG^S,1)
      {^NOONED
      s%xbar(ixG^S,2)= s%x(ixG^S,2)}
      {^IFTHREED
      s%xbar(ixG^S,3)= s%x(ixG^S,3)}
    case (spherical)
      x(ixG^S,1)=s%x(ixG^S,1)
      drs(ixG^S)=s%dx(ixG^S,1)
      s%xbar(ixG^S,1)= s%x(ixG^S,1)
      !s%xbar(ixG^S,1)= (x(ixG^S,1)**2 - 0.25d0 * drs(ixG^S)**2) * x(ixG^S,1) &
      !              / (x(ixG^S,1)**2 +drs(ixG^S)**2/12.0d0)
      !s%xbar(ixG^S,1)= x(ixG^S,1) - 4.0d0 * x(ixG^S,1)*drs(ixG^S)**2 &
      !              / (12.0d0*x(ixG^S,1)**2 +drs(ixG^S)**2)
      {^NOONED
      s%xbar(ixG^S,2)= s%x(ixG^S,2)}
      {^IFTHREED
      s%xbar(ixG^S,3)= s%x(ixG^S,3)}
      !call mpistop("Sorry, spherical not yet")
    case default
      call mpistop("Sorry, coordinate unknown")
    end select
  end subroutine get_xbar

  !> calculate area of surfaces of cells
  subroutine get_surface_area(s,ixG^L)
    use mod_global_parameters

    type(mesh_t) :: s
    integer, intent(in) :: ixG^L

    double precision :: x(ixG^S,ndim), drs(ixG^S), dx2(ixG^S), dx3(ixG^S)

    select case (coordinate)
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

      s%surfaceC(ixG^S,1)=(x(ixG^S,1)+0.5d0*drs(ixG^S))**2 {^NOONED &
           *2.0d0*dsin(x(ixG^S,2))*dsin(0.5d0*dx2(ixG^S))}{^IFTHREED*dx3(ixG^S)}

      {^NOONED
      s%surfaceC(ixG^S,2)=(x(ixG^S,1)**2 +drs(ixG^S)**2/12.0d0)*drs(ixG^S)&
           *dsin(x(ixG^S,2)+0.5d0*dx2(ixG^S))}{^IFTHREED*dx3(ixG^S)}

      {^IFTHREED
      s%surfaceC(ixG^S,3)=x(ixG^S,1)*drs(ixG^S)*dx2(ixG^S)
      }

      {^IFONED
      s%surfaceC(0,1)=dabs(x(1,1)-0.5d0*drs(1))**2
      }
      {^IFTWOD
      s%surfaceC(0,ixGmin2:ixGmax2,1)=(x(1,ixGmin2:ixGmax2,1)-0.5d0*drs(1,ixGmin2:ixGmax2))**2 &
         *2.0d0*dsin(x(1,ixGmin2:ixGmax2,2))*dsin(0.5d0*dx2(1,ixGmin2:ixGmax2))

      s%surfaceC(ixGmin1:ixGmax1,0,2)=(x(ixGmin1:ixGmax1,1,1)**2 +drs(ixGmin1:ixGmax1,1)**2/12.0d0) *drs(ixGmin1:ixGmax1,1)&
                       *dsin(x(ixGmin1:ixGmax1,1,2)-0.5d0*dx2(ixGmin1:ixGmax1,1))
      }
      {^IFTHREED
      s%surfaceC(0,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1)=(x(1,ixGmin2:ixGmax2,&
         ixGmin3:ixGmax3,1)-0.5d0*drs(1,ixGmin2:ixGmax2,&
         ixGmin3:ixGmax3))**2*2.0d0*dsin(x(1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
         2))*dsin(0.5d0*dx2(1,ixGmin2:ixGmax2,ixGmin3:ixGmax3))*dx3(1,&
         ixGmin2:ixGmax2,ixGmin3:ixGmax3)
      s%surfaceC(ixGmin1:ixGmax1,0,ixGmin3:ixGmax3,2)=x(ixGmin1:ixGmax1,1,&
         ixGmin3:ixGmax3,1)*drs(ixGmin1:ixGmax1,1,&
         ixGmin3:ixGmax3)*dsin(x(ixGmin1:ixGmax1,1,ixGmin3:ixGmax3,&
         2)-0.5d0*dx2(ixGmin1:ixGmax1,1,ixGmin3:ixGmax3))*dx3(ixGmin1:ixGmax1,1,&
         ixGmin3:ixGmax3)
      s%surfaceC(ixGmin1:ixGmax1,ixGmin2:ixGmax2,0,3)=&
         s%surfaceC(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1,3)
      }

      s%surface(ixG^S,1)=x(ixG^S,1)**2 {^NOONED &
           *2.0d0*dsin(x(ixG^S,2))*dsin(0.5d0*dx2(ixG^S))}{^IFTHREED*dx3(ixG^S)}
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

      s%surfaceC(ixG^S,1)=dabs(x(ixG^S,1)+0.5d0*drs(ixG^S)) {^DE&*dx^DE(ixG^S) }
      {^NOONED
      if (z_==2) s%surfaceC(ixG^S,2)=x(ixG^S,1)*drs(ixG^S){^IFTHREED*dx3(ixG^S)}
      if (phi_==2) s%surfaceC(ixG^S,2)=drs(ixG^S){^IFTHREED*dx3(ixG^S)}
      }
      {^IFTHREED
      if (z_==3) s%surfaceC(ixG^S,3)=x(ixG^S,1)*drs(ixG^S)*dx2(ixG^S)
      if (phi_==3) s%surfaceC(ixG^S,3)=drs(ixG^S)*dx2(ixG^S)
      }

      {^IFONED
      s%surfaceC(0,1)=dabs(x(1,1)-0.5d0*drs(1))
      }
      {^IFTWOD
      s%surfaceC(0,ixGmin2:ixGmax2,1)=dabs(x(1,ixGmin2:ixGmax2,1)-0.5d0*drs(1,&
         ixGmin2:ixGmax2))*dx2(1,ixGmin2:ixGmax2)
      }
      {^IFTHREED
      s%surfaceC(0,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1)=dabs(x(1,ixGmin2:ixGmax2,&
         ixGmin3:ixGmax3,1)-0.5d0*drs(1,ixGmin2:ixGmax2,ixGmin3:ixGmax3))*dx2(1,&
         ixGmin2:ixGmax2,ixGmin3:ixGmax3)*dx3(1,ixGmin2:ixGmax2,&
         ixGmin3:ixGmax3)
      }
      {s%surfaceC(0^DE%ixG^S,^DE)=s%surfaceC(1^DE%ixG^S,^DE);\}

      s%surface(ixG^S,1)=dabs(x(ixG^S,1)){^DE&*dx^DE(ixG^S) }
      {^NOONED
      s%surface(ixG^S,2)=x(ixG^S,1)*drs(ixG^S){^IFTHREED*dx3(ixG^S)}
      }
      {^IFTHREED
      s%surface(ixG^S,3)=x(ixG^S,1)*drs(ixG^S)*dx2(ixG^S)
      }

    case default
      call mpistop("Sorry, coordinate unknown")
    end select

  end subroutine get_surface_area

  !> calculate christoffel symbols of cells after calculating the surface area and volume
  subroutine get_christoffel(s,ixG^L)
    use mod_global_parameters

    type(mesh_t) :: s
    integer, intent(in) :: ixG^L

    double precision :: x(ixG^S,ndim), drs(ixG^S), dx2(ixG^S), dx3(ixG^S)
    integer          :: h1x^L{^NOONED, h2x^L}

    s%christoffel=0.0d0
    select case (coordinate)
    case (Cartesian,Cartesian_stretched)
       return ! nothing to do here
    case (spherical)
      h1x^L=ixG^L-kr(1,^D); {^NOONED h2x^L=ixG^L-kr(2,^D);}
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

      ! for (2,1,2), (2,2,1), (3,3,1), (3,1,3)
      s%christoffel(ixG^S,1)=0.5d0 * (s%surfaceC(ixG^S,1) - s%surfaceC(h1x^S,1)) / s%dvolume(ixG^S)

      ! for (1,2,2)
      s%christoffel(ixG^S,2)=-0.25d0 * ((x(ixG^S,1)+0.5d0*drs(ixG^S))**2*s%surfaceC(ixG^S,1) &
                                          - (x(ixG^S,1)-0.5d0*drs(ixG^S))**2*s%surfaceC(h1x^S,1)) / s%dvolume(ixG^S)

      ! for (1,3,3)
      s%christoffel(ixG^S,3)=-0.25d0 / s%dvolume(ixG^S) &
                           * ((x(ixG^S,1)+0.5d0*drs(ixG^S))**4 - (x(ixG^S,1)-0.5d0*drs(ixG^S))**4){^NOONED &
           *( (dcos(x(ixG^S,2)+0.5d0*dx2(ixG^S))**3/3.0d0 - dcos(x(ixG^S,2)+0.5d0*dx2(ixG^S))) &
             -(dcos(x(ixG^S,2)-0.5d0*dx2(ixG^S))**3/3.0d0 - dcos(x(ixG^S,2)-0.5d0*dx2(ixG^S))) )}{^IFTHREED*dx3(ixG^S)}

      {^NOONED
      ! for (2,3,3)
      s%christoffel(ixG^S,4)= - ( &
            (dsin(x(ixG^S,2)+0.5d0*dx2(ixG^S))**2)*s%surfaceC(ixG^S,2) &
          - (dsin(x(ixG^S,2)-0.5d0*dx2(ixG^S))**2)*s%surfaceC(h2x^S,2) &
              ) / 3.0d0 / s%dvolume(ixG^S) 

      ! for (3,2,3), (3,3,2)
      s%christoffel(ixG^S,5)=dcos(x(ixG^S,2))/dsin(x(ixG^S,2))
      }

      i_chris => i_chris_sph
    case (cylindrical)
      x(ixG^S,1)=s%x(ixG^S,1)
      drs(ixG^S)=s%dx(ixG^S,1)
      {^NOONED
      dx2(ixG^S)=s%dx(ixG^S,2)}
      {^IFTHREED
      dx3(ixG^S)=s%dx(ixG^S,3)}

      ! for (1,3,3)
      s%christoffel(ixG^S,1)= - (x(ixG^S,r_)**2 +drs(ixG^S)**2/12.0d0) / x(ixG^S,r_)

      ! for (3,1,3), (3,3,1)
      s%christoffel(ixG^S,2)= 1.0d0 / x(ixG^S,r_)

      i_chris => i_chris_cyl
    case default
      call mpistop("Sorry, coordinate unknown")
    end select

  end subroutine get_christoffel

  integer function i_chris_sph(i,j,k)
     integer, intent(in) :: i,j,k
     if ( i==theta_ .and. j==r_ .and. k==theta_) then
        i_chris_sph = 1
     else if ( i==theta_ .and. j==theta_ .and. k==r_) then
        i_chris_sph = 1
     else if ( i==phi_ .and. j==phi_ .and. k==r_) then
        i_chris_sph = 1
     else if ( i==phi_ .and. j==r_ .and. k==phi_) then
        i_chris_sph = 1
     else if ( i==r_ .and. j==theta_ .and. k==theta_) then
        i_chris_sph = 2
     else if ( i==r_ .and. j==phi_ .and. k==phi_) then
        i_chris_sph = 3
     else if ( i==theta_ .and. j==phi_ .and. k==phi_) then
        i_chris_sph = 4
     else if ( i==phi_ .and. j==theta_ .and. k==phi_) then
        i_chris_sph = 5
     else if ( i==phi_ .and. j==phi_ .and. k==theta_) then
        i_chris_sph = 5
     else
        i_chris_sph = 0
     end if
  end function i_chris_sph

  integer function i_chris_cyl(i,j,k)
     integer, intent(in) :: i,j,k
     if ( i==r_ .and. j==phi_ .and. k==phi_) then
        i_chris_cyl = 1
     else if ( i==phi_ .and. j==r_ .and. k==phi_) then
        i_chris_cyl = 2
     else if ( i==phi_ .and. j==phi_ .and. k==r_) then
        i_chris_cyl = 2
     else
        i_chris_cyl = 0
     end if
  end function i_chris_cyl

  !> get the 3-metric in flat space sqrt_gamma_hat 
  subroutine get_sqrt_gamma_hat(x, ixI^L, ixO^L, sqrt_gamma)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(out) :: sqrt_gamma(ixI^S)
    double precision, intent(in)  :: x(ixI^S, 1:^ND)

    select case (coordinate)
    case (Cartesian)
      sqrt_gamma(ixO^S) = 1.0d0
    case (cylindrical)
      sqrt_gamma(ixO^S) = dabs( x(ixO^S, 1) )
    case (spherical)
      sqrt_gamma(ixO^S) = dabs( x(ixO^S, 1)**2 {^NOONED * dsin(x(ixO^S, 2)) } )
    case default
      call mpistop("Sorry, coordinate unknown")
    end select
  end subroutine get_sqrt_gamma_hat

  !> get the 3-metric in flat space inv_sqrt_gamma_hat 
  subroutine get_inv_sqrt_gamma_hat(x, ixI^L, ixO^L, inv_sqrt_gamma)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(out) :: inv_sqrt_gamma(ixI^S)
    double precision, intent(in)  :: x(ixI^S, 1:^ND)

    select case (coordinate)
    case (Cartesian)
      inv_sqrt_gamma(ixO^S) = 1.0d0
    case (cylindrical)
      where (x(ixO^S, r_) /= 0.0d0)
         inv_sqrt_gamma(ixO^S) = 1.0d0 / dabs( x(ixO^S, 1) )
      else where
         inv_sqrt_gamma(ixO^S) = 0.0d0 
      end where
    case (spherical)
      where (x(ixO^S, r_) {^NOONED * dsin(x(ixO^S, 2)) } /= 0.0d0)
         inv_sqrt_gamma(ixO^S) = 1.0d0 / dabs( x(ixO^S, 1)**2 {^NOONED * dsin(x(ixO^S, 2)) })
      else where
         inv_sqrt_gamma(ixO^S) = 0.0d0 
      end where
    case default
      call mpistop("Sorry, coordinate unknown")
    end select
  end subroutine get_inv_sqrt_gamma_hat

  !> get the 3-metric in flat space gamma_ij_hat 
  subroutine get_gamma_ij_hat(x, ixI^L, ixO^L, gamma)
    ! Since the flat metric is always diagonal, 
    ! we store only one index
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(out)   :: gamma(ixI^S, 1:3)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)
    select case (coordinate)
    case (Cartesian)
      gamma(ixO^S,1) = 1.0d0
      gamma(ixO^S,2) = 1.0d0
      gamma(ixO^S,3) = 1.0d0
    case (cylindrical)
      gamma(ixO^S,r_)   = 1.0d0
      gamma(ixO^S,z_)   = 1.0d0
      gamma(ixO^S,phi_) = x(ixO^S, 1)**2
    case (spherical)
      gamma(ixO^S,r_)     = 1.0d0
      gamma(ixO^S,theta_) = x(ixO^S, 1)**2
      gamma(ixO^S,phi_)   = gamma(ixO^S,2) {^NOONED * dsin(x(ixO^S, 2))**2 }
    case default
      call mpistop("Sorry, coordinate unknown")
    end select
  end subroutine get_gamma_ij_hat

  !> get the 3-metric in flat space gammaij_hat 
  subroutine get_gammaij_hat(x, ixI^L, ixO^L, gamma)
    ! Since the flat metric is always diagonal, 
    ! we store only one index
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(out)   :: gamma(ixI^S, 1:3)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)
    select case (coordinate)
    case (Cartesian)
      gamma(ixO^S,1) = 1.0d0
      gamma(ixO^S,2) = 1.0d0
      gamma(ixO^S,3) = 1.0d0
    case (cylindrical)
      gamma(ixO^S,r_) = 1.0d0
      gamma(ixO^S,z_) = 1.0d0
      where (x(ixO^S, r_) /= 0.0d0)
         gamma(ixO^S,phi_) = 1.0d0 / x(ixO^S, 1)**2
      else where
         gamma(ixO^S,phi_) = 0.0d0 
      end where
    case (spherical)
      gamma(ixO^S,r_)     = 1.0d0
      where (x(ixO^S, r_) /= 0.0d0)
         gamma(ixO^S,theta_) = 1.0d0 / x(ixO^S, 1)**2
      else where
         gamma(ixO^S,theta_) = 0.0d0
      end where
      {^IFONED 
      gamma(ixO^S,phi_) = gamma(ixO^S,2)
      }
      {^NOONED 
      where (x(ixO^S, theta_) /= 0.0d0 )
         gamma(ixO^S,phi_) = gamma(ixO^S,2) / dsin(x(ixO^S, 2))**2
      else where
         gamma(ixO^S,phi_) = 0.0d0
      end where
      }
    case default
      call mpistop("Sorry, coordinate unknown")
    end select
  end subroutine get_gammaij_hat

  !> get the partial dervitives 3-metric in flat space gamma_ij_hat 
  subroutine get_dgamma_ij_hat(x, ixI^L, ixO^L, dgamma)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(out)   :: dgamma(ixI^S, 1:3, 1:3, 1:3)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)

    dgamma = 0.0d0

    select case (coordinate)
    case (cylindrical)
      ! for partial_r
      if ( z_ == 2 ) then
         dgamma(ixO^S,3,3,r_) = 2.0d0 * x(ixO^S, 1)
      else
         dgamma(ixO^S,2,2,r_) = 2.0d0 * x(ixO^S, 1)
      end if
    case (spherical)
      ! note: r_ = 1,theta_=2, phi_ = 3.
      ! for partial_r
      dgamma(ixO^S,2,2,r_) = 2.0d0 * x(ixO^S, 1)
      dgamma(ixO^S,3,3,r_) = dgamma(ixO^S,2,2,r_) {^NOONED * dsin(x(ixO^S, 2))**2 }

      {^NOONED 
      ! for partial_theta
      dgamma(ixO^S,3,3,theta_) = 2.0d0 * x(ixO^S, 1)**2 &
               * dsin(x(ixO^S, 2)) * dcos(x(ixO^S, 2)) 
      }
    case default
      call mpistop("Sorry, coordinate unknown")
    end select
  end subroutine get_dgamma_ij_hat

  !> transform vectors between natural and orthonomal basis 
  subroutine get_natural2orthonormal(x, ixI^L, ixO^L, N2R)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S, 1:^ND)
    double precision, intent(out)   :: N2R(ixI^S, 1:3)

    N2R(ixO^S,1:3) = 1.0d0

    select case (coordinate)
    case (Cartesian)
      ! nothing to do here
    case (cylindrical)
      if ( z_ == 2 ) then
         N2R(ixO^S,3) = x(ixO^S, 1)
      else
         N2R(ixO^S,2) = x(ixO^S, 1)
      end if
    case (spherical)
      ! note: r_ = 1,theta_=2, phi_ = 3.
      N2R(ixO^S,2) = x(ixO^S, 1)
      N2R(ixO^S,3) = x(ixO^S, 1) {^NOONED * dsin(x(ixO^S, 2)) }
    case default
      call mpistop("Sorry, coordinate unknown")
    end select
  end subroutine get_natural2orthonormal

  subroutine center_to_faces(q,ixI^L,ixC^L,idims,qC)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixC^L,idims
    double precision, intent(in)    :: q(ixI^S)
    double precision, intent(out)   :: qC(ixI^S)
    integer                         :: jxC^L
    jxC^L=ixC^L+kr(idims,^D);
    ! linear interpolation at cell interface
    qC(ixC^S)=0.5d0*(q(ixC^S)+q(jxC^S))
    !qC(ixC^S)=(q(ixC^S)+0.5d0*mesh_in%dx(ixC^S,idims)*&
    !           (q(jxC^S)-q(ixC^S))/(mesh_in%x(jxC^S,idims)-mesh_in%x(ixC^S,idims)))
  end subroutine center_to_faces

  !> Calculate dervitives of a scalar q within ixL in direction idir
  subroutine partial_d(q,ixI^L,ixO^L,idir,mesh_in,gradq)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: q(ixI^S)
    type(mesh_t), intent(in)        :: mesh_in
    double precision, intent(inout) :: gradq(ixI^S)
    integer                         :: jxO^L, hxO^L

    associate( x => mesh_in%x )

    hxO^L=ixO^L-kr(idir,^D);
    jxO^L=ixO^L+kr(idir,^D);
    select case(coordinate)
    case(Cartesian)
      gradq(ixO^S)=0.5d0*(q(jxO^S)-q(hxO^S))/dxlevel(idir)
    case(Cartesian_stretched)
      gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,idir)-x(hxO^S,idir))
    case(spherical)
      select case(idir)
      case(1)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,1)-x(hxO^S,1)))
        {^NOONED
      case(2)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,2)-x(hxO^S,2))
        }
        {^IFTHREED
      case(3)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,3)-x(hxO^S,3))
        }
      end select
    case(cylindrical)
      if(idir==phi_) then
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,phi_)-x(hxO^S,phi_))!*x(ixO^S,r_))
      else
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,idir)-x(hxO^S,idir))
      end if
    case default
      call mpistop('Unknown geometry')
    end select

    end associate
  end subroutine partial_d

  !> Calculate gradient of a scalar q within ixL in direction idir
  subroutine gradient(q,ixI^L,ixO^L,idir,mesh_in,gradq)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: q(ixI^S)
    type(mesh_t), intent(in)        :: mesh_in
    double precision, intent(inout) :: gradq(ixI^S)
    integer                         :: jxO^L, hxO^L

    associate( x => mesh_in%x )

    hxO^L=ixO^L-kr(idir,^D);
    jxO^L=ixO^L+kr(idir,^D);
    select case(coordinate)
    case(Cartesian)
      gradq(ixO^S)=half*(q(jxO^S)-q(hxO^S))/dxlevel(idir)
    case(Cartesian_stretched)
      gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,idir)-x(hxO^S,idir))
    case(spherical)
      select case(idir)
      case(1)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,1)-x(hxO^S,1)))
        {^NOONED
      case(2)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,2)-x(hxO^S,2))*x(ixO^S,1))
        }
        {^IFTHREED
      case(3)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,3)-x(hxO^S,3))*x(ixO^S,1)*dsin(x(ixO^S,2)))
        }
      end select
    case(cylindrical)
      if(idir==phi_) then
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,phi_)-x(hxO^S,phi_))*x(ixO^S,r_))
      else
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,idir)-x(hxO^S,idir))
      end if
    case default
      call mpistop('Unknown geometry')
    end select

    end associate
  end subroutine gradient

  !> Calculate gradient of a scalar q in direction idir at cell interfaces
  subroutine gradientx(q,x,ixI^L,ixO^L,idir,gradq,fourth_order)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: q(ixI^S), x(ixI^S,1:ndim)
    double precision, intent(inout) :: gradq(ixI^S)
    logical, intent(in)             :: fourth_order
    integer                         :: jxO^L, hxO^L, kxO^L

    if(fourth_order) then
      ! Fourth order, stencil width is two
      kxO^L=ixO^L^LADD2;
      if(ixImin^D>kxOmin^D.or.ixImax^D<kxOmax^D|.or.) &
           call mpistop("Error in gradientx: Non-conforming input limits")
      hxO^L=ixO^L-kr(idir,^D);
      jxO^L=ixO^L+kr(idir,^D);
      kxO^L=ixO^L+kr(idir,^D)*2;
    else
      hxO^L=ixO^L;
    end if
    jxO^L=ixO^L+kr(idir,^D);
    select case(coordinate)
    case(Cartesian)
      if(fourth_order) then
        gradq(ixO^S)=(27.d0*(q(jxO^S)-q(ixO^S))-q(kxO^S)+q(hxO^S))/24.d0/dxlevel(idir)
      else
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/dxlevel(idir)
      end if
    case(Cartesian_stretched)
      gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,idir)-x(hxO^S,idir))
    case(spherical)
      select case(idir)
      case(1)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,1)-x(hxO^S,1)))
        {^NOONED
      case(2)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,2)-x(hxO^S,2))*x(ixO^S,1))
        }
        {^IFTHREED
      case(3)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,3)-x(hxO^S,3))*x(ixO^S,1)*dsin(x(ixO^S,2)))
        }
      end select
    case(cylindrical)
      if(idir==phi_) then
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,phi_)-x(hxO^S,phi_))*x(ixO^S,r_))
      else
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,idir)-x(hxO^S,idir))
      end if
    case default
      call mpistop('Unknown geometry')
    end select

  end subroutine gradientx

  !> Calculate gradient of a scalar q within ixL in direction idir
  !> first use limiter to go from cell center to edge
  subroutine gradientS(q,ixI^L,ixO^L,idir,mesh_in,gradq)
    use mod_global_parameters
    use mod_limiter
    use mod_ppm

    integer, intent(in)                :: ixI^L, ixO^L, idir
    double precision, intent(in)       :: q(ixI^S)
    type(mesh_t), intent(in)           :: mesh_in
    double precision, intent(inout)    :: gradq(ixI^S)
    double precision ,dimension(ixI^S) :: qC,qL,qR,dqC,ldq,rdq

    double precision :: invdx
    integer          :: hxO^L,ixC^L,jxC^L,gxC^L,hxC^L

    associate( x => mesh_in%x )

    invdx=1.d0/dxlevel(idir)
    hxO^L=ixO^L-kr(idir,^D);
    ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
    jxC^L=ixC^L+kr(idir,^D);
    gxCmin^D=ixCmin^D-kr(idir,^D);gxCmax^D=jxCmax^D;
    hxC^L=gxC^L+kr(idir,^D);

    ! set the gradient limiter here
    qR(gxC^S) = q(hxC^S)
    qL(gxC^S) = q(gxC^S)
    if (type_gradient_limiter(mesh_in%level)/=limiter_ppm) then
      dqC(gxC^S)= qR(gxC^S)-qL(gxC^S)
      call dwlimiter2(dqC,ixI^L,gxC^L,idir,type_gradient_limiter(mesh_in%level),ldw=ldq,rdw=rdq)
      qL(ixC^S) = qL(ixC^S) + 0.5d0*ldq(ixC^S)
      qR(ixC^S) = qR(ixC^S) - 0.5d0*rdq(jxC^S)
    else
      call PPMlimiter(ixI^L,ixO^L,idir,q,qL,qR,PPM_extrema)
    endif

    select case(coordinate)
    case(Cartesian)
      gradq(ixO^S)=(qR(ixO^S)-qL(hxO^S))*invdx
    case(Cartesian_stretched)
      gradq(ixO^S)=(qR(ixO^S)-qL(hxO^S))/mesh_in%dx(ixO^S,idir)
    case(spherical)
      gradq(ixO^S)=(qR(ixO^S)-qL(hxO^S))/mesh_in%dx(ixO^S,idir)
      select case(idir)
      case(2)
        gradq(ixO^S)=gradq(ixO^S)/x(ixO^S,1)
        {^IFTHREED
      case(3)
        gradq(ixO^S)=gradq(ixO^S)/(x(ixO^S,1)*dsin(x(ixO^S,2)))
        }
      end select
    case(cylindrical)
      gradq(ixO^S)=(qR(ixO^S)-qL(hxO^S))/mesh_in%dx(ixO^S,idir)
      if(idir==phi_) gradq(ixO^S)=gradq(ixO^S)/x(ixO^S,1)
    end select

    end associate
  end subroutine gradientS

  !> Calculate divergence of a vector qvec within ixL
  subroutine divvector(qvec,ixI^L,ixO^L,mesh_in,divq,fourthorder,sixthorder)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qvec(ixI^S,1:ndir)
    type(mesh_t), intent(in)        :: mesh_in
    double precision, intent(inout) :: divq(ixI^S)
    logical, intent(in), optional   :: fourthorder !< Default: false
    logical, intent(in), optional   :: sixthorder !< Default: false
    logical                         :: use_4th_order
    logical                         :: use_6th_order
    double precision                :: qC(ixI^S), invdx(1:ndim)
    integer                         :: jxO^L, hxO^L, ixC^L, jxC^L
    integer                         :: idims, ix^L, gxO^L, kxO^L
    integer                         :: lxO^L, fxO^L

    use_4th_order = .false.
    use_6th_order = .false.
    if (present(fourthorder)) use_4th_order = fourthorder
    if (present(sixthorder))  use_6th_order = sixthorder
    if(use_4th_order .and. use_6th_order) &
      call mpistop("divvector: using 4th and 6th order at the same time")

    if(use_4th_order) then
      if (.not. slab_uniform) &
           call mpistop("divvector: 4th order only supported for slab geometry")
      ! Fourth order, stencil width is two
      ix^L=ixO^L^LADD2;
    else if(use_6th_order) then
      ! Sixth order, stencil width is three
      if (.not. slab_uniform) &
           call mpistop("divvector: 6th order only supported for slab geometry")
      ix^L=ixO^L^LADD3;
    else
      ! Second order, stencil width is one
      ix^L=ixO^L^LADD1;
    end if

    if (ixImin^D>ixmin^D.or.ixImax^D<ixmax^D|.or.) &
         call mpistop("Error in divvector: Non-conforming input limits")

    invdx=1.d0/dxlevel
    divq(ixO^S)=0.0d0

    if (slab_uniform) then
      do idims=1,ndim
        if(use_6th_order) then
          lxO^L=ixO^L+3*kr(idims,^D);
          kxO^L=ixO^L+2*kr(idims,^D);
          jxO^L=ixO^L+kr(idims,^D);
          hxO^L=ixO^L-kr(idims,^D);
          gxO^L=ixO^L-2*kr(idims,^D);
          fxO^L=ixO^L-3*kr(idims,^D);
          divq(ixO^S)=divq(ixO^S)+&
                    (-qvec(fxO^S,idims)+9.d0*qvec(gxO^S,idims)-45.d0*qvec(hxO^S,idims)&
                     +qvec(lxO^S,idims)-9.d0*qvec(kxO^S,idims)+45.d0*qvec(jxO^S,idims))&
                     /(60.d0*dxlevel(idims))
        else if(use_4th_order) then
          ! Use fourth order scheme
          kxO^L=ixO^L+2*kr(idims,^D);
          jxO^L=ixO^L+kr(idims,^D);
          hxO^L=ixO^L-kr(idims,^D);
          gxO^L=ixO^L-2*kr(idims,^D);
          divq(ixO^S)=divq(ixO^S)+&
                    (-qvec(kxO^S,idims)+8.d0*qvec(jxO^S,idims)-8.d0*&
                      qvec(hxO^S,idims)+qvec(gxO^S,idims))/(12.d0*dxlevel(idims))
        else
          ! Use second order scheme
          jxO^L=ixO^L+kr(idims,^D);
          hxO^L=ixO^L-kr(idims,^D);
          divq(ixO^S)=divq(ixO^S)+half*(qvec(jxO^S,idims) &
               - qvec(hxO^S,idims))*invdx(idims)
        end if
      end do
    else
      do idims=1,ndim
        hxO^L=ixO^L-kr(idims,^D);
        ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
        jxC^L=ixC^L+kr(idims,^D);
        if(stretched_dim(idims) .and. stretch_uncentered) then
          ! linear interpolation at cell interface along stretched dimension
          qC(ixC^S)=mesh_in%surfaceC(ixC^S,idims)*(qvec(ixC^S,idims)+0.5d0*mesh_in%dx(ixC^S,idims)*&
               (qvec(jxC^S,idims)-qvec(ixC^S,idims))/(mesh_in%x(jxC^S,idims)-mesh_in%x(ixC^S,idims)))
        else
          qC(ixC^S)=mesh_in%surfaceC(ixC^S,idims)*half*(qvec(ixC^S,idims)+qvec(jxC^S,idims))
        end if
        divq(ixO^S)=divq(ixO^S)+qC(ixO^S)-qC(hxO^S)
      end do
      divq(ixO^S)=divq(ixO^S)/mesh_in%dvolume(ixO^S)
    end if
  end subroutine divvector

  !> Calculate divergence of a vector qvec within ixL
  !> using limited extrapolation to cell edges
  subroutine divvectorS(qvec,ixI^L,ixO^L,mesh_in,divq)
    use mod_global_parameters
    use mod_limiter
    use mod_ppm

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: qvec(ixI^S,1:ndir)
    type(mesh_t), intent(in)           :: mesh_in
    double precision, intent(inout)    :: divq(ixI^S)
    double precision, dimension(ixI^S) :: qL,qR,dqC,ldq,rdq

    double precision :: invdx(1:ndim)
    integer          :: hxO^L,ixC^L,jxC^L,idims,ix^L,gxC^L,hxC^L

    ix^L=ixO^L^LADD2;

    if (ixImin^D>ixmin^D.or.ixImax^D<ixmax^D|.or.) &
         call mpistop("Error in divvectorS: Non-conforming input limits")

    invdx=1.d0/dxlevel
    divq(ixO^S)=zero
    do idims=1,ndim
      hxO^L=ixO^L-kr(idims,^D);
      ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
      jxC^L=ixC^L+kr(idims,^D);
      gxCmin^D=ixCmin^D-kr(idims,^D);gxCmax^D=jxCmax^D;
      hxC^L=gxC^L+kr(idims,^D);
      qR(gxC^S) = qvec(hxC^S,idims)
      qL(gxC^S) = qvec(gxC^S,idims)
      if(type_gradient_limiter(mesh_in%level)/=limiter_ppm) then
        dqC(gxC^S)= qR(gxC^S)-qL(gxC^S)
        call dwlimiter2(dqC,ixI^L,gxC^L,idims,type_gradient_limiter(mesh_in%level),ldw=ldq,rdw=rdq)
        qL(ixC^S) = qL(ixC^S) + 0.5d0*ldq(ixC^S)
        qR(ixC^S) = qR(ixC^S) - 0.5d0*rdq(jxC^S)
      else
        dqC(ixI^S)=qvec(ixI^S,idims)
        call PPMlimiter(ixI^L,ixO^L,idims,dqC,qL,qR,PPM_extrema)
      endif

      if (slab_uniform) then
        divq(ixO^S)=divq(ixO^S)+0.5d0*(qR(ixO^S)-qL(hxO^S))*invdx(idims)
      else
        qR(ixC^S)=mesh_in%surfaceC(ixC^S,idims)*qR(ixC^S)
        qL(ixC^S)=mesh_in%surfaceC(ixC^S,idims)*qL(ixC^S)
        divq(ixO^S)=divq(ixO^S)+qR(ixO^S)-qL(hxO^S)
      end if
    end do
    if(.not.slab_uniform) divq(ixO^S)=divq(ixO^S)/mesh_in%dvolume(ixO^S)

  end subroutine divvectorS

  !> Calculate curl of a vector qvec within ixL
  !> Options to 
  !>        employ standard second order CD evaluations
  !>        use Gauss theorem for non-Cartesian grids
  !>        use Stokes theorem for non-Cartesian grids
  subroutine curlvector(qvec,ixI^L,ixO^L,mesh_in,curlvec,idirmin,idirmin0,ndir0,fourthorder)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    integer, intent(in)             :: ndir0, idirmin0
    integer, intent(inout)          :: idirmin
    type(mesh_t), intent(in)        :: mesh_in
    double precision, intent(in)    :: qvec(ixI^S,1:ndir0)
    double precision, intent(inout) :: curlvec(ixI^S,idirmin0:3)
    logical, intent(in), optional   :: fourthorder !< Default: false

    integer          :: ixA^L,ixC^L,jxC^L,idir,jdir,kdir,hxO^L,jxO^L,kxO^L,gxO^L
    double precision :: invdx(1:ndim)
    double precision :: tmp(ixI^S),tmp2(ixI^S),xC(ixI^S),surface(ixI^S)
    logical          :: use_4th_order

    ! Calculate curl within ixL: CurlV_i=eps_ijk*d_j V_k
    ! Curl can have components (idirmin:3)
    ! Determine exact value of idirmin while doing the loop.

    use_4th_order = .false.
    if (present(fourthorder)) use_4th_order = fourthorder

    if (use_4th_order) then
      if (.not. slab_uniform) &
           call mpistop("curlvector: 4th order only supported for slab geometry")
      ! Fourth order, stencil width is two
      ixA^L=ixO^L^LADD2;
    else
      ! Second order, stencil width is one
      ixA^L=ixO^L^LADD1;
    end if

    if (ixImin^D>ixAmin^D.or.ixImax^D<ixAmax^D|.or.) &
         call mpistop("Error in curlvector: Non-conforming input limits")

    idirmin=4
    curlvec(ixO^S,idirmin0:3)=zero

    ! all non-Cartesian cases
    select case(coordinate)
      case(Cartesian) ! Cartesian grids
        invdx=1.d0/dxlevel
        do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
          if(lvc(idir,jdir,kdir)/=0)then
            if (.not. use_4th_order) then
              ! Use second order scheme
              jxO^L=ixO^L+kr(jdir,^D);
              hxO^L=ixO^L-kr(jdir,^D);
              tmp(ixO^S)=half*(qvec(jxO^S,kdir) &
                   - qvec(hxO^S,kdir))*invdx(jdir)
            else
              ! Use fourth order scheme
              kxO^L=ixO^L+2*kr(jdir,^D);
              jxO^L=ixO^L+kr(jdir,^D);
              hxO^L=ixO^L-kr(jdir,^D);
              gxO^L=ixO^L-2*kr(jdir,^D);
              tmp(ixO^S)=(-qvec(kxO^S,kdir) + 8.0d0 * qvec(jxO^S,kdir) - 8.0d0 * &
                   qvec(hxO^S,kdir) + qvec(gxO^S,kdir))/(12.0d0 * dxlevel(jdir))
            end if
            if(lvc(idir,jdir,kdir)==1)then
              curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp(ixO^S)
            else
              curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp(ixO^S)
            endif
            if(idir<idirmin)idirmin=idir
          endif
        enddo; enddo; enddo;
      case(Cartesian_stretched) ! stretched Cartesian grids
        do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
          if(lvc(idir,jdir,kdir)/=0)then
            select case(type_curl)
              case(central)
                tmp(ixA^S)=qvec(ixA^S,kdir)
                hxO^L=ixO^L-kr(jdir,^D);
                jxO^L=ixO^L+kr(jdir,^D);
                ! second order centered differencing
                tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(mesh_in%x(jxO^S,jdir)-mesh_in%x(hxO^S,jdir))
              case(Gaussbased)
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                tmp(ixC^S)=mesh_in%surfaceC(ixC^S,jdir)*(qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir)))
                tmp2(ixO^S)=(tmp(ixO^S)-tmp(hxO^S))/mesh_in%dvolume(ixO^S)
              case(Stokesbased)
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                if(kdir<=ndim)then
                  tmp(ixC^S)=mesh_in%ds(ixO^S,kdir)*(qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir)))
                else
                  tmp(ixC^S)=(qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir)))
                endif
                if(idir<=ndim)then
                  tmp2(ixO^S)=(tmp(ixO^S)-tmp(hxO^S))/mesh_in%surface(ixO^S,idir)
                else ! essentially for 2.5D case, idir=3 and jdir,kdir<=2
                  tmp2(ixO^S)=(tmp(ixO^S)-tmp(hxO^S))/(mesh_in%ds(ixO^S,jdir)*mesh_in%ds(ixO^S,kdir))
                endif
              case default
                call mpistop('no such curl evaluator')
            end select
            if(lvc(idir,jdir,kdir)==1)then
              curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
            else
              curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
            endif
            if(idir<idirmin)idirmin=idir
          endif
        enddo; enddo; enddo;
      case(spherical) ! possibly stretched spherical grids
        select case(type_curl)
          case(central) ! ok for any dimensionality
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                tmp(ixA^S)=qvec(ixA^S,kdir)
                hxO^L=ixO^L-kr(jdir,^D);
                jxO^L=ixO^L+kr(jdir,^D);
                select case(jdir)
                case(1)
                tmp(ixA^S)=tmp(ixA^S)*mesh_in%x(ixA^S,1)
                tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((mesh_in%x(jxO^S,1)-mesh_in%x(hxO^S,1))*mesh_in%x(ixO^S,1))
                {^NOONED    case(2)
                if(idir==1) tmp(ixA^S)=tmp(ixA^S)*dsin(mesh_in%x(ixA^S,2))
                tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((mesh_in%x(jxO^S,2)-mesh_in%x(hxO^S,2))*mesh_in%x(ixO^S,1))
                if(idir==1) tmp2(ixO^S)=tmp2(ixO^S)/dsin(mesh_in%x(ixO^S,2))
                }
                {^IFTHREED  case(3)
                tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((mesh_in%x(jxO^S,3)-mesh_in%x(hxO^S,3))*mesh_in%x(ixO^S,1)*dsin(mesh_in%x(ixO^S,2)))
                }
                end select
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
                else
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
                endif
                if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case(Gaussbased)
            do idir=idirmin0,3;
            do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                tmp(ixC^S)=mesh_in%surfaceC(ixC^S,jdir)*(qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir)))
                tmp2(ixO^S)=(tmp(ixO^S)-tmp(hxO^S))/mesh_in%dvolume(ixO^S)
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
                else
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo;
            ! geometric terms
            if(idir==2.and.phi_>0) curlvec(ixO^S,2)=curlvec(ixO^S,2)+qvec(ixO^S,phi_)/mesh_in%x(ixO^S,r_)
            {^NOONED
            if(idir==phi_) curlvec(ixO^S,phi_)=curlvec(ixO^S,phi_)-qvec(ixO^S,2)/mesh_in%x(ixO^S,r_) &
                 +qvec(ixO^S,r_)*dcos(mesh_in%x(ixO^S,2))/(mesh_in%x(ixO^S,r_)*dsin(mesh_in%x(ixO^S,2)))
            }
            enddo;
          case(Stokesbased)
            !if(ndim<3) call mpistop("Stokesbased for 3D spherical only")
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
             if(lvc(idir,jdir,kdir)/=0)then
              select case(idir)
              case(1)
                if(jdir<kdir) then
                  ! idir=1,jdir=2,kdir=3
                  !! integral along 3rd dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(3) at cell interface along 2nd dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  ! 2nd coordinate at cell interface along 2nd dimension
                  xC(ixC^S)=mesh_in%x(ixC^S,jdir)+0.5d0*mesh_in%dx(ixC^S,jdir)
                  curlvec(ixO^S,idir)=(dsin(xC(ixO^S))*tmp(ixO^S)-&
                       dsin(xC(hxO^S))*tmp(hxO^S))*mesh_in%dx(ixO^S,kdir)
                  !! integral along 2nd dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(2) at cell interface along 3rd dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*mesh_in%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(mesh_in%x(jxC^S,kdir)-mesh_in%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*mesh_in%dx(ixO^S,jdir))&
                       /mesh_in%surface(ixO^S,idir)*mesh_in%x(ixO^S,idir)
                end if
              case(2)
                if(jdir<kdir) then
                  ! idir=2,jdir=1,kdir=3
                  !! integral along 1st dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(1) at cell interface along 3rd dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*mesh_in%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(mesh_in%x(jxC^S,kdir)-mesh_in%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*mesh_in%dx(ixO^S,1)
                  !! integral along 3rd dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(3) at cell interface along 1st dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  ! 1st coordinate at cell interface along 1st dimension
                  xC(ixC^S)=mesh_in%x(ixC^S,jdir)+0.5d0*mesh_in%dx(ixC^S,jdir)
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(xC(hxO^S)*tmp(hxO^S)-xC(ixO^S)*tmp(ixO^S))*&
                       dsin(mesh_in%x(ixO^S,idir))*mesh_in%dx(ixO^S,kdir))/mesh_in%surface(ixO^S,idir)
                end if
              case(3)
                if(jdir<kdir) then
                  ! idir=3,jdir=1,kdir=2
                  !! integral along 1st dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(1) at cell interface along 2nd dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*mesh_in%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(mesh_in%x(jxC^S,kdir)-mesh_in%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(hxO^S)-tmp(ixO^S))*mesh_in%dx(ixO^S,jdir)
                  !! integral along 2nd dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(2) at cell interface along 1st dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  ! 1st coordinate at cell interface along 1st dimension
                  xC(ixC^S)=mesh_in%x(ixC^S,jdir)+0.5d0*mesh_in%dx(ixC^S,jdir)
                  if(ndim==3) then
                    surface(ixO^S)=mesh_in%surface(ixO^S,idir)
                  else
                    surface(ixO^S)=mesh_in%x(ixO^S,jdir)*mesh_in%dx(ixO^S,kdir)*mesh_in%dx(ixO^S,jdir)
                  end if
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(xC(ixO^S)*tmp(ixO^S)-xC(hxO^S)*tmp(hxO^S))*mesh_in%dx(ixO^S,kdir))&
                       /surface(ixO^S)
                end if
              end select
              if(idir<idirmin)idirmin=idir
             endif
            enddo; enddo; enddo;
          case default
            call mpistop('no such curl evaluator')
        end select
      case(cylindrical) ! possibly stretched cylindrical grids
        select case(type_curl)
          case(central)  ! works for any dimensionality, polar/cylindrical
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                tmp(ixA^S)=qvec(ixA^S,kdir)
                hxO^L=ixO^L-kr(jdir,^D);
                jxO^L=ixO^L+kr(jdir,^D);
                if(z_==3.or.z_==-1) then
                  ! Case Polar_2D, Polar_2.5D or Polar_3D, i.e. R,phi,Z
                  select case(jdir)
                  case(1)
                  if(idir==z_) tmp(ixA^S)=tmp(ixA^S)*mesh_in%x(ixA^S,1) ! R V_phi
                  ! computes d(R V_phi)/dR or d V_Z/dR
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(mesh_in%x(jxO^S,1)-mesh_in%x(hxO^S,1))
                  if(idir==z_) tmp2(ixO^S)=tmp2(ixO^S)/mesh_in%x(ixO^S,1) ! (1/R)*d(R V_phi)/dR
                  {^NOONED      case(2)
                  ! handles (1/R)d V_Z/dphi or (1/R)d V_R/dphi
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((mesh_in%x(jxO^S,2)-mesh_in%x(hxO^S,2))*mesh_in%x(ixO^S,1))
                  }
                  {^IFTHREED    case(3)
                  ! handles d V_phi/dZ or d V_R/dZ
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(mesh_in%x(jxO^S,3)-mesh_in%x(hxO^S,3))
                  }
                  end select
                end if
                if(phi_==3.or.phi_==-1) then
                  ! Case Cylindrical_2D, Cylindrical_2.5D or Cylindrical_3D, i.e. R,Z,phi
                  select case(jdir)
                  case(1)
                  if(idir==z_) tmp(ixA^S)=tmp(ixA^S)*mesh_in%x(ixA^S,1) ! R V_phi
                  ! computes d(R V_phi)/dR or d V_Z/dR
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(mesh_in%x(jxO^S,1)-mesh_in%x(hxO^S,1))
                  if(idir==z_) tmp2(ixO^S)=tmp2(ixO^S)/mesh_in%x(ixO^S,1) ! (1/R)*d(R V_phi)/dR
                  {^NOONED      case(2)
                  ! handles d V_phi/dZ or d V_R/dZ
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(mesh_in%x(jxO^S,2)-mesh_in%x(hxO^S,2))
                  }
                  {^IFTHREED    case(3)
                  ! handles (1/R)d V_Z/dphi or (1/R)d V_R/dphi
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((mesh_in%x(jxO^S,3)-mesh_in%x(hxO^S,3))*mesh_in%x(ixO^S,1))
                  }
                  end select
                end if
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
                else
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case(Gaussbased) ! works for any dimensionality, polar/cylindrical
            if(ndim<2) call mpistop("Gaussbased for 2D, 2.5D or 3D polar or cylindrical only")
            do idir=idirmin0,3;
            do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                tmp(ixC^S)=mesh_in%surfaceC(ixC^S,jdir)*(qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir)))
                tmp2(ixO^S)=(tmp(ixO^S)-tmp(hxO^S))/mesh_in%dvolume(ixO^S)
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
                else
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo;
            ! geometric term from d e_R/d phi= e_phi for unit vectors e_R, e_phi
            !       but minus sign appears due to R,Z,phi ordering (?)
            ! note that in cylindrical 2D (RZ), phi_ is -1
            ! note that in polar 2D     (Rphi), z_ is -1
            if((idir==phi_.or.(phi_==-1.and.idir==3)).and.z_>0) then
              ! cylindrical
              if(  z_==2) curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-qvec(ixO^S,z_)/mesh_in%x(ixO^S,r_)
              ! polar
              if(phi_==2) curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+qvec(ixO^S,z_)/mesh_in%x(ixO^S,r_)
            endif
            enddo;
          case(Stokesbased)
            if(ndim<3) call mpistop("Stokesbased for 3D cylindrical only")
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
               if(idir==r_) then
                if(jdir==phi_) then
                  ! idir=r,jdir=phi,kdir=z
                  !! integral along z dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(z) at cell interface along phi dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*mesh_in%dx(ixO^S,kdir)
                  !! integral along phi dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(phi) at cell interface along z dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*mesh_in%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(mesh_in%x(jxC^S,kdir)-mesh_in%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*mesh_in%x(ixO^S,idir)*mesh_in%dx(ixO^S,jdir))&
                       /mesh_in%surface(ixO^S,idir)
                end if
               else if(idir==phi_) then
                if(jdir<kdir) then
                  ! idir=phi,jdir=r,kdir=z
                  !! integral along r dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(r) at cell interface along z dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*mesh_in%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(mesh_in%x(jxC^S,kdir)-mesh_in%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*mesh_in%dx(ixO^S,jdir)
                  !! integral along z dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(z) at cell interface along r dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*mesh_in%dx(ixO^S,kdir))&
                       /mesh_in%surface(ixO^S,idir)
                end if
               else ! idir==z_
                if(jdir<kdir) then
                  ! idir=z,jdir=r,kdir=phi
                  !! integral along r dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(r) at cell interface along phi dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*mesh_in%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(mesh_in%x(jxC^S,kdir)-mesh_in%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(hxO^S)-tmp(ixO^S))*mesh_in%dx(ixO^S,jdir)
                  !! integral along phi dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(phi) at cell interface along r dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  ! r coordinate at cell interface along r dimension
                  xC(ixC^S)=mesh_in%x(ixC^S,jdir)+0.5d0*mesh_in%dx(ixC^S,jdir)
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(xC(ixO^S)*tmp(ixO^S)-xC(hxO^S)*tmp(hxO^S))*mesh_in%dx(ixO^S,kdir))&
                       /mesh_in%surface(ixO^S,idir)
                end if
               end if
               if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case default
            call mpistop('no such curl evaluator')
        end select
        case default
          call mpistop('not possible to calculate curl')
    end select

  end subroutine curlvector

  !> Calculate idim transverse components of curl of a vector qvec within ixL
  !> Options to
  !>        employ standard second order CD evaluations
  !>        use Gauss theorem for non-Cartesian grids
  !>        use Stokes theorem for non-Cartesian grids
  subroutine curlvector_trans(qvec,qvecc,ixI^L,ixO^L,mesh_in,curlvec,idim,idirmin,idirmin0,ndir0)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    integer, intent(in)             :: idim, ndir0, idirmin0
    integer, intent(inout)          :: idirmin
    type(mesh_t), intent(in)        :: mesh_in
    double precision, intent(in)    :: qvec(ixI^S,1:ndir0),qvecc(ixI^S,1:ndir0)
    double precision, intent(inout) :: curlvec(ixI^S,idirmin0:3)

    integer          :: ixA^L,ixC^L,jxC^L,idir,jdir,kdir,hxO^L,jxO^L
    double precision :: invdx(1:ndim)
    double precision :: tmp(ixI^S),tmp2(ixI^S),xC(ixI^S),surface(ixI^S)

    ! Calculate curl within ixL: CurlV_i=eps_ijk*d_j V_k
    ! Curl can have components (idirmin:3)
    ! Determine exact value of idirmin while doing the loop.

    idirmin=4
    curlvec(ixO^S,idirmin0:3)=zero
    ! Second order, stencil width is one
    ixA^L=ixO^L^LADD1;

    ! all non-Cartesian cases
    select case(coordinate)
      case(Cartesian) ! Cartesian grids
        invdx=1.d0/dxlevel
        do idir=idirmin0,3
          do jdir=1,ndim; do kdir=1,ndir0
            if(lvc(idir,jdir,kdir)/=0)then
              if(jdir/=idim) then
                tmp(ixI^S)=qvec(ixI^S,kdir)
                hxO^L=ixO^L-kr(jdir,^D);
                jxO^L=ixO^L+kr(jdir,^D);
              else
                ! use two cell-center values for gradient at interface of the two cells
                ! because left and right neighbor interface values is unavailable at block boundary faces
                tmp(ixI^S)=qvecc(ixI^S,kdir)
                hxO^L=ixO^L;
                jxO^L=ixO^L+kr(jdir,^D);
              end if
              ! second order centered differencing
              tmp(ixO^S)=half*(tmp(jxO^S)-tmp(hxO^S))*invdx(jdir)
              !> \todo allow for 4th order CD evaluation here as well
              if(lvc(idir,jdir,kdir)==1)then
                curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp(ixO^S)
              else
                curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp(ixO^S)
              endif
              if(idir<idirmin)idirmin=idir
            endif
          enddo; enddo;
        end do
      case(Cartesian_stretched) ! stretched Cartesian grids
        do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
          if(lvc(idir,jdir,kdir)/=0)then
            select case(type_curl)
              case(central)
                tmp(ixI^S)=qvec(ixI^S,kdir)
                hxO^L=ixO^L-kr(jdir,^D);
                jxO^L=ixO^L+kr(jdir,^D);
                ! second order centered differencing
                tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(mesh_in%x(jxO^S,jdir)-mesh_in%x(hxO^S,jdir))
              case(Gaussbased)
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                tmp(ixC^S)=mesh_in%surfaceC(ixC^S,jdir)*(qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir)))
                tmp2(ixO^S)=(tmp(ixO^S)-tmp(hxO^S))/mesh_in%dvolume(ixO^S)
              case(Stokesbased)
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                if(kdir<=ndim)then
                  tmp(ixC^S)=mesh_in%ds(ixO^S,kdir)*(qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir)))
                else
                  tmp(ixC^S)=(qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir)))
                endif
                if(idir<=ndim)then
                  tmp2(ixO^S)=(tmp(ixO^S)-tmp(hxO^S))/mesh_in%surface(ixO^S,idir)
                else ! essentially for 2.5D case, idir=3 and jdir,kdir<=2
                  tmp2(ixO^S)=(tmp(ixO^S)-tmp(hxO^S))/(mesh_in%ds(ixO^S,jdir)*mesh_in%ds(ixO^S,kdir))
                endif
              case default
                call mpistop('no such curl evaluator')
            end select
            if(lvc(idir,jdir,kdir)==1)then
              curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
            else
              curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
            endif
            if(idir<idirmin)idirmin=idir
          endif
        enddo; enddo; enddo;
      case(spherical) ! possibly stretched spherical grids
        select case(type_curl)
          case(central) ! ok for any dimensionality
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                tmp(ixI^S)=qvec(ixI^S,kdir)
                hxO^L=ixO^L-kr(jdir,^D);
                jxO^L=ixO^L+kr(jdir,^D);
                select case(jdir)
                case(1)
                tmp(ixA^S)=tmp(ixA^S)*mesh_in%x(ixA^S,1)
                tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((mesh_in%x(jxO^S,1)-mesh_in%x(hxO^S,1))*mesh_in%x(ixO^S,1))
                {^NOONED    case(2)
                if(idir==1) tmp(ixA^S)=tmp(ixA^S)*dsin(mesh_in%x(ixA^S,2))
                tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((mesh_in%x(jxO^S,2)-mesh_in%x(hxO^S,2))*mesh_in%x(ixO^S,1))
                if(idir==1) tmp2(ixO^S)=tmp2(ixO^S)/dsin(mesh_in%x(ixO^S,2))
                }
                {^IFTHREED  case(3)
                tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((mesh_in%x(jxO^S,3)-mesh_in%x(hxO^S,3))*mesh_in%x(ixO^S,1)*dsin(mesh_in%x(ixO^S,2)))
                }
                end select
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
                else
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
                endif
                if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case(Gaussbased)
            do idir=idirmin0,3;
            do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                tmp(ixC^S)=mesh_in%surfaceC(ixC^S,jdir)*(qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir)))
                tmp2(ixO^S)=(tmp(ixO^S)-tmp(hxO^S))/mesh_in%dvolume(ixO^S)
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
                else
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo;
            ! geometric terms
            if(idir==2.and.phi_>0) curlvec(ixO^S,2)=curlvec(ixO^S,2)+qvec(ixO^S,phi_)/mesh_in%x(ixO^S,r_)
            {^NOONED
            if(idir==phi_) curlvec(ixO^S,phi_)=curlvec(ixO^S,phi_)-qvec(ixO^S,2)/mesh_in%x(ixO^S,r_) &
                 +qvec(ixO^S,r_)*dcos(mesh_in%x(ixO^S,2))/(mesh_in%x(ixO^S,r_)*dsin(mesh_in%x(ixO^S,2)))
            }
            enddo;
          case(Stokesbased)
            !if(ndim<3) call mpistop("Stokesbased for 3D spherical only")
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
             if(lvc(idir,jdir,kdir)/=0)then
              select case(idir)
              case(1)
                if(jdir<kdir) then
                  ! idir=1,jdir=2,kdir=3
                  !! integral along 3rd dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(3) at cell interface along 2nd dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  ! 2nd coordinate at cell interface along 2nd dimension
                  xC(ixC^S)=mesh_in%x(ixC^S,jdir)+0.5d0*mesh_in%dx(ixC^S,jdir)
                  curlvec(ixO^S,idir)=(dsin(xC(ixO^S))*tmp(ixO^S)-&
                       dsin(xC(hxO^S))*tmp(hxO^S))*mesh_in%dx(ixO^S,kdir)
                  !! integral along 2nd dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(2) at cell interface along 3rd dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*mesh_in%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(mesh_in%x(jxC^S,kdir)-mesh_in%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*mesh_in%dx(ixO^S,jdir))&
                       /mesh_in%surface(ixO^S,idir)*mesh_in%x(ixO^S,idir)
                end if
              case(2)
                if(jdir<kdir) then
                  ! idir=2,jdir=1,kdir=3
                  !! integral along 1st dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(1) at cell interface along 3rd dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*mesh_in%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(mesh_in%x(jxC^S,kdir)-mesh_in%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*mesh_in%dx(ixO^S,1)
                  !! integral along 3rd dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(3) at cell interface along 1st dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  ! 1st coordinate at cell interface along 1st dimension
                  xC(ixC^S)=mesh_in%x(ixC^S,jdir)+0.5d0*mesh_in%dx(ixC^S,jdir)
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(xC(hxO^S)*tmp(hxO^S)-xC(ixO^S)*tmp(ixO^S))*&
                       dsin(mesh_in%x(ixO^S,idir))*mesh_in%dx(ixO^S,kdir))/mesh_in%surface(ixO^S,idir)
                end if
              case(3)
                if(jdir<kdir) then
                  ! idir=3,jdir=1,kdir=2
                  !! integral along 1st dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(1) at cell interface along 2nd dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*mesh_in%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(mesh_in%x(jxC^S,kdir)-mesh_in%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(hxO^S)-tmp(ixO^S))*mesh_in%dx(ixO^S,jdir)
                  !! integral along 2nd dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(2) at cell interface along 1st dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  ! 1st coordinate at cell interface along 1st dimension
                  xC(ixC^S)=mesh_in%x(ixC^S,jdir)+0.5d0*mesh_in%dx(ixC^S,jdir)
                  if(ndim==3) then
                    surface(ixO^S)=mesh_in%surface(ixO^S,idir)
                  else
                    surface(ixO^S)=mesh_in%x(ixO^S,jdir)*mesh_in%dx(ixO^S,kdir)*mesh_in%dx(ixO^S,jdir)
                  end if
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(xC(ixO^S)*tmp(ixO^S)-xC(hxO^S)*tmp(hxO^S))*mesh_in%dx(ixO^S,kdir))&
                       /surface(ixO^S)
                end if
              end select
              if(idir<idirmin)idirmin=idir
             endif
            enddo; enddo; enddo;
          case default
            call mpistop('no such curl evaluator')
        end select
      case(cylindrical) ! possibly stretched cylindrical grids
        select case(type_curl)
          case(central)  ! works for any dimensionality, polar/cylindrical
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                tmp(ixI^S)=qvec(ixI^S,kdir)
                hxO^L=ixO^L-kr(jdir,^D);
                jxO^L=ixO^L+kr(jdir,^D);
                if(z_==3.or.z_==-1) then
                  ! Case Polar_2D, Polar_2.5D or Polar_3D, i.e. R,phi,Z
                  select case(jdir)
                  case(1)
                  if(idir==z_) tmp(ixA^S)=tmp(ixA^S)*mesh_in%x(ixA^S,1) ! R V_phi
                  ! computes d(R V_phi)/dR or d V_Z/dR
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(mesh_in%x(jxO^S,1)-mesh_in%x(hxO^S,1))
                  if(idir==z_) tmp2(ixO^S)=tmp2(ixO^S)/mesh_in%x(ixO^S,1) ! (1/R)*d(R V_phi)/dR
                  {^NOONED      case(2)
                  ! handles (1/R)d V_Z/dphi or (1/R)d V_R/dphi
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((mesh_in%x(jxO^S,2)-mesh_in%x(hxO^S,2))*mesh_in%x(ixO^S,1))
                  }
                  {^IFTHREED    case(3)
                  ! handles d V_phi/dZ or d V_R/dZ
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(mesh_in%x(jxO^S,3)-mesh_in%x(hxO^S,3))
                  }
                  end select
                end if
                if(phi_==3.or.phi_==-1) then
                  ! Case Cylindrical_2D, Cylindrical_2.5D or Cylindrical_3D, i.e. R,Z,phi
                  select case(jdir)
                  case(1)
                  if(idir==z_) tmp(ixA^S)=tmp(ixA^S)*mesh_in%x(ixA^S,1) ! R V_phi
                  ! computes d(R V_phi)/dR or d V_Z/dR
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(mesh_in%x(jxO^S,1)-mesh_in%x(hxO^S,1))
                  if(idir==z_) tmp2(ixO^S)=tmp2(ixO^S)/mesh_in%x(ixO^S,1) ! (1/R)*d(R V_phi)/dR
                  {^NOONED      case(2)
                  ! handles d V_phi/dZ or d V_R/dZ
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(mesh_in%x(jxO^S,2)-mesh_in%x(hxO^S,2))
                  }
                  {^IFTHREED    case(3)
                  ! handles (1/R)d V_Z/dphi or (1/R)d V_R/dphi
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((mesh_in%x(jxO^S,3)-mesh_in%x(hxO^S,3))*mesh_in%x(ixO^S,1))
                  }
                  end select
                end if
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
                else
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case(Gaussbased) ! works for any dimensionality, polar/cylindrical
            if(ndim<2) call mpistop("Gaussbased for 2D, 2.5D or 3D polar or cylindrical only")
            do idir=idirmin0,3;
            do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                tmp(ixC^S)=mesh_in%surfaceC(ixC^S,jdir)*(qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir)))
                tmp2(ixO^S)=(tmp(ixO^S)-tmp(hxO^S))/mesh_in%dvolume(ixO^S)
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
                else
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo;
            ! geometric term from d e_R/d phi= e_phi for unit vectors e_R, e_phi
            !       but minus sign appears due to R,Z,phi ordering (?)
            ! note that in cylindrical 2D (RZ), phi_ is -1
            ! note that in polar 2D     (Rphi), z_ is -1
            if((idir==phi_.or.(phi_==-1.and.idir==3)).and.z_>0) then
              ! cylindrical
              if(  z_==2) curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-qvec(ixO^S,z_)/mesh_in%x(ixO^S,r_)
              ! polar
              if(phi_==2) curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+qvec(ixO^S,z_)/mesh_in%x(ixO^S,r_)
            endif
            enddo;
          case(Stokesbased)
            if(ndim<3) call mpistop("Stokesbased for 3D cylindrical only")
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
               if(idir==r_) then
                if(jdir==phi_) then
                  ! idir=r,jdir=phi,kdir=z
                  !! integral along z dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(z) at cell interface along phi dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*mesh_in%dx(ixO^S,kdir)
                  !! integral along phi dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(phi) at cell interface along z dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*mesh_in%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(mesh_in%x(jxC^S,kdir)-mesh_in%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*mesh_in%x(ixO^S,idir)*mesh_in%dx(ixO^S,jdir))&
                       /mesh_in%surface(ixO^S,idir)
                end if
               else if(idir==phi_) then
                if(jdir<kdir) then
                  ! idir=phi,jdir=r,kdir=z
                  !! integral along r dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(r) at cell interface along z dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*mesh_in%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(mesh_in%x(jxC^S,kdir)-mesh_in%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*mesh_in%dx(ixO^S,jdir)
                  !! integral along z dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(z) at cell interface along r dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*mesh_in%dx(ixO^S,kdir))&
                       /mesh_in%surface(ixO^S,idir)
                end if
               else ! idir==z_
                if(jdir<kdir) then
                  ! idir=z,jdir=r,kdir=phi
                  !! integral along r dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(r) at cell interface along phi dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*mesh_in%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(mesh_in%x(jxC^S,kdir)-mesh_in%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(hxO^S)-tmp(ixO^S))*mesh_in%dx(ixO^S,jdir)
                  !! integral along phi dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(phi) at cell interface along r dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*mesh_in%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(mesh_in%x(jxC^S,jdir)-mesh_in%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  ! r coordinate at cell interface along r dimension
                  xC(ixC^S)=mesh_in%x(ixC^S,jdir)+0.5d0*mesh_in%dx(ixC^S,jdir)
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(xC(ixO^S)*tmp(ixO^S)-xC(hxO^S)*tmp(hxO^S))*mesh_in%dx(ixO^S,kdir))&
                       /mesh_in%surface(ixO^S,idir)
                end if
               end if
               if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case default
            call mpistop('no such curl evaluator')
        end select
        case default
          call mpistop('not possible to calculate curl')
    end select

  end subroutine curlvector_trans

  double precision function get_radii_pt(x)
    use mod_global_parameters
    double precision, intent(in) :: x(1:^ND)
    select case (coordinate)
    case (spherical)
       get_radii_pt = x(r_)
    case (cylindrical)
       get_radii_pt = dsqrt( x(r_)**2 {^NOONED + x(z_)**2 } )
    case default
       get_radii_pt = dsqrt( sum(x(:)**2) )
    end select
  end function get_radii_pt

  ! ------------------------------
  ! Metric related subroutine
  ! fixme: move this to somewhere else
  ! ------------------------------

  !> Set the coordinate system to be used
  subroutine set_metric()
    use mod_global_parameters
    use mod_variables
    integer :: idir

    ! create metric variables
    ! Set index of metric variables
    alp_ = var_set_metricvar("alp")
    psi_ = var_set_metricvar("psi")
    allocate(beta(ndir))
    do idir=1,ndir
       beta(idir) = var_set_metricvar("beta", idir)
    end do
    allocate(vecX(ndir))
    do idir=1,ndir
       vecX(idir) = var_set_metricvar("X", idir)
    end do
  end subroutine set_metric

  subroutine get_3metric(ixI^L, ixO^L, metric_in, mesh_in, &
        sqrt_gma_over_gma_hat, sqrt_gamma, &
        gamma_ij, gammaij )
    use mod_global_parameters
    integer, intent(in)                     :: ixI^L, ixO^L
    type(metric_t), intent(in)              :: metric_in
    type(mesh_t), intent(in)                :: mesh_in
    double precision, intent(out), optional :: sqrt_gma_over_gma_hat(ixI^S)
    double precision, intent(out), optional :: sqrt_gamma(ixI^S)
    double precision, intent(out), optional :: gamma_ij(ixI^S, 1:3, 1:3)
    double precision, intent(out), optional :: gammaij(ixI^S, 1:3, 1:3)

    double precision :: flat_metric(ixI^S, 1:3)
    integer          :: idir, jdir

    associate(m => metric_in%vars)

    ! so far we assume CFC, extend this later
    if (present(sqrt_gma_over_gma_hat)) &
       sqrt_gma_over_gma_hat(ixO^S) = m(ixO^S, psi_)**6

    ! fixme: optimise me!
    if (present(sqrt_gamma)) then
       call get_sqrt_gamma_hat(mesh_in%x, ixI^L, ixO^L, sqrt_gamma)
       sqrt_gamma(ixO^S) = m(ixO^S, psi_)**6 * sqrt_gamma(ixO^S)
    end if

    if (present(gamma_ij)) then
       call get_gamma_ij_hat(mesh_in%x, ixI^L, ixO^L, flat_metric)
       do idir = 1, 3; do jdir = 1, 3
          if (idir==jdir) then
             gamma_ij(ixO^S, idir, idir) = flat_metric(ixO^S, idir) &
                      * m(ixO^S, psi_)**4
          else
             gamma_ij(ixO^S, idir, jdir) = 0.0d0
          end if
       end do; end do
    end if

    if (present(gammaij)) then
       call get_gammaij_hat(mesh_in%x, ixI^L, ixO^L, flat_metric)
       do idir = 1, 3; do jdir = 1, 3
          if (idir==jdir) then
             gammaij(ixO^S, idir, idir) = flat_metric(ixO^S, idir) &
                      / m(ixO^S, psi_)**4
          else
             gammaij(ixO^S, idir, jdir) = 0.0d0
          end if
       end do; end do
    end if

    end associate
  end subroutine get_3metric

  !> get the gammaii
  subroutine get_gammaii(ixI^L, ixO^L, idim, metric_in, mesh_in, gammaii)
    ! Since the flat metric is always diagonal, 
    ! we store only one index
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    integer, intent(in)           :: idim
    type(metric_t), intent(in)    :: metric_in
    type(mesh_t), intent(in)      :: mesh_in
    double precision, intent(out) :: gammaii(ixI^S)

    double precision :: gammaij_hat(ixI^S, 1:3)

    ! so far we assume CFC, extend this later
    call get_gammaij_hat(mesh_in%x, ixI^L, ixO^L, gammaij_hat)
    gammaii(ixO^S) = gammaij_hat(ixO^S,idim) / metric_in%vars(ixO^S, psi_)**4

  end subroutine get_gammaii

  !> get the gamma_ii
  subroutine get_gamma_ii(ixI^L, ixO^L, idim, metric_in, mesh_in, gamma_ii)
    ! Since the flat metric is always diagonal, 
    ! we store only one index
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    integer, intent(in)           :: idim
    type(metric_t), intent(in)    :: metric_in
    type(mesh_t), intent(in)      :: mesh_in
    double precision, intent(out) :: gamma_ii(ixI^S)

    double precision :: gamma_ij_hat(ixI^S, 1:3)

    ! so far we assume CFC, extend this later
    call get_gamma_ij_hat(mesh_in%x, ixI^L, ixO^L, gamma_ij_hat)
    gamma_ii(ixO^S) = gamma_ij_hat(ixO^S,idim) * metric_in%vars(ixO^S, psi_)**4

  end subroutine get_gamma_ii

  subroutine contract_vec(ixI^L, ixO^L, n_dir, metric_in, mesh_in, ai, bj, a_dot_b)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    integer, intent(in)           :: n_dir
    type(metric_t), intent(in)    :: metric_in
    type(mesh_t), intent(in)      :: mesh_in
    double precision, intent(in)  :: ai(ixI^S,n_dir), bj(ixI^S,n_dir)
    double precision, intent(out) :: a_dot_b(ixI^S)

    double precision :: gamma_ij_hat(ixI^S, 1:3)
    integer          :: idir

    ! fixme: so far we assume CFC, extend this later
    call get_gamma_ij_hat(mesh_in%x, ixI^L, ixO^L, gamma_ij_hat)
    do idir = 1, n_dir
       gamma_ij_hat(ixO^S,idir) = metric_in%vars(ixO^S, psi_)**4 * gamma_ij_hat(ixO^S,idir)
    end do

    a_dot_b(ixO^S) = 0.0d0
    do idir = 1, n_dir
       a_dot_b(ixO^S) = a_dot_b(ixO^S) &
            + gamma_ij_hat(ixO^S,idir) * ai(ixO^S,idir) * bj(ixO^S,idir)
    end do
  end subroutine contract_vec

  subroutine contract_covec(ixI^L, ixO^L, n_dir, metric_in, mesh_in, a_i, b_j, a_dot_b)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    integer, intent(in)           :: n_dir
    type(metric_t), intent(in)    :: metric_in
    type(mesh_t), intent(in)      :: mesh_in
    double precision, intent(in)  :: a_i(ixI^S,n_dir), b_j(ixI^S,n_dir)
    double precision, intent(out) :: a_dot_b(ixI^S)

    double precision :: gammaij_hat(ixI^S, 1:3)
    integer          :: idir

    ! fixme: so far we assume CFC, extend this later
    call get_gammaij_hat(mesh_in%x, ixI^L, ixO^L, gammaij_hat)
    do idir = 1, n_dir
       gammaij_hat(ixO^S,idir) = gammaij_hat(ixO^S,idir) / metric_in%vars(ixO^S, psi_)**4 
    end do

    a_dot_b(ixO^S) = 0.0d0
    do idir = 1, n_dir
       a_dot_b(ixO^S) = a_dot_b(ixO^S) &
            + gammaij_hat(ixO^S,idir) * a_i(ixO^S,idir) * b_j(ixO^S,idir)
    end do
  end subroutine contract_covec

  subroutine contract_vec_pt(n_dir, gamma_ij, ai, bj, a_dot_b)
    use mod_global_parameters
    integer, intent(in)           :: n_dir
    double precision, intent(in)  :: gamma_ij(1:3,1:3)
    double precision, intent(in)  :: ai(1:n_dir), bj(1:n_dir)
    double precision, intent(out) :: a_dot_b

    integer          :: idir

    ! fixme: so far we assume CFC, extend this later
    ! where gamma_ij is diagonal
    a_dot_b = 0.0d0
    do idir = 1, n_dir
       a_dot_b = a_dot_b &
            + gamma_ij(idir,idir) * ai(idir) * bj(idir)
    end do
  end subroutine contract_vec_pt

  ! fixme: since the metric gammaij is given,
  ! this is actually the same as contract_vec_pt.
  ! maybe combine these two?
  subroutine contract_covec_pt(n_dir, gammaij, a_i, b_j, a_dot_b)
    use mod_global_parameters
    integer, intent(in)           :: n_dir
    double precision, intent(in)  :: gammaij(1:3,1:3)
    double precision, intent(in)  :: a_i(1:n_dir), b_j(1:n_dir)
    double precision, intent(out) :: a_dot_b

    integer          :: idir

    ! fixme: so far we assume CFC, extend this later
    ! where gamma_ij is diagonal
    a_dot_b = 0.0d0
    do idir = 1, n_dir
       a_dot_b = a_dot_b &
            + gammaij(idir,idir) * a_i(idir) * b_j(idir)
    end do
  end subroutine contract_covec_pt

  subroutine lower_indices(ixI^L, ixO^L, n_dir, metric_in, mesh_in, ai)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    integer, intent(in)             :: n_dir
    type(metric_t), intent(in)      :: metric_in
    type(mesh_t), intent(in)        :: mesh_in
    double precision, intent(inout) :: ai(ixI^S,1:n_dir)

    double precision :: gamma_ij(ixI^S, 1:3, 1:3)
    integer          :: idir

    ! fixme: so far we assume CFC, extend this later
    call get_3metric(ixI^L, ixO^L, metric_in, mesh_in, &
         gamma_ij=gamma_ij )
    do idir = 1, n_dir
       ai(ixO^S,idir) = ai(ixO^S,idir) * gamma_ij(ixO^S,idir,idir)
    end do
  end subroutine lower_indices

  subroutine lower_indices_pt(n_dir, gamma_ij, ai)
    use mod_global_parameters
    integer, intent(in)             :: n_dir
    double precision, intent(in)    :: gamma_ij(1:3,1:3)
    double precision, intent(inout) :: ai(1:n_dir)

    integer          :: idir

    ! fixme: so far we assume CFC, extend this later
    do idir = 1, n_dir
       ai(idir) = ai(idir) * gamma_ij(idir,idir)
    end do

  end subroutine lower_indices_pt

  subroutine raise_i(ixI^L, ixO^L, n_dir, metric_in, mesh_in, a_i, ai, idim)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    integer, intent(in)           :: n_dir
    type(metric_t), intent(in)    :: metric_in
    type(mesh_t), intent(in)      :: mesh_in
    double precision, intent(in)  :: a_i(ixI^S,n_dir)
    double precision, intent(out) :: ai(ixI^S)
    integer, intent(in)           :: idim

    ! fixme: so far we assume CFC, extend this later
    call get_gammaii(ixI^L, ixO^L, idim, metric_in, mesh_in, ai)
    ai(ixO^S) = ai(ixO^S) * a_i(ixO^S,idim)

  end subroutine raise_i

  subroutine raise_indices(ixI^L, ixO^L, n_dir, metric_in, mesh_in, a_i)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    integer, intent(in)             :: n_dir
    type(metric_t), intent(in)      :: metric_in
    type(mesh_t), intent(in)        :: mesh_in
    double precision, intent(inout) :: a_i(ixI^S,1:n_dir)

    double precision :: gammaij(ixI^S, 1:3, 1:3)
    integer          :: idir

    ! fixme: so far we assume CFC, extend this later
    call get_3metric(ixI^L, ixO^L, metric_in, mesh_in, &
         gammaij=gammaij )
    do idir = 1, n_dir
       a_i(ixO^S,idir) = a_i(ixO^S,idir) * gammaij(ixO^S,idir,idir)
    end do

  end subroutine raise_indices

  subroutine raise_indices_pt(n_dir, gammaij, a_i)
    use mod_global_parameters
    integer, intent(in)             :: n_dir
    double precision, intent(in)    :: gammaij(1:3,1:3)
    double precision, intent(inout) :: a_i(1:n_dir)

    integer          :: idir

    ! fixme: so far we assume CFC, extend this later
    do idir = 1, n_dir
       a_i(idir) = a_i(idir) * gammaij(idir,idir)
    end do

  end subroutine raise_indices_pt

  subroutine get_4metric(ixI^L, ixO^L, metric_in, mesh_in, gmunu, g_munu)
    use mod_global_parameters
    integer, intent(in)                     :: ixI^L, ixO^L
    type(metric_t), intent(in)              :: metric_in
    type(mesh_t), intent(in)                :: mesh_in
    double precision, intent(out), optional :: gmunu(ixI^S, 0:3, 0:3)
    double precision, intent(out), optional :: g_munu(ixI^S, 0:3, 0:3)
    integer                                 :: mu, nu
    double precision                        :: gamma(ixI^S, 1:3, 1:3)
    associate( m => metric_in%vars)
    if (present(gmunu)) then
       gmunu(ixO^S, 0:3,0:3) = 0.0d0
       gmunu(ixO^S, 0,0) = - 1.0d0 / m(ixO^S, alp_)**2
       do mu = 1, ndir
          gmunu(ixO^S, 0, mu) = - m(ixO^S, beta(mu)) * gmunu(ixO^S, 0, 0)
          gmunu(ixO^S, mu, 0) = gmunu(ixO^S, 0, mu)  
       end do
       ! get \gamma^{ij}
       call get_3metric(ixI^L, ixO^L, metric_in, mesh_in, gammaij=gamma)
       ! this is actually \gamma^{ij}
       gmunu(ixO^S, 1:3, 1:3) = gamma(ixO^S, 1:3, 1:3)
       ! g^{ij} = \gamma^{ij} - \beta^i \beta^j / \alpha^2
       do mu = 1, ndir; do nu = 1, ndir
          gmunu(ixO^S, mu, nu) = gmunu(ixO^S, mu, nu) &
               + m(ixO^S, beta(mu)) * m(ixO^S, beta(nu)) * gmunu(ixO^S, 0, 0)
       end do; end do
    end if

    if (present(g_munu)) then
       g_munu(ixO^S, 0:3,0:3) = 0.0d0
       g_munu(ixO^S, 0,0) = - m(ixO^S, alp_)**2
       ! fixme: use raise/lower indices
       ! get \gamma_{ij}
       call get_3metric(ixI^L, ixO^L, metric_in, mesh_in, gamma_ij=gamma)
       do mu = 1, ndir
          g_munu(ixO^S, 0, 0) = g_munu(ixO^S, 0, 0) + m(ixO^S, beta(mu))**2 * gamma(ixO^S, mu, mu)
       end do
       do mu = 1, ndir
          g_munu(ixO^S, 0, mu) = m(ixO^S, beta(mu)) * gamma(ixO^S, mu, mu)
          g_munu(ixO^S, mu, 0) = g_munu(ixO^S, 0, mu)  
       end do
       ! this is actually \gamma_{ij}
       g_munu(ixO^S, 1:3, 1:3) = gamma(ixO^S, 1:3, 1:3)
    end if
    end associate
  end subroutine get_4metric

  subroutine get_metric_dervitives(ixI^L, ixO^L, metric_in, mesh_in, &
        dalp_out, dpsi_out, dbeta_out, &
        Dgamma_out, D_beta_out, &
        K_ij )
    use mod_global_parameters
    integer, intent(in)                     :: ixI^L, ixO^L
    type(metric_t), intent(in)              :: metric_in
    type(mesh_t), intent(in)                :: mesh_in
    ! dervitaves of the metric variables
    double precision, intent(out), optional :: dalp_out(ixI^S, 1:3)
    double precision, intent(out), optional :: dpsi_out(ixI^S, 1:3)
    double precision, intent(out), optional :: dbeta_out(ixI^S, 1:3, 1:3)
    ! covariant dervitaves of beta D_i beta^k
    double precision, intent(out), optional :: D_beta_out(ixI^S, 1:3, 1:3)
    ! dervitaves of the 3-metric gamma_{jk},i
    double precision, intent(out), optional :: Dgamma_out(ixI^S, 1:3, 1:3, 1:3)
    ! extrinsic curvature K_ij
    double precision, intent(out), optional :: K_ij(ixI^S, 1:3, 1:3)

    {^NOONED
    double precision :: cot_theta(ixI^S)
    }
    integer          :: idir, jdir, kdir, i_ch
    double precision :: gamma_ij(ixI^S, 1:3, 1:3)
    double precision :: dalp(ixI^S, 1:3)
    double precision :: dpsi(ixI^S, 1:3)
    double precision :: dbeta(ixI^S, 1:3, 1:3)
    ! covariant dervitaves of beta D_i beta^k
    double precision :: D_beta(ixI^S, 1:3, 1:3)
    ! dervitaves of the 3-metric gamma_{jk},i
    double precision :: Dgamma(ixI^S, 1:3, 1:3, 1:3)

    associate( m => metric_in%vars, &
               x => mesh_in%x, chris_hat => mesh_in%christoffel )

    ! fixme: so far we assume CFC, extend this later

    call get_3metric(ixI^L, ixO^L, metric_in, mesh_in, &
           gamma_ij=gamma_ij )

    ! calculate derivitives of the metric variables
    dalp(ixI^S,1:3)      = 0.0d0
    dbeta(ixI^S,1:3,1:3) = 0.0d0
    dpsi(ixI^S,1:3)      = 0.0d0
    do idir = 1, ndim
       call partial_d( m(ixI^S,alp_) ,ixI^L,ixO^L,idir,mesh_in,dalp(ixI^S,idir) )
       call partial_d( m(ixI^S,psi_) ,ixI^L,ixO^L,idir,mesh_in,dpsi(ixI^S,idir) )
       dpsi(ixO^S,idir) = dpsi(ixO^S,idir) / m(ixO^S,psi_)
       do jdir = 1, ndir
          call partial_d( m(ixI^S,beta(jdir)) ,ixI^L,ixO^L,idir,mesh_in,dbeta(ixI^S,jdir,idir) )
       end do
    end do
 
    ! covariant derivative of beta: partial_i beta^k + Gamma^k_{ij} beta^j
    D_beta(ixO^S, 1:3, 1:3) = dbeta(ixO^S, 1:3, 1:3)
    if ( coordinate /= cartesian ) then
       do idir = 1, 3; do kdir = 1, 3
          do jdir = 1, ndir
             i_ch = i_chris(kdir,idir,jdir)
             if (i_ch==0) cycle
             D_beta(ixO^S, kdir, idir) = D_beta(ixO^S, kdir, idir) &
                        + chris_hat(ixO^S,i_ch) * m(ixO^S,beta(jdir))
          end do
       end do; end do
    end if

    ! dervitaves of metric D_i gamma_jk
    Dgamma(ixI^S,1:3,1:3,1:3) = 0.0d0
    do kdir = 1, ndim
       do idir = 1, 3
          Dgamma(ixO^S,idir,idir,kdir) = 4.0d0 * gamma_ij(ixO^S,idir,idir) * dpsi(ixO^S,kdir)
       end do
    end do

    if (present(dalp_out   )) dalp_out   = dalp
    if (present(dpsi_out   )) dpsi_out   = dpsi
    if (present(dbeta_out  )) dbeta_out  = dbeta
    if (present(D_beta_out )) D_beta_out = D_beta
    if (present(Dgamma_out )) Dgamma_out = Dgamma

    if (present(K_ij)) then
       ! extrinsic curvature K_{ij}
       K_ij(ixO^S,:,:) = 0.0d0
       select case (coordinate)
       case (cartesian)
          K_ij(ixO^S,1,1) = gamma_ij(ixO^S,1,1)/(3.0d0*m(ixO^S,alp_)) &
                       * ( 2.0d0*dbeta(ixO^S,1,1) & 
                         {^NOONED - dbeta(ixO^S,2,2) } &
                         {^IFTHREED - dbeta(ixO^S,3,3) } )
   
          K_ij(ixO^S,2,2) = gamma_ij(ixO^S,2,2)/(3.0d0*m(ixO^S,alp_)) &
                       * ( - dbeta(ixO^S,1,1) &
                         {^NOONED + 2.0d0*dbeta(ixO^S,2,2) } &
                         {^IFTHREED - dbeta(ixO^S,3,3) } )
   
          K_ij(ixO^S,3,3) = gamma_ij(ixO^S,3,3)/(3.0d0*m(ixO^S,alp_)) &
                       * ( - dbeta(ixO^S,1,1) &
                         {^NOONED - dbeta(ixO^S,2,2) } &
                         {^IFTHREED + 2.0d0 * dbeta(ixO^S,3,3) } )
       case (cylindrical)
          K_ij(ixO^S,1,1) = gamma_ij(ixO^S,1,1)/(3.0d0*m(ixO^S,alp_)) &
                       * ( 2.0d0*dbeta(ixO^S,1,1) - m(ixO^S,beta(1))/x(ixO^S,r_) & 
                         {^NOONED - dbeta(ixO^S,2,2) } &
                         {^IFTHREED - dbeta(ixO^S,3,3) } )
   
          K_ij(ixO^S,2,2) = gamma_ij(ixO^S,2,2)/(3.0d0*m(ixO^S,alp_)) &
                       * ( - dbeta(ixO^S,1,1) - m(ixO^S,beta(1))/x(ixO^S,r_) &
                         {^NOONED + 2.0d0 * dbeta(ixO^S,2,2) } &
                         {^IFTHREED - dbeta(ixO^S,3,3) } )
   
          K_ij(ixO^S,3,3) = gamma_ij(ixO^S,3,3)/(3.0d0*m(ixO^S,alp_)) &
                       * ( - dbeta(ixO^S,1,1) + 2.0d0 * m(ixO^S,beta(1))/x(ixO^S,r_) &
                         {^NOONED - dbeta(ixO^S,2,2) } &
                         {^IFTHREED + 2.0d0 * dbeta(ixO^S,3,3) } )
       case (spherical)
          {^NOONED
          cot_theta(ixO^S) = dcos(x(ixO^S,theta_))/dsin(x(ixO^S,theta_))
          }

          K_ij(ixO^S,1,1) = gamma_ij(ixO^S,1,1)/(3.0d0*m(ixO^S,alp_)) &
                       * ( 2.0d0*dbeta(ixO^S,1,1) - 2.0d0*m(ixO^S,beta(1))/x(ixO^S,r_) & 
                         {^NOONED - dbeta(ixO^S,2,2) - m(ixO^S,beta(2))*cot_theta(ixO^S) } &
                         {^IFTHREED - dbeta(ixO^S,3,3) } )
   
          K_ij(ixO^S,2,2) = gamma_ij(ixO^S,2,2)/(3.0d0*m(ixO^S,alp_)) &
                       * ( -dbeta(ixO^S,1,1) + m(ixO^S,beta(1))/x(ixO^S,r_) &
                         {^NOONED + 2.0d0*dbeta(ixO^S,2,2) - m(ixO^S,beta(2))*cot_theta(ixO^S) } &
                         {^IFTHREED - dbeta(ixO^S,3,3) } )
   
          K_ij(ixO^S,3,3) = gamma_ij(ixO^S,3,3)/(3.0d0*m(ixO^S,alp_)) &
                       * ( -dbeta(ixO^S,1,1) + m(ixO^S,beta(1))/x(ixO^S,r_) &
                         {^NOONED - dbeta(ixO^S,2,2) + 2.0d0*m(ixO^S,beta(2))*cot_theta(ixO^S) } &
                         {^IFTHREED + 2.0d0 * dbeta(ixO^S,3,3) } )
       end select
       {^NOONED
       K_ij(ixO^S,1,2) = ( dbeta(ixO^S,1,2)*gamma_ij(ixO^S,1,1) + dbeta(ixO^S,2,1)*gamma_ij(ixO^S,2,2) ) / (2.0d0*m(ixO^S,alp_))
       K_ij(ixO^S,1,3) = ( dbeta(ixO^S,3,1)*gamma_ij(ixO^S,3,3) + dbeta(ixO^S,1,3)*gamma_ij(ixO^S,1,1) ) / (2.0d0*m(ixO^S,alp_))
       K_ij(ixO^S,2,3) = ( dbeta(ixO^S,3,2)*gamma_ij(ixO^S,3,3) + dbeta(ixO^S,2,3)*gamma_ij(ixO^S,2,2) ) / (2.0d0*m(ixO^S,alp_))
       ! K_ji=K_ij
       do idir=1,2
          do jdir=idir+1,3
             K_ij(ixO^S,jdir,idir) = K_ij(ixO^S,idir,jdir)
          end do
       end do
       }
    end if

    end associate
  end subroutine get_metric_dervitives

  ! ------------------------------

end module mod_geometry
