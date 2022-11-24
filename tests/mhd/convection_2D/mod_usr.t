module mod_usr
  use mod_mhd
  use mod_viscosity, only: vc_mu
  implicit none
  double precision :: usr_grav,bstr,temptop

contains

  subroutine usr_init()
    call set_coordinate_system("Cartesian")

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_gravity         => gravity
    usr_refine_grid     => specialrefine_grid
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output 

    call mhd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr
    integer:: mpoly
    double precision:: zeta0,rhat,sigma,zz0,qchand
    double precision:: gamma, qchi, qmpoly, eta2

    ! This setup relates to Hurlburt and Toomre
    ! and calculates the eqpar array from problem-specific input parameters
    ! ---------------------------------------------------------------------
    !
    ! problem specific input parameters for
    ! Hurlburt and Toomre magnetoconvection problem
    !
    ! we need (apart from the aspect ratio which is set in par-file):
    ! i)equilibrium parameters: (1) ratio of specific heats (gamma)
    ! ------------------------- (2) the polytropic index mpoly (m)
    !                            --> sets gravity g
    !                           (3) the density contrast chi
    !                            --> sets top temperature zz0 (z0)
    !                           (4) the Chandrasekhar number qchand (Q)
    !                            --> initial field strength B (need nu and eta)
    ! ii)parameters setting the dissipative coefficients:
    ! ---------------------------------------------------
    !    (1) Prandtl number sigma (viscous/thermal conduction)
    !    (2) the Rayleigh number at half-depth rhat (degree of instability)
    !    (3) the magnetic Prandtl number zeta0
    !      (magnetic diffusion/thermal conduction)
    !
    ! hence fixing   i) gamma, mpoly, chi, qchand
    !               ii) sigma, rhat, zeta0
    !
    ! allows calculation of all equation parameters, namely
    ! the general purpose ones:
    ! ------------------------
    !  mhd_gamma mhd_eta tc_k_para eqpar(grav1.2_) vc_mu
    !
    ! and the problem specific one:
    ! ----------------------------
    !  temptop for setting the top temperature
    !  bstr for setting the B field
    !
    
    ! equilibrium parameters
    gamma=1.66666667d0
    mpoly=1
    qchi=11.0d0
    qchand=72.0d0
    
    qmpoly=dble(mpoly)
    zz0=one/(qchi**(one/qmpoly)-one)
    
    if (mype==0) then
      write(*,*)'gamma, mpoly, qchi,qchand'
      write(*,*)gamma, mpoly, qchi, qchand
      write(*,*)'Deduced top dimensionless temperature z0=',zz0
    endif
    
    temptop=zz0
    mhd_gamma=gamma
    usr_grav=-(qmpoly+one)
    
    ! dissipative parameters
    sigma=1.0d0
    rhat=1.0d5
    zeta0=0.25d0
    
    if (mype==0) then
      write(*,*) 'sigma,rhat,zeta0:'
      write(*,*) sigma,rhat,zeta0
    endif
    
    eta2=(qmpoly+one)*(gamma/(gamma-1)-(qmpoly+one))*((gamma-one)/gamma)*&
         ((zz0+half)**(two*qmpoly-one)/zz0**(two*qmpoly))*&
         (zeta0**2/(sigma*rhat))
    
    if(eta2<=smalldouble)then
       if(mype==0) write(*,*) 'eta2=',eta2
       call mpistop("Negative or too small value for eta**2")
    endif
    
    mhd_eta=dsqrt(eta2)
    vc_mu=mhd_eta*sigma/zeta0
    tc_fl%tc_k_para=(gamma/(gamma-one))*mhd_eta/zeta0
    bstr=dsqrt(qchand*vc_mu*mhd_eta)
    
    if (mype==0) then
      write(*,*)'dimensionless values for dissipative coefficients:'
      write(*,*)'resistivity          eta=',mhd_eta
      write(*,*)'viscosity             mu=',vc_mu
      write(*,*)'thermal conduction tc_k_para=',tc_fl%tc_k_para
      write(*,*)'dimensionless magnetic field strength:',bstr
    endif

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    ! initialize one grid
    integer, intent(in) :: ixG^L,ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision :: dvx,dvy,dvz,nkx,nky,nkz,zz0,qmpoly
    logical, save :: first=.true.

    ! velocity perturbation parameters
    dvx=1.0d-4
    dvy=1.0d-4
    nkx=two*dpi/3.0d0
    nky=two*dpi
    {^IFTHREED
    dvz=1.0d-4
    nkz=two*two*dpi
    }
    
    if(first)then
      if(mype==0) then
         write(*,*)'Simulating mhd convection: Hurlburt-Toomre'
         write(*,*) 'velocity perturbation: dvx,dvy,nkx,nky'
         write(*,*) dvx,dvy,nkx,nky
         {^IFTHREED
         write(*,*) '3D setup!'
         write(*,*) dvz,nkz
         }
      endif
      first=.false.
    endif
    
    zz0=temptop
    qmpoly=-one-usr_grav
    
    w(ix^S,rho_)=((zz0+one-x(ix^S,2))/zz0)**qmpoly
    w(ix^S,mag(1))=zero
    w(ix^S,mag(2))=bstr
    {^IFTHREED
    w(ix^S,mag(3))=zero
    }
    
    w(ix^S,p_)= zz0*(((zz0+one-x(ix^S,2))/zz0)**(qmpoly+one))
    
    w(ix^S,mom(1))=dvx*dsin(x(ix^S,1)*nkx)*dsin(x(ix^S,2)*nky){^IFTHREED *dsin(x(ix^S,3)*nkz)}
    w(ix^S,mom(2))=dvy*dsin(x(ix^S,1)*nkx)*dsin(x(ix^S,2)*nky){^IFTHREED *dsin(x(ix^S,3)*nkz)}
    {^IFTHREED
    w(ix^S,mom(3))=dvz*dsin(x(ix^S,1)*nkx)*dsin(x(ix^S,2)*nky)*dsin(x(ix^S,3)*nkz)
    }

    call mhd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixG^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixG^L
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision :: tempb(ixG^S)
    double precision :: delydelx,delydelz
    integer :: ix^D,ixIM^L

    select case(iB)
    case(3)
      ! special bottom boundary
      ! ensure fixed temperature gradient in ghost layers
      ! use asymm/symm for B-components (ensuring vertical field)
      ! use divB correction from central difference for vertical field
      ! use symm/asymm for v-components (ensuring no flow-through)
      ! density profile: use symmetry
    
      {^IFTWOD
      ! in nghostcells rows above the bottom boundary: switch to primitive
      ixIMmin2=ixOmax2+1;ixIMmax2=ixOmax2+nghostcells;
      ixIMmin1=ixOmin1;ixIMmax1=ixOmax1;
      call mhd_to_primitive(ixG^L,ixIM^L,w,x)
      do ix2=ixOmax2,ixOmin2,-1
         w(ixOmin1:ixOmax1,ix2,rho_)= w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,rho_)
         w(ixOmin1:ixOmax1,ix2,mom(1)) = w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,mom(1))
         w(ixOmin1:ixOmax1,ix2,mom(2)) =-w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,mom(2))
         w(ixOmin1:ixOmax1,ix2,mag(1)) =-w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,mag(1))
         w(ixOmin1:ixOmax1,ix2,mag(2)) = w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,mag(2))
      enddo
      ! fill temperature array: extrapolate linearly with fixed dT/dy=-1, from bottom row
      do ix2=ixOmin2,ixOmax2
        tempb(ixOmin1:ixOmax1,ix2)=w(ixOmin1:ixOmax1,ixOmax2+1,p_)/w(ixOmin1:ixOmax1,ixOmax2+1,rho_) &
                                 -(x(ixOmin1:ixOmax1,ix2,2)-x(ixOmin1:ixOmax1,ixOmax2+1,2))
      enddo
      ! determine pressure
      w(ixO^S,p_)=w(ixO^S,rho_)*tempb(ixO^S)
      !delydelx=(x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2))/(x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1))
      !do ix2=ixOmax2,ixOmin2,-1
      !   do ix1=ixOmin1+1,ixOmax1-1
      !       w(ix1,ix2,mag(2))=w(ix1,ix2+2,mag(2))+delydelx*(w(ix1+1,ix2+1,mag(1))-w(ix1-1,ix2+1,mag(1)))
      !   enddo
      !enddo
      }
      {^IFTHREED
      ! in nghostcells rows above the bottom boundary: switch to primitive
      ixIMmin2=ixOmax2+1;ixIMmax2=ixOmax2+nghostcells;
      ixIMmin1=ixOmin1;ixIMmax1=ixOmax1;
      ixIMmin3=ixOmin3;ixIMmax3=ixOmax3;
      call mhd_to_primitive(ixG^L,ixIM^L,w,x)
      do ix2=ixOmax2,ixOmin2,-1
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,rho_)= w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,ixOmin3:ixOmax3,rho_)
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mom(1)) = w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,ixOmin3:ixOmax3,mom(1))
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mom(2)) =-w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,ixOmin3:ixOmax3,mom(2))
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mom(3)) = w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,ixOmin3:ixOmax3,mom(3))
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1)) =-w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,ixOmin3:ixOmax3,mag(1))
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(2)) = w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,ixOmin3:ixOmax3,mag(2))
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(3)) =-w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,ixOmin3:ixOmax3,mag(3))
      enddo
      ! fill temperature array: extrapolate linearly with fixed dT/dy=-1, from bottom row
      do ix2=ixOmin2,ixOmax2
        tempb(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,p_) &
                                                  /w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,rho_) &
               -(x(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,2)-x(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,2))
      enddo
      ! determine pressure
      w(ixO^S,p_)=w(ixO^S,rho_)*tempb(ixO^S)
      !delydelx=(x(ixOmin1,ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2)) &
      !        /(x(ixOmin1+1,ixOmin2,ixOmin3,1)-x(ixOmin1,ixOmin2,ixOmin3,1))
      !delydelz=(x(ixOmin1,ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2)) &
      !        /(x(ixOmin1,ixOmin2,ixOmin3+1,3)-x(ixOmin1,ixOmin2,ixOmin3,3))
      !do ix2=ixOmax2,ixOmin2,-1
      !   do ix1=ixOmin1+1,ixOmax1-1
      !   do ix3=ixOmin3+1,ixOmax3-1
      !       w(ix1,ix2,ix3,mag(2))=w(ix1,ix2+2,ix3,mag(2))+delydelx*(w(ix1+1,ix2+1,ix3,mag(1))-w(ix1-1,ix2+1,ix3,mag(1))) &
      !                                              +delydelz*(w(ix1,ix2+1,ix3+1,mag(3))-w(ix1,ix2+1,ix3-1,mag(3)))
      !   enddo
      !   enddo
      !enddo
      }
      ! now reset the inner mesh values to conservative
      call mhd_to_conserved(ixG^L,ixIM^L,w,x)
      ! now switch to conservative in full bottom ghost layer
      call mhd_to_conserved(ixG^L,ixO^L,w,x)
    
    case(4)
      ! special top boundary
      ! ensure the fixed temperature in ghost layers
      ! use asymm/symm for B-components (ensuring vertical field)
      ! use divB correction from central difference for vertical field
      ! use symm/asymm for v-components (ensuring no flow-through)
      ! density profile: use symmetry
    
      {^IFTWOD
      ! in nghostcells rows below top boundary: switch to primitive
      ixIMmin2=ixOmin2-nghostcells;ixIMmax2=ixOmin2-1;
      ixIMmin1=ixOmin1;ixIMmax1=ixOmax1;
      call mhd_to_primitive(ixG^L,ixIM^L,w,x)
      do ix2=ixOmin2,ixOmax2,+1
          w(ixOmin1:ixOmax1,ix2,rho_)= w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,rho_)
          w(ixOmin1:ixOmax1,ix2,mom(1)) = w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,mom(1))
          w(ixOmin1:ixOmax1,ix2,mom(2)) =-w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,mom(2))
          w(ixOmin1:ixOmax1,ix2,mag(1)) =-w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,mag(1))
          w(ixOmin1:ixOmax1,ix2,mag(2)) = w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,mag(2))
      enddo
      w(ixO^S,p_)=w(ixO^S,rho_)*temptop
      !delydelx=(x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2))/(x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1))
      !do ix2=ixOmin2,ixOmax2,+1
      !   do ix1=ixOmin1+1,ixOmax1-1
      !       w(ix1,ix2,mag(2))=w(ix1,ix2-2,mag(2))-delydelx*(w(ix1+1,ix2-1,mag(1))-w(ix1-1,ix2-1,mag(1)))
      !   enddo
      !enddo
      }
      {^IFTHREED
      ! in nghostcells rows below top boundary: switch to primitive
      ixIMmin2=ixOmin2-nghostcells;ixIMmax2=ixOmin2-1;
      ixIMmin1=ixOmin1;ixIMmax1=ixOmax1;
      ixIMmin3=ixOmin3;ixIMmax3=ixOmax3;
      call mhd_to_primitive(ixG^L,ixIM^L,w,x)
      do ix2=ixOmin2,ixOmax2,+1
          w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,rho_)= w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,ixOmin3:ixOmax3,rho_)
          w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mom(1)) = w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,ixOmin3:ixOmax3,mom(1))
          w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mom(2)) =-w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,ixOmin3:ixOmax3,mom(2))
          w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mom(3)) = w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,ixOmin3:ixOmax3,mom(3))
          w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1)) =-w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,ixOmin3:ixOmax3,mag(1))
          w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(2)) = w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,ixOmin3:ixOmax3,mag(2))
          w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(3)) =-w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,ixOmin3:ixOmax3,mag(3))
      enddo
      w(ixO^S,p_)=w(ixO^S,rho_)*temptop
      !delydelx=(x(ixOmin1,ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2)) &
      !         /(x(ixOmin1+1,ixOmin2,ixOmin3,1)-x(ixOmin1,ixOmin2,ixOmin3,1))
      !delydelz=(x(ixOmin1,ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2)) &
      !         /(x(ixOmin1,ixOmin2,ixOmin3+1,3)-x(ixOmin1,ixOmin2,ixOmin3,3))
      !do ix2=ixOmin2,ixOmax2,+1
      !   do ix1=ixOmin1+1,ixOmax1-1
      !   do ix3=ixOmin3+1,ixOmax3-1
      !       w(ix1,ix2,ix3,mag(2))=w(ix1,ix2-2,ix3,mag(2))-delydelx*(w(ix1+1,ix2-1,ix3,mag(1))-w(ix1-1,ix2-1,ix3,mag(1))) &
      !                                              -delydelz*(w(ix1,ix2-1,ix3+1,mag(3))-w(ix1,ix2-1,ix3-1,mag(3)))
      !   enddo
      !   enddo
      !enddo
      }
      ! now reset the inner mesh values to conservative
      call mhd_to_conserved(ixG^L,ixIM^L,w,x)
      ! now switch to conservative in full bottom ghost layer
      call mhd_to_conserved(ixG^L,ixO^L,w,x)
    case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

  ! Calculate gravitational acceleration in each dimension
  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)

    gravity_field(ixO^S,:)=0.d0
    gravity_field(ixO^S,2)=usr_grav

  end subroutine gravity

  subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    integer, intent(in) :: igrid, level, ix^L, ixG^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    if (any(x(ix^S,2)<=xprobmin2+0.05d0)) then
      refine=1
      coarsen=-1
    endif

  end subroutine specialrefine_grid

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    ! this subroutine can be used in convert, to add auxiliary variables to the
    ! converted output file, for further analysis using tecplot, paraview, ....
    ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
    !
    ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
    ! corresponding normalization values (default value 1)
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: tmp(ixI^S),divb(ixI^S),current(ixI^S,7-2*ndir:3)
    integer          :: idirmin

    ! output Te
    call mhd_get_pthermal(w,x,ixI^L,ixO^L,tmp)
    w(ixO^S,nw+1)=tmp(ixO^S)/w(ixO^S,rho_)
    ! output B 
    w(ixO^S,nw+2)=dsqrt(sum(w(ixO^S,mag(:))**2,dim=ndim+1))
    ! output divB1
    call get_divb(w,ixI^L,ixO^L,divb)
    w(ixO^S,nw+3)=divb(ixO^S)
    ! output the plasma beta p*2/B**2
    w(ixO^S,nw+4)=tmp(ixO^S)*two/sum(w(ixO^S,mag(:))**2,dim=ndim+1)
    ! store current
    call get_current(w,ixI^L,ixO^L,idirmin,current)
    w(ixO^S,nw+5)=current(ixO^S,3)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
  character(len=*) :: varnames
  varnames='Te B divB beta j3'

  end subroutine specialvarnames_output

end module mod_usr
