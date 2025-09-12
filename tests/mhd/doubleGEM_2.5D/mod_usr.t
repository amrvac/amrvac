! Double GEM problem: resistive to Hall MHD
module mod_usr
  use mod_mhd
  implicit none
  double precision :: sheetl, rhorat, BB0, llx, lly, psi0bot, psi0top

contains

  subroutine usr_init()
    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 
    usr_refine_grid   => specialrefine_grid
    usr_init_vector_potential=>initvecpot_usr

    call set_coordinate_system("Cartesian_2.5D")
    call mhd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    sheetl=1.0d0
    rhorat=0.1d0
    BB0=1.0d0
    llx=xprobmax1-xprobmin1
    lly=xprobmax2-xprobmin2

    select case(iprob)
     case(1)
      psi0bot=0.1d0
      psi0top=0.1d0
      mhd_eta=1.0d-3
      mhd_etah=dsqrt(rhorat)
     case(10)
      psi0bot=0.0d0
      psi0top=0.0d0
      mhd_eta=1.0d-3
      mhd_etah=dsqrt(rhorat)
     case(2)
      psi0bot=0.1d0
      psi0top=0.1d0
      mhd_eta=1.0d-3
      mhd_etah=0.0d0
     case(20)
      psi0bot=0.0d0
      psi0top=0.0d0
      mhd_eta=1.0d-3
      mhd_etah=0.0d0
     case default
      call mpistop('iprob not given')
    end select

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
  ! initialize one grid
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    
    double precision:: xmid,ysh1,ysh2,fkx,fky
    logical, save :: first=.true.

    if(mype==0.and.first)then
       write(*,*)'Doing 2.5D double GEM challenge, resistive Hall-MHD'
       write(*,*)'iprob=',iprob
       write(*,*)'resistivity equal to=',mhd_eta
       write(*,*)'hall parameter equal to=',mhd_etah
       write(*,*)'hyperresistivity equal to=',mhd_eta_hyper
       write(*,*)'perturbation amplitudes=',psi0bot,psi0top
       first=.false.
    endif
    
    ! no initial velocity 
    w(ixO^S,mom(1))  =zero
    w(ixO^S,mom(2))  =zero
    w(ixO^S,mom(3))  =zero
    
    xmid=xprobmin1+0.5d0*llx
    ysh1=xprobmin2+0.25d0*lly
    ysh2=xprobmin2+0.75d0*lly
    fkx=two*dpi/llx
    fky=two*dpi/lly

    if(stagger_grid) then
      call b_from_vector_potential(ixGs^LL,ixI^L,ixO^L,block%ws,x)
      call mhd_face_to_center(ixO^L,block)
    else
      ! set up the 1D equilibrium variation
      w(ixI^S,mag(1))  =BB0*(-one+dtanh((x(ixI^S,2)-ysh1)/sheetl)  &
                                      +dtanh((ysh2-x(ixI^S,2))/sheetl))
      w(ixI^S,mag(2))  =zero
      
      ! add the perturbation
      w(ixI^S,mag(1))= w(ixI^S,mag(1))-psi0bot*fky &
                    *dcos(fkx*(x(ixI^S,1)-xmid))              &
                    *(dsin(fky*(x(ixI^S,2)-ysh1))+two*(x(ixI^S,2)-ysh1)*dcos(fky*(x(ixI^S,2)-ysh1))) &
                    *dexp(-fkx*(x(ixI^S,1)-xmid)**2-fky*(x(ixI^S,2)-ysh1)**2) &
                               +psi0top*fky &
                    *dcos(fkx*(x(ixI^S,1)-xmid))              &
                    *(dsin(fky*(x(ixI^S,2)-ysh2))+two*(x(ixI^S,2)-ysh2)*dcos(fky*(x(ixI^S,2)-ysh2))) &
                    *dexp(-fkx*(x(ixI^S,1)-xmid)**2-fky*(x(ixI^S,2)-ysh2)**2)
      w(ixI^S,mag(2))= w(ixI^S,mag(2))+psi0bot*fkx &
                    *dcos(fky*(x(ixI^S,2)-ysh1))              &
                    *(dsin(fkx*(x(ixI^S,1)-xmid))+two*(x(ixI^S,1)-xmid)*dcos(fkx*(x(ixI^S,1)-xmid))) &
                    *dexp(-fkx*(x(ixI^S,1)-xmid)**2-fky*(x(ixI^S,2)-ysh1)**2) &
                               -psi0top*fkx &
                    *dcos(fky*(x(ixI^S,2)-ysh2))              &
                    *(dsin(fkx*(x(ixI^S,1)-xmid))+two*(x(ixI^S,1)-xmid)*dcos(fkx*(x(ixI^S,1)-xmid))) &
                    *dexp(-fkx*(x(ixI^S,1)-xmid)**2-fky*(x(ixI^S,2)-ysh2)**2)
    end if
      
    w(ixI^S,rho_)=(rhorat+one/(dcosh((x(ixI^S,2)-ysh1)/sheetl)**2) &
                  +one/(dcosh((x(ixI^S,2)-ysh2)/sheetl)**2))
    
    w(ixO^S,mag(3))=zero

    w(ixO^S,p_)=half*BB0**2*w(ixO^S,rho_)
    
    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)
    ! initialize the vectorpotential on the corners
    ! used by b_from_vectorpotential()
    integer, intent(in)                :: ixI^L, ixC^L, idir
    double precision, intent(in)       :: xC(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)

    double precision                   :: xmid,ysh1,ysh2,fkx,fky

    A = zero

    xmid=xprobmin1+0.5d0*llx
    ysh1=xprobmin2+0.25d0*lly
    ysh2=xprobmin2+0.75d0*lly
    fkx=two*dpi/llx
    fky=two*dpi/lly

    if (idir .eq. 3) then
      A(ixC^S) = -BB0*(xC(ixC^S,2) - &
      sheetl*log(cosh((xC(ixC^S,2)-ysh1)/sheetl)) + &
      sheetl*log(cosh((ysh2-xC(ixC^S,2))/sheetl))) &
      + psi0bot*cos(fkx*(xC(ixC^S,1)-( 0.d0*llx)/(16.d0)))*cos(fky*(xC(ixC^S,2)-ysh1))*exp(-fkx*(xC(ixC^S,1)-( 0.d0*llx)/(16.d0))**2-fky*(xC(ixC^S,2)-ysh1)**2)&
      - psi0top*cos(fkx*(xC(ixC^S,1)-( 0.d0*llx)/(16.d0)))*cos(fky*(xC(ixC^S,2)-ysh2))*exp(-fkx*(xC(ixC^S,1)-( 0.d0*llx)/(16.d0))**2-fky*(xC(ixC^S,2)-ysh2)**2)
      !+ psi0bot*cos(fkx*(xC(ixC^S,1)-( 1.d0*llx)/(16.d0)))*cos(fky*(xC(ixC^S,2)-ysh1))*exp(-fkx*(xC(ixC^S,1)-( 1.d0*llx)/(16.d0))**2-fky*(xC(ixC^S,2)-ysh1)**2)&
      !- psi0top*cos(fkx*(xC(ixC^S,1)-( 1.d0*llx)/(16.d0)))*cos(fky*(xC(ixC^S,2)-ysh2))*exp(-fkx*(xC(ixC^S,1)-( 1.d0*llx)/(16.d0))**2-fky*(xC(ixC^S,2)-ysh2)**2)&
      !+ psi0bot*cos(fkx*(xC(ixC^S,1)-( 3.d0*llx)/(16.d0)))*cos(fky*(xC(ixC^S,2)-ysh1))*exp(-fkx*(xC(ixC^S,1)-( 3.d0*llx)/(16.d0))**2-fky*(xC(ixC^S,2)-ysh1)**2)&
      !- psi0top*cos(fkx*(xC(ixC^S,1)-( 3.d0*llx)/(16.d0)))*cos(fky*(xC(ixC^S,2)-ysh2))*exp(-fkx*(xC(ixC^S,1)-( 3.d0*llx)/(16.d0))**2-fky*(xC(ixC^S,2)-ysh2)**2)&
      !+ psi0bot*cos(fkx*(xC(ixC^S,1)-( 5.d0*llx)/(16.d0)))*cos(fky*(xC(ixC^S,2)-ysh1))*exp(-fkx*(xC(ixC^S,1)-( 5.d0*llx)/(16.d0))**2-fky*(xC(ixC^S,2)-ysh1)**2)&
      !- psi0top*cos(fkx*(xC(ixC^S,1)-( 5.d0*llx)/(16.d0)))*cos(fky*(xC(ixC^S,2)-ysh2))*exp(-fkx*(xC(ixC^S,1)-( 5.d0*llx)/(16.d0))**2-fky*(xC(ixC^S,2)-ysh2)**2)&
      !+ psi0bot*cos(fkx*(xC(ixC^S,1)-( 7.d0*llx)/(16.d0)))*cos(fky*(xC(ixC^S,2)-ysh1))*exp(-fkx*(xC(ixC^S,1)-( 7.d0*llx)/(16.d0))**2-fky*(xC(ixC^S,2)-ysh1)**2)&
      !- psi0top*cos(fkx*(xC(ixC^S,1)-( 7.d0*llx)/(16.d0)))*cos(fky*(xC(ixC^S,2)-ysh2))*exp(-fkx*(xC(ixC^S,1)-( 7.d0*llx)/(16.d0))**2-fky*(xC(ixC^S,2)-ysh2)**2)&
      !+ psi0bot*cos(fkx*(xC(ixC^S,1)-(-1.d0*llx)/(16.d0)))*cos(fky*(xC(ixC^S,2)-ysh1))*exp(-fkx*(xC(ixC^S,1)-(-1.d0*llx)/(16.d0))**2-fky*(xC(ixC^S,2)-ysh1)**2)&
      !- psi0top*cos(fkx*(xC(ixC^S,1)-(-1.d0*llx)/(16.d0)))*cos(fky*(xC(ixC^S,2)-ysh2))*exp(-fkx*(xC(ixC^S,1)-(-1.d0*llx)/(16.d0))**2-fky*(xC(ixC^S,2)-ysh2)**2)&
      !+ psi0bot*cos(fkx*(xC(ixC^S,1)-(-3.d0*llx)/(16.d0)))*cos(fky*(xC(ixC^S,2)-ysh1))*exp(-fkx*(xC(ixC^S,1)-(-3.d0*llx)/(16.d0))**2-fky*(xC(ixC^S,2)-ysh1)**2)&
      !- psi0top*cos(fkx*(xC(ixC^S,1)-(-3.d0*llx)/(16.d0)))*cos(fky*(xC(ixC^S,2)-ysh2))*exp(-fkx*(xC(ixC^S,1)-(-3.d0*llx)/(16.d0))**2-fky*(xC(ixC^S,2)-ysh2)**2)&
      !+ psi0bot*cos(fkx*(xC(ixC^S,1)-(-5.d0*llx)/(16.d0)))*cos(fky*(xC(ixC^S,2)-ysh1))*exp(-fkx*(xC(ixC^S,1)-(-5.d0*llx)/(16.d0))**2-fky*(xC(ixC^S,2)-ysh1)**2)&
      !- psi0top*cos(fkx*(xC(ixC^S,1)-(-5.d0*llx)/(16.d0)))*cos(fky*(xC(ixC^S,2)-ysh2))*exp(-fkx*(xC(ixC^S,1)-(-5.d0*llx)/(16.d0))**2-fky*(xC(ixC^S,2)-ysh2)**2)&
      !+ psi0bot*cos(fkx*(xC(ixC^S,1)-(-7.d0*llx)/(16.d0)))*cos(fky*(xC(ixC^S,2)-ysh1))*exp(-fkx*(xC(ixC^S,1)-(-7.d0*llx)/(16.d0))**2-fky*(xC(ixC^S,2)-ysh1)**2)&
      !- psi0top*cos(fkx*(xC(ixC^S,1)-(-7.d0*llx)/(16.d0)))*cos(fky*(xC(ixC^S,2)-ysh2))*exp(-fkx*(xC(ixC^S,1)-(-7.d0*llx)/(16.d0))**2-fky*(xC(ixC^S,2)-ysh2)**2)
      !nasty implementation, do loop for every pinchpoint
    end if

  end subroutine initvecpot_usr

  subroutine specialrefine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    double precision:: ysh1,ysh2,shwidth

    ysh1=xprobmin2+0.25d0*lly
    ysh2=xprobmin2+0.75d0*lly
    shwidth=sheetl

    if(any(dabs(x(ixO^S,2)-ysh1)<shwidth).or.any(dabs(x(ixO^S,2)-ysh2)<shwidth)) then
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

    double precision :: current(ixI^S,7-2*ndir:3)
    double precision :: curlv(ixI^S,7-2*ndir:3),vvec(ixI^S,1:ndir)
    double precision :: wloc(ixI^S,1:nw)
    double precision                   :: gradrho(ixI^S),rho(ixI^S),drho(ixI^S)
    double precision                   :: kk,kk0,grhomax,kk1
    integer :: idirmin,idir,idims,idirmin0

    idirmin0=7-2*ndir
    wloc(ixI^S,1:nw)=w(ixI^S,1:nw)
    call get_current(wloc,ixI^L,ixO^L,idirmin,current)
    w(ixO^S,nw+1)=current(ixO^S,3)

    rho(ixI^S)=w(ixI^S,rho_)
    gradrho(ixO^S)=zero
    do idims=1,ndim
       select case(typegrad)
       case("central")
         call gradient(rho,ixI^L,ixO^L,idims,drho)
       case("limited")
         call gradientL(rho,ixI^L,ixO^L,idims,drho)
       end select
       gradrho(ixO^S)=gradrho(ixO^S)+drho(ixO^S)**2.0d0
    enddo
    gradrho(ixO^S)=dsqrt(gradrho(ixO^S))
    kk=5.0d0
    kk0=0.01d0
    kk1=1.0d0
    grhomax=1000.0d0
    w(ixO^S,nw+2)=dexp(-kk*(gradrho(ixO^S)-kk0*grhomax)/(kk1*grhomax-kk0*grhomax))

    do idir=1,ndir
      vvec(ixI^S,idir)=w(ixI^S,mom(idir))/w(ixI^S,rho_)
    enddo
    call curlvector(vvec,ixI^L,ixO^L,curlv,idirmin,idirmin0,ndir)
    w(ixO^S,nw+3)=curlv(ixO^S,3)
    
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables to be concatenated with the primnames/w_names string
    character(len=*) :: varnames
    varnames='jz schlier curlvz'

  end subroutine specialvarnames_output

end module mod_usr
