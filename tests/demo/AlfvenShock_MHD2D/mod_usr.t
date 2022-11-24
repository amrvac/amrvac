!============================================================================2
! MHD Alfven shock setup
!=============================================================================
module mod_usr
  use mod_mhd
  implicit none

  ! this is input
  double precision :: delta, betab, thetab, xshock, scalefactor
  integer :: type
  logical :: shiftframe
  ! this is computed from input
  double precision :: thetabrad,Bxb,Byb,Bxa,Bya,vxb,vyb,vxa,vya,pa,thetaarad,thetaa
  
  type t_cubic_solution
    logical :: allreal, coincident
    double precision :: x1,x2,x3
    double complex :: z2
    double complex :: z3
  end type t_cubic_solution

contains

  ! this solves a cubic of form x^3+a1 x^2+ a2 x+ a3 = 0 according to Cardano's formulae
  pure type(t_cubic_solution) function solve(a1, a2, a3) result(res)
        double precision, intent(in)    :: a1, a2, a3
        double precision                :: s, t, q, r, rpart, ipart, temp, theta, root
        double precision                :: minx,maxx,xx1,xx2,xx3

        q      = (3.0d0*a2-a1**2)/9.0d0
        r      = (-2.0d0*a1**3+9.0d0*a1*a2-27.0d0*a3)/54.0d0
       
        if (q**3+r**2==0.0d0) then 
           res%coincident=.true.
           res%allreal=.true.
           ! one real and two identical real roots
           temp    =  r 
           s       = dsign(1.0d0, temp) * dabs(temp)**(1.0/3.0)
           temp    =  r 
           t       = dsign(1.0d0, temp) * dabs(temp)**(1.0/3.0)
           res%x1  = s + t - a1/3.0d0
           rpart   = -(s+t)/2.0d0 - a1/3.0d0
           res%x2  = rpart
           res%x3  = rpart
           res%allreal=.false.
        endif
        if (q**3+r**2<0.0d0) then
           ! three real and distinct roots
           root=dsqrt(-q)
           theta=dacos(r/(dabs(q)*root))
           xx1=2.0d0*root*dcos(theta/3.0d0)-a1/3.0d0
           xx2=2.0d0*root*dcos((theta-2.0d0*dpi)/3.0d0)-a1/3.0d0
           minx=min(xx1,xx2)
           maxx=max(xx1,xx2)
           xx3=2.0d0*root*dcos((theta-4.0d0*dpi)/3.0d0)-a1/3.0d0
           if(xx3>maxx)then
             xx2=maxx
             maxx=xx3
           else
              if(xx3>minx)then
                 xx2=xx3
              else
                 minx=xx3
              endif
           endif 
           res%x1=minx
           res%x2=xx2
           res%x3=maxx
           res%allreal=.true.
           res%coincident=.false.
        else
           ! one real and two complex roots
           temp    =  r + dsqrt(q**3 + r**2)
           s       = dsign(1.0d0, temp) * dabs(temp)**(1.0/3.0)
           temp    =  r - dsqrt(q**3 + r**2)
           t       = dsign(1.0d0, temp) * dabs(temp)**(1.0/3.0)
           res%x1  = s + t - a1/3.0d0
           rpart   = -(s+t)/2.0d0 - a1/3.0d0
           ipart   = (dsqrt(3.0d0)/2.0d0)*(s-t) 
           res%z2  = dcmplx( rpart, ipart )
           res%z3  = dcmplx( rpart, -ipart )
           res%allreal=.false.
           res%coincident=.false.
        endif

  end function solve

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ delta, betab, thetab, xshock, type, shiftframe, scalefactor

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do
 
    thetabrad=thetab*dpi/180.0d0

  end subroutine usr_params_read


  subroutine usr_init()

    call usr_params_read(par_files)

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 
    usr_init_vector_potential => initvecpot_usr

    call set_coordinate_system("Cartesian")

    call mhd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    character(len=20) :: printsettingformat
    type(t_cubic_solution) :: p
    double precision:: vAb2, vSb2, cos2theb, sin2theb, aval, bval,cval, dval, eval, a1, a2, a3
    double precision:: MAb2, MAa2, MArat, Bb, deltamax

    printsettingformat='(1x,A50,ES15.7,A7)'

    mhd_gamma=1.66666667d0
    mhd_eta=zero
    mhd_etah=zero
    deltamax=(mhd_gamma+1.0d0)/(mhd_gamma-1.0d0)

    Bb=dsqrt(2.0d0/(mhd_gamma*betab))
    Bxb=-Bb*dcos(thetabrad)
    Byb=Bb*dsin(thetabrad)
    if(mype==0) then
      write(*,*) "MHD alfven shock setup:"
      write(*,printsettingformat) "density contrast ",delta," input"
      if(delta>deltamax) call mpistop("density contrast too high")
      if(delta<1.0d0)    call mpistop("density contrast too low")
      write(*,printsettingformat) "plasma beta before shock ",betab," input"
      write(*,printsettingformat) "angle theta in degrees before shock ",thetab," input"
      write(*,printsettingformat) "angle theta in radians before shock ",thetabrad," input"
      write(*,printsettingformat) "Bx before shock ",Bxb," output"
      write(*,printsettingformat) "By before shock ",Byb," output"
      write(*,*) "Computing the three roots of the shock adiabat"
    endif

    vAb2=2.0d0/(mhd_gamma*betab)
    vSb2=1.0d0
    cos2theb=(dcos(thetabrad))**2
    sin2theb=1.0d0-cos2theb
    aval=delta*vAb2*cos2theb
    bval=2.0d0*delta*vSb2/(delta+1.0d0-mhd_gamma*(delta-1.0d0))
    cval=delta*sin2theb*vAb2
    dval=(2.0d0*delta-mhd_gamma*(delta-1.0d0))/(delta+1.0d0-mhd_gamma*(delta-1.0d0))
    eval=delta*vAb2*cos2theb
    a1=-bval-2.0d0*aval-cval*dval
    a2=2.0d0*aval*bval+aval**2+cval*eval
    a3=-aval**2*bval
    p=solve(a1,a2,a3)
    if(mype==0) write(*,*)"real solution or not: ALL REAL",p%allreal," COINCIDENT ",p%coincident
    if (p%allreal.and..not.p%coincident) then
         if(mype==0)then
            write(*,*)"real solutions:",p%x1,p%x2,p%x3
            write(*,*)"check zero of",p%x1**3+a1*p%x1**2+a2*p%x1+a3
            write(*,*)"check zero of",p%x2**3+a1*p%x2**2+a2*p%x2+a3
            write(*,*)"check zero of",p%x3**3+a1*p%x3**2+a2*p%x3+a3
            write(*,*) "selecting the root number:", type
         endif
         select case(type)
            case(1)
              MAb2=p%x1/Bxb**2
              MAa2=p%x1/(delta*Bxb**2)
              MArat=(MAb2-1.0d0)/(MAa2-1.0d0)
              vxb=-dsqrt(p%x1)
              if(mype==0)then
                  write(*,*)"check neg of M_A^2-1=",MAb2-1.0d0,MAa2-1.0d0
                  write(*,*)"this is the SLOW root"
                  if(MArat<1.0d0)then
                     write(*,*)'OK: bend towards shock normal for this slow root'
                  else
                     call mpistop("error for slow root")
                  endif
              endif
            case(3)
              MAb2=p%x3/Bxb**2
              MAa2=p%x3/(delta*Bxb**2)
              MArat=(MAb2-1.0d0)/(MAa2-1.0d0)
              vxb=-dsqrt(p%x3)
              if(mype==0)then
                  write(*,*)"check pos of M_A^2-1=",MAb2-1.0d0,MAa2-1.0d0
                  write(*,*)"this is the FAST root"
                  if(MArat>1.0d0)then
                     write(*,*)'OK: bend away from shock normal for this fast root'
                  else
                     call mpistop("error for fast root")
                  endif
              endif
            case(2)
              MAb2=p%x2/Bxb**2
              MAa2=p%x2/(delta*Bxb**2)
              MArat=(MAb2-1.0d0)/(MAa2-1.0d0)
              vxb=-dsqrt(p%x2)
              if(mype==0)then
                  write(*,*)"check neg to pos change of M_A^2-1=",MAb2-1.0d0,MAa2-1.0d0
                  write(*,*)"this is the ALFVEN-INTERMEDIATE root"
                  if(MArat<0.0d0)then
                    write(*,*)'OK: flips over shock normal for this alfven root'
                  else
                     call mpistop("error for alfven root")
                  endif
              endif
             case default
                 call mpistop("no type selected")
         end select
         vyb=-vxb*dtan(thetabrad)
         if(mype==0)then
            write(*,printsettingformat) "vx before shock ",vxb," computed"
            write(*,printsettingformat) "vy before shock ",vyb," computed"
         endif
         vxa=vxb/delta
         if(mype==0) write(*,printsettingformat) "vx after shock ",vxa," output"
         vya=-vxa*dtan(thetabrad)*MArat
         if(mype==0) write(*,printsettingformat) "vy after shock ",vya," output"
         Bxa=Bxb
         Bya=Byb*MArat
         thetaarad=datan(-Bya/Bxa)
         thetaa=thetaarad*180.0/dpi
         pa=(Byb**2/2.0d0)*(1.0d0-MArat**2)+1.0d0/mhd_gamma+vxb**2*(delta-1.0d0)/delta
         if(mype==0)then
            write(*,printsettingformat) "Bx before shock ",Bxb," input"
            write(*,printsettingformat) "By before shock ",Byb," input"
            write(*,printsettingformat) "Bx after shock ",Bxa," output"
            write(*,printsettingformat) "By after shock ",Bya," output"
            write(*,*) "hence angle changes from   pre-shock upstream (in degrees):", thetab
            write(*,*) "hence angle changes to  post-shock downstream (in degrees):", thetaa
            write(*,printsettingformat) "pressure p after shock ",pa," output"
            write(*,*) "check on entropy ",pa," must be bigger than ",(delta**mhd_gamma)/mhd_gamma
         endif
         if (pa<(delta**mhd_gamma)/mhd_gamma) then
             call mpistop('no entropy condition fulfilled')
         endif 
    else
       if(mype==0) write(*,*)"real and complex solutions:",p%x1,p%z2,p%z3
       call mpistop('no 3 real solutions')
    endif

    if(shiftframe)then
       if(mype==0) write(*,*)"shifting to the upstream rest frame"
       vxa=vxa-vxb
       vya=vya-vyb
       vxb=0.0d0
       vyb=0.0d0
    endif
 
    if(usr_filename/='')then
       if(mype==0) write(*,*)"inserting picture in the upstream state:",usr_filename
    endif
 
  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    use mod_global_parameters
    use mod_init_datafromfile, only: read_data_set
    integer, intent(in) :: ixG^L,ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision :: val(ixG^S,1)

    logical, save :: first=.true.

    if(usr_filename/='')then
       call read_data_set(ixG^L,ix^L,x,val,1)
       ! the data of the Alfven picture is on the range x [0,1] and y [0,1.5]
       ! rescale: original data between 0 and 256
       val(ix^S,1)=val(ix^S,1)*scalefactor
       where(x(ix^S,1)<0.0d0)
         val(ix^S,1)=0.0d0
       endwhere
       where(x(ix^S,1)>1.0d0)
         val(ix^S,1)=0.0d0
       endwhere
       ! to ensure periodicity
       where(x(ix^S,2)<0.0d0)
         val(ix^S,1)=0.0d0
       endwhere
       where(x(ix^S,2)>1.5d0)
         val(ix^S,1)=0.0d0
       endwhere
    endif

    ! set up the shock positioned at xshock
    where(x(ix^S,1)<xshock)
        w(ix^S,rho_)  =delta
        w(ix^S,mom(1))=vxa
        w(ix^S,mom(2))=vya
        w(ix^S,p_)    =pa
    elsewhere
        w(ix^S,rho_)  =1.0d0+val(ix^S,1)
        w(ix^S,mom(1))=vxb
        w(ix^S,mom(2))=vyb
        w(ix^S,p_)    =1.0d0/mhd_gamma
    endwhere
  
    if(stagger_grid)then
       call b_from_vector_potential(block%ixGs^L,ixG^L,ix^L,block%ws,x)
       call mhd_face_to_center(ix^L,block)
    else
       where(x(ix^S,1)<xshock)
           w(ix^S,mag(1))=Bxa
           w(ix^S,mag(2))=Bya
       elsewhere
           w(ix^S,mag(1))=Bxb
           w(ix^S,mag(2))=Byb
       endwhere
    endif

    if(mype==0.and.first)then
       write(*,*)'Doing alfven shock challenge, ideal MHD, shock at xshock=',xshock,' gamma=',mhd_gamma
       write(*,*)'plasma beta after shock (downstream) is=',2.0d0*pa/(Bxa**2+Bya**2)
       write(*,*)'plasma beta before shock (upstream)  is=',2.0d0/(mhd_gamma*(Bxb**2+Byb**2))
       first=.false.
    endif

    call mhd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine initonegrid_usr

  subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)
    ! initialize the vectorpotential on the corners
    ! used by b_from_vectorpotential()
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixC^L,idir
    double precision, intent(in)       :: xC(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)

    if (idir.eq.3) then
      where(xC(ixC^S,1)<xshock)
        A(ixC^S)=Bxa*xC(ixC^S,2)-Bya*(xC(ixC^S,1)-xshock)
      elsewhere
        A(ixC^S)=Bxb*xC(ixC^S,2)-Byb*(xC(ixC^S,1)-xshock)
      endwhere
    else 
      A(ixC^S) = zero
    end if

  end subroutine initvecpot_usr

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

    double precision                   :: wlocal(ixI^S,nw)
    double precision                   :: tmp(ixI^S) 

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    call mhd_get_pthermal(wlocal,x,ixI^L,ixO^L,tmp)
    ! output the temperature p/rho
    w(ixO^S,nw+1)=tmp(ixO^S)/wlocal(ixO^S,rho_)
    ! output the plasma beta p*2/B**2
    w(ixO^S,nw+2)=tmp(ixO^S)*two/(wlocal(ixO^S,mag(1))**2+wlocal(ixO^S,mag(2))**2)
    ! output the z-component of electric field
    w(ixO^S,nw+3)=(-wlocal(ixO^S,mom(1))*wlocal(ixO^S,mag(2))+wlocal(ixO^S,mom(2))*wlocal(ixO^S,mag(1)))/wlocal(ixO^S,rho_)
    
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='Te beta Ez'

  end subroutine specialvarnames_output

end module mod_usr
