! TODO: divide into smaller modules(?), make sure these don't do anything by
! default.
module mod_usr_dummy

  implicit none
  public

contains

  !> special boundary types, user defined user must assign conservative
  !> variables in bounderies
  subroutine specialbound_usr(qt,ixI^L,ixO^L,iw,iB,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, iw, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    call mpistop("specialbound not defined")
  end subroutine specialbound_usr

  !> internal boundary, user defined
  !
  !> This subroutine can be used to artificially overwrite ALL conservative 
  !> variables in a user-selected region of the mesh, and thereby act as
  !> an internal boundary region. It is called just before external (ghost cell)
  !> boundary regions will be set by the BC selection. Here, you could e.g. 
  !> want to introduce an extra variable (nwextra, to be distinguished from nwaux)
  !> which can be used to identify the internal boundary region location.
  !> Its effect should always be local as it acts on the mesh.
  subroutine bc_int(level,qt,ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    call mpistop("bc_int not defined")
  end subroutine bc_int

  subroutine printlog_special()
    use mod_global_parameters

    call mpistop("special log file undefined")
  end subroutine printlog_special

  !> this subroutine is ONLY to be used for computing auxiliary variables
  !> which happen to be non-local (like div v), and are in no way used for
  !> flux computations. As auxiliaries, they are also not advanced
  subroutine process_grid_usr(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    !-----------------------------------------------------------------------------

  end subroutine process_grid_usr

  !> This subroutine is called at the beginning of each time step 
  !> by each processor. No communication is specified, so the user
  !> has to implement MPI routines if information has to be shared
  subroutine process_global_usr(iit,qt)
    use mod_global_parameters
    integer, intent(in)          :: iit
    double precision, intent(in) :: qt

  end subroutine process_global_usr

  !> this subroutine can be used in convert, to add auxiliary variables to the
  !> converted output file, for further analysis using tecplot, paraview, ....
  !> these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  !> the array normconv can be filled in the (nw+1:nw+nwauxio) range with 
  !> corresponding normalization values (default value 1)
  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision             :: w(ixI^S,nw+nwauxio)
    double precision             :: normconv(0:nw+nwauxio)

    call mpistop("special output file undefined")

    ! Example: assuming nwauxio=1 at convert stage and desire to see -w(1)
    ! w(ixO^S,nw+1)=-w(ixO^S,1)
  end subroutine specialvar_output

  !> newly added variables to be concatenated with the primnames/wnames string
  subroutine specialvarnames_output()
    use mod_global_parameters

    call mpistop("special wnames and primnames undefined")

    !> Example : as above in specialvar_output, assuming relativistic HD here...
    !> primnames= TRIM(primnames)//' '//'-rho'
    !> wnames=TRIM(wnames)//' '//'-d'
  end subroutine specialvarnames_output

  !> Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
  !> iw=iwmin...iwmax.  wCT is at time qCT
  subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)
  end subroutine specialsource

  !> Limit "dt" further if necessary, e.g. due to the special source terms.
  !> The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
  !> module have already been called.
  subroutine getdt_special(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw), dtnew

    dtnew=bigdouble
  end subroutine getdt_special

  !> Set the "eta" array for resistive MHD based on w or the
  !> "current" variable which has components between idirmin and 3.
  subroutine specialeta(w,ixI^L,ixO^L,idirmin,x,current,eta)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L, idirmin
    double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision             :: current(ixG^T,7-2*ndir:3), eta(ixG^T)

    call mpistop("specialeta is not defined")
  end subroutine specialeta

  !> Enforce additional refinement or coarsening
  !> One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
  !> you must set consistent values for integers refine/coarsen:
  !> refine = -1 enforce to not refine
  !> refine =  0 doesn't enforce anything
  !> refine =  1 enforce refinement
  !> coarsen = -1 enforce to not coarsen
  !> coarsen =  0 doesn't enforce anything
  !> coarsen =  1 enforce coarsen
  subroutine specialrefine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters
    integer, intent(in)          :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout)       :: refine, coarsen

    ! e.g. refine for negative first coordinate x < 0 as
    ! if (any(x(ix^S,1) < zero)) refine=1
  end subroutine specialrefine_grid

  !> this is the place to compute a local auxiliary variable to be used
  !> as refinement criterion for the Lohner error estimator only
  !>  -->it is then requiring and iflag>nw
  !> note that ixO=ixI=ixG, hence the term local (gradients need special attention!)
  subroutine specialvarforerrest(ixI^L,ixO^L,iflag,w,var)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L,ixO^L,iflag
    double precision, intent(in)  :: w(ixI^S,1:nw)
    double precision, intent(out) :: var(ixG^T)

    if (iflag > nw) call mpistop('iflag > nw, make change in parfile or in user file')
    var(ixI^S) = zero
  end subroutine specialvarforerrest

  !> Here one can add a steady (time-independent) potential background field
  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixG^T,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    call mpistop('specialset_B0 not implemented')
  end subroutine specialset_B0

  subroutine specialsource_impl(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw), wCT(ixI^S,1:nw)
  end subroutine specialsource_impl

  subroutine getdt_impl(w,ixG^L,ix^L,dtnew,dx^D,x)
    use mod_global_parameters
    integer, intent(in)             :: ixG^L,ix^L
    double precision, intent(in)    :: dx^D, x(ixG^S,1:ndim)
    !> note that depending on strictsmall etc, w values may change
    double precision, intent(inout) :: w(ixG^S,1:nw), dtnew

    dtnew = bigdouble
  end subroutine getdt_impl

  !> to allow use defined global process before main loop evolution
  subroutine special_process_usr()
    use mod_global_parameters

    ! call mpistop("special_process_usr not implemented")
    ! use mod_input_output

    ! {#IFDEF MAGNETOFRICTION
    ! if(itmaxmf>0) then
    !    time_in=MPI_WTIME()
    !    call magnetofriction
    !    if(mype==0) write(*,*) 'Magnetofriction phase took : ',MPI_WTIME()-time_in,' sec'
    ! endif
    ! }
    ! !> write transformed data file
    if(nwtf>0 .and. neqpartf>0) then
       call mpistop("special_process_usr not implemented")
    !    call write_snapshot_tf
    end if
  end subroutine special_process_usr

  !> regenerate w and eqpar arrays to output into *tf.dat
  subroutine transformw_usr(w,wtf,eqpar_tf,ixI^L,ixO^L)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw)
    double precision, intent(out) :: wtf(ixI^S,1:nwtf)
    double precision, intent(out) :: eqpar_tf(neqpartf)
    double precision              :: gamma_usr
    integer                       :: iwup_,diw,e_usr
    !-----------------------------------------------------------------------------

    {^IFMHD
    gamma_usr=5.d0/3.d0
    e_usr=m^NC_+1
    wtf(ixO^S,1:m^NC_)=w(ixO^S,1:m^NC_)
    wtf(ixO^S,e_usr)=w(ixO^S,rho_)*1.d0/(gamma_usr-one)+&
         half*((^C&w(ixO^S,m^C_)**2.0d0+)/w(ixO^S,rho_)+(^C&w(ixO^S,b^C_)**2.0d0+))
    wtf(ixO^S,e_usr+1:e_usr+^NC)=w(ixO^S,b1_:b^NC_)
    iwup_=b^NC_
    if(iwup_<nw .and. nw<nwtf) then
       diw=nw-iwup_
       wtf(ixO^S,e_usr+^NC+1:e_usr+^NC+diw)=w(ixO^S,iwup_+1:iwup_+diw)
    endif
    eqpar_tf(1:4)=eqpar(1:4)
    }
  end subroutine transformw_usr

  !> use different tolerance in special regions for AMR to
  !> reduce/increase resolution there where nothing/something interesting happens.
  subroutine special_tolerance(wlocal,xlocal,tolerance,qt)
    use mod_global_parameters
    double precision, intent(in)    :: wlocal(1:nw),xlocal(1:ndim),qt
    double precision, intent(inout) :: tolerance
    double precision                :: bczone^D,addtol,tol_add
  end subroutine special_tolerance

  !> Allow user to use their own data-postprocessing procedures
  subroutine userspecialconvert(qunitconvert)
    use mod_global_parameters
    integer, intent(in) :: qunitconvert
    character(len=20)   :: userconvert_type
  end subroutine userspecialconvert

  subroutine fixp_usr(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(inout)    :: w(ixI^S,1:nw)
    double precision, intent(in)       :: x(ixI^S,1:ndim)
  end subroutine fixp_usr

  !> flag=-1 : Treat all cells active, omit deactivation (onentry, default)
  !> flag=0  : Treat as normal domain
  !> flag=1  : Treat as passive, but reduce by safety belt
  !> flag=2  : Always treat as passive
  subroutine flag_grid_usr(qt,ixG^L,ixO^L,w,x,flag)
    use mod_global_parameters
    integer, intent(in)             :: ixG^L, ixO^L
    integer, intent(inout)          :: flag
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision, intent(in)    :: x(ixG^S,1:ndim)
  end subroutine flag_grid_usr

end module mod_usr_dummy
