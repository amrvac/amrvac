module mod_usr
  use mod_mhd
  use mod_eos, only: eos
  implicit none

  double precision :: rho0=100.d0
  double precision :: p0=1.d0
  double precision :: perturbation=1.d-6
  double precision :: bdir1=2.d0/dsqrt(5.d0)
  double precision :: bdir2=1.d0/dsqrt(5.d0)
  double precision :: bmag=1.d-6
  double precision :: kappa_parallel_usr=150.d0
  double precision :: kappa_perp_ratio_usr=0.1d0
  double precision :: tau_target=1.d-3
  integer :: kx_mode=2
  integer :: ky_mode=1
  logical :: fourier_consistent_flux=.true.

contains

  subroutine usr_init()
    usr_init_one_grid => initonegrid_usr
    usr_print_log     => print_fourier_log
    usr_set_parameters => set_usr_parameters

    call set_coordinate_system("Cartesian_2D")
    call mhd_activate()
  end subroutine usr_init

  subroutine set_usr_parameters()
    integer :: n
    namelist /usr_list/ rho0,p0,perturbation,bdir1,bdir2,bmag, &
         kappa_parallel_usr,kappa_perp_ratio_usr,tau_target, &
         kx_mode,ky_mode,fourier_consistent_flux

    do n=1,size(par_files)
      open(unitpar,file=trim(par_files(n)),status='old')
      read(unitpar,usr_list,end=111)
111   close(unitpar)
    end do
    if(rho0<=zero .or. p0<=zero) call mpistop('rho0 and p0 must be positive')
    if(perturbation<=zero .or. perturbation>=1.d-2) &
      call mpistop('perturbation must be positive and remain in the linear regime')
    if(kx_mode==0 .and. ky_mode==0) call mpistop('Fourier mode cannot be zero')
    if(kappa_parallel_usr<=zero) call mpistop('kappa_parallel_usr must be positive')
    if(kappa_perp_ratio_usr<zero) call mpistop('kappa_perp_ratio_usr cannot be negative')
  end subroutine set_usr_parameters

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    use mod_physics, only: phys_to_conserved
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: phase(ixI^S),temperature(ixI^S)
    double precision :: sinphase(ixI^S)
    double precision :: t0,kvec1,kvec2,kmag,kdotb,kperpmag

    kvec1=2.d0*dpi*dble(kx_mode)
    kvec2=2.d0*dpi*dble(ky_mode)
    kmag=dsqrt(kvec1**2+kvec2**2)
    kdotb=kvec1*bdir1+kvec2*bdir2
    kperpmag=dsqrt(max(kmag**2-kdotb**2,zero))
    t0=p0/rho0
    phase(ixO^S)=kvec1*x(ixO^S,1)+kvec2*x(ixO^S,2)
    sinphase(ixO^S)=dsin(phase(ixO^S))
    temperature(ixO^S)=t0*(one+perturbation*dcos(phase(ixO^S)))

    w(ixO^S,1:nw)=zero
    ! An isobaric entropy mode delays hydrodynamic contamination.  The large
    ! rho0 makes the sound-crossing time much longer than this benchmark.
    w(ixO^S,p_)=p0
    w(ixO^S,rho_)=p0/temperature(ixO^S)
    w(ixO^S,mom(:))=zero
    w(ixO^S,mag(1))=bmag*bdir1
    w(ixO^S,mag(2))=bmag*bdir2
    if(mhd_hyperbolic_tc) then
      w(ixO^S,qpar_)=zero
      w(ixO^S,qperp_)=zero
      if(fourier_consistent_flux) then
        w(ixO^S,qpar_)=kappa_parallel_usr*t0*perturbation*kdotb*sinphase(ixO^S)
        w(ixO^S,qperp_)=-kappa_parallel_usr*kappa_perp_ratio_usr* &
             t0*perturbation*kperpmag*dabs(sinphase(ixO^S))
      end if
    end if

    call phys_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine initonegrid_usr

  subroutine print_fourier_log()
    integer, parameter :: log_unit=123
    integer :: iigrid,igrid
    double precision :: local_sum(8),global_sum(8),dvolume(ixG^T)
    double precision :: wprim(ixG^T,1:nw),phase(ixG^T),temperature(ixG^T)
    double precision :: sinphase(ixG^T),qpar_k(ixG^T),qperp_k(ixG^T)
    double precision :: t0,kvec1,kvec2,kmag,costheta,sintheta
    double precision :: amplitude_T,amplitude_qpar,amplitude_qperp
    double precision :: amplitude_qtotal,mean_temperature,mean_energy,vrms
    double precision :: gamma_pred,kappa_eff
    logical, save :: first_call=.true.
    character(len=std_len) :: filename

    kvec1=2.d0*dpi*dble(kx_mode)
    kvec2=2.d0*dpi*dble(ky_mode)
    kmag=dsqrt(kvec1**2+kvec2**2)
    costheta=(kvec1*bdir1+kvec2*bdir2)/kmag
    costheta=max(-one,min(one,costheta))
    sintheta=dsqrt(max(one-costheta**2,zero))
    t0=p0/rho0
    kappa_eff=kappa_parallel_usr*(costheta**2+kappa_perp_ratio_usr*sintheta**2)
    gamma_pred=eos%gamma_minus_1*kappa_eff*kmag**2/rho0

    local_sum=zero
    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      block=>ps(igrid)
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      if(slab) then
        dvolume(ixM^T)={rnode(rpdx^D_,igrid)|*}
      else
        dvolume(ixM^T)=block%dvolume(ixM^T)
      end if

      wprim(ixG^T,1:nw)=ps(igrid)%w(ixG^T,1:nw)
      call eos%to_primitive(ixG^LL,ixM^LL,wprim,ps(igrid)%x)
      temperature(ixM^T)=wprim(ixM^T,p_)/wprim(ixM^T,rho_)
      phase(ixM^T)=kvec1*ps(igrid)%x(ixM^T,1)+kvec2*ps(igrid)%x(ixM^T,2)
      sinphase(ixM^T)=dsin(phase(ixM^T))
      qpar_k(ixM^T)=zero
      qperp_k(ixM^T)=zero
      if(mhd_hyperbolic_tc) then
        qpar_k(ixM^T)=wprim(ixM^T,qpar_)*costheta
        qperp_k(ixM^T)=-wprim(ixM^T,qperp_)* &
             sign(one,sinphase(ixM^T))*sintheta
      end if

      local_sum(1)=local_sum(1)+sum(dvolume(ixM^T))
      local_sum(2)=local_sum(2)+sum((temperature(ixM^T)/t0-one)* &
           dcos(phase(ixM^T))*dvolume(ixM^T))
      local_sum(3)=local_sum(3)+sum(qpar_k(ixM^T)*sinphase(ixM^T)*dvolume(ixM^T))
      local_sum(4)=local_sum(4)+sum(qperp_k(ixM^T)*sinphase(ixM^T)*dvolume(ixM^T))
      local_sum(5)=local_sum(5)+sum(temperature(ixM^T)*dvolume(ixM^T))
      local_sum(6)=local_sum(6)+sum(ps(igrid)%w(ixM^T,e_)*dvolume(ixM^T))
      local_sum(7)=local_sum(7)+sum(sum(wprim(ixM^T,mom(:))**2,dim=ndim+1)*dvolume(ixM^T))
      local_sum(8)=local_sum(8)+sum(wprim(ixM^T,rho_)*dvolume(ixM^T))
    end do
    call MPI_ALLREDUCE(local_sum,global_sum,8,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)

    amplitude_T=two*global_sum(2)/global_sum(1)
    amplitude_qpar=two*global_sum(3)/global_sum(1)
    amplitude_qperp=two*global_sum(4)/global_sum(1)
    amplitude_qtotal=amplitude_qpar+amplitude_qperp
    mean_temperature=global_sum(5)/global_sum(1)
    mean_energy=global_sum(6)/global_sum(1)
    vrms=dsqrt(global_sum(7)/global_sum(1))

    if(mype==0) then
      filename=trim(base_filename)//'.log'
      if(first_call .and. restart_from_file==undefined) then
        open(log_unit,file=trim(filename),status='replace')
        write(log_unit,'(a)') &
          '# time it dt A_T Aq_par_k Aq_perp_k Aq_total_k ' // &
          'mean_T mean_energy vrms tau_target gamma_pred'
        close(log_unit)
      end if
      open(log_unit,file=trim(filename),status='old',position='append')
      write(log_unit,'(es18.10,1x,i10,1x,10(es18.10,1x))') global_time,it,dt, &
           amplitude_T,amplitude_qpar,amplitude_qperp,amplitude_qtotal, &
           mean_temperature,mean_energy,vrms,tau_target,gamma_pred
      close(log_unit)
      first_call=.false.
    end if
  end subroutine print_fourier_log

end module mod_usr
