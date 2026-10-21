module mod_usr
  use mod_mhd
  use mod_eos, only: eos
  implicit none

contains

  subroutine usr_init()
    usr_init_one_grid => initonegrid_usr
    usr_print_log => uawsom_log
    call set_coordinate_system('Cartesian_1D')
    call mhd_activate()
  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: pulse(ixI^S)

    w(ixO^S,rho_)=one
    w(ixO^S,mom(:))=zero
    w(ixO^S,p_)=one
    w(ixO^S,mag(:))=zero
    w(ixO^S,mag(1))=one
    w(ixO^S,wAplus_)=1.d-8
    w(ixO^S,wAminus_)=1.d-8
    w(ixO^S,wkplus_)=1.d-8
    w(ixO^S,wkminus_)=1.d-8

    select case(iprob)
    case(1)
      pulse(ixO^S)=dexp(-((x(ixO^S,1)-0.35d0)/0.06d0)**2)
      w(ixO^S,wAplus_)=1.d-3*pulse(ixO^S)+1.d-8
      w(ixO^S,wAminus_)=5.d-4*dexp(-((x(ixO^S,1)-0.65d0)/0.06d0)**2)+1.d-8
      w(ixO^S,wkplus_)=7.d-4*pulse(ixO^S)+1.d-8
      w(ixO^S,wkminus_)=3.d-4*dexp(-((x(ixO^S,1)-0.65d0)/0.06d0)**2)+1.d-8
    case(2)
      ! A density gradient gives dv_A/dx while total wave energy starts uniform.
      w(ixO^S,rho_)=one+0.2d0*dsin(two*dpi*x(ixO^S,1))
      w(ixO^S,wAplus_)=2.d-3
      w(ixO^S,wAminus_)=1.d-3
    case(3)
      w(ixO^S,wAplus_)=2.d-3
      w(ixO^S,wAminus_)=1.d-3
      w(ixO^S,wkplus_)=2.d-3
      w(ixO^S,wkminus_)=1.d-3
    case(4)
      ! The ordinary small-value replacement path must repair this value.
      w(ixO^S,wAplus_)=-1.d-6
      w(ixO^S,wAminus_)=1.d-3
    case default
      call mpistop('UAWSoM_1D: iprob must be 1, 2, 3, or 4')
    end select

    call eos%to_conserved(ixI^L,ixO^L,w,x)
  end subroutine initonegrid_usr

  subroutine uawsom_log()
    use mod_global_parameters
    use mod_input_output, only: get_global_minima
    integer, parameter :: log_unit=731
    integer :: iigrid,igrid,ix1
    double precision :: sums(10),global_sums(10),wmin(nw),volume
    character(len=std_len) :: filename
    logical, save :: first=.true.

    sums=zero
    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      do ix1=ixMlo1,ixMhi1
        volume=ps(igrid)%dvolume(ix1)
        sums(1)=sums(1)+volume*ps(igrid)%w(ix1,wAplus_)
        sums(2)=sums(2)+volume*ps(igrid)%w(ix1,wAminus_)
        sums(3)=sums(3)+volume*ps(igrid)%w(ix1,wkplus_)
        sums(4)=sums(4)+volume*ps(igrid)%w(ix1,wkminus_)
        sums(5)=sums(5)+volume*ps(igrid)%w(ix1,e_)
        sums(6)=sums(6)+volume*ps(igrid)%x(ix1,1)*ps(igrid)%w(ix1,wAplus_)
        sums(7)=sums(7)+volume*ps(igrid)%x(ix1,1)*ps(igrid)%w(ix1,wAminus_)
        sums(8)=sums(8)+volume*(ps(igrid)%w(ix1,wAplus_)+&
             ps(igrid)%w(ix1,wAminus_))
        sums(9)=sums(9)+volume*mhd_uawsom_wave_pressure_cell(&
             ps(igrid)%w(ix1,:),5.d0)
        sums(10)=sums(10)+volume*mhd_uawsom_rho2_factor_cell(5.d0)
      end do
    end do
    call MPI_ALLREDUCE(sums,global_sums,10,MPI_DOUBLE_PRECISION,&
         MPI_SUM,icomm,ierrmpi)
    call get_global_minima(wmin,ps)

    if(mype==0) then
      filename=trim(base_filename)//'_uawsom.log'
      if(first .and. restart_from_file==undefined) then
        open(log_unit,file=trim(filename),status='replace')
        write(log_unit,'(a)') '# time WA+ WA- Wk+ Wk- Etot xWA+ xWA- WA_sum Pwave rho2_factor minW'
      else
        open(log_unit,file=trim(filename),status='old',position='append')
      end if
      write(log_unit,'(12(es16.8,1x))') global_time,global_sums(1:10),&
           min(wmin(wAplus_),wmin(wAminus_),wmin(wkplus_),wmin(wkminus_))
      close(log_unit)
    end if
    first=.false.
  end subroutine uawsom_log

end module mod_usr
