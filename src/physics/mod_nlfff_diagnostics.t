!> Common user-facing diagnostics for NLFFF relaxation and extrapolation.
module mod_nlfff_diagnostics
  implicit none
  private

  type, public :: nlfff_physical_metrics
    double precision :: cw_sin_theta=0.d0
    double precision :: epsilon_force=0.d0
    double precision :: epsilon_div=0.d0
    double precision :: magnetic_energy=0.d0
  end type nlfff_physical_metrics

{^IFTHREED
  public :: evaluate_nlfff_metrics_amrvac
  public :: evaluate_nlfff_metrics_dense
  public :: write_nlfff_metrics_header
  public :: write_nlfff_metrics_row
}

contains

{^IFTHREED
  !> Evaluate the common metrics on all active AMRVAC cells.  The local cell
  !> size makes epsilon_force and epsilon_div dimensionless on AMR meshes.
  subroutine evaluate_nlfff_metrics_amrvac(iw_b,metrics)
    use mpi
    use mod_comm_lib, only: mpistop
    use mod_global_parameters
    use mod_geometry, only: curlvector,divvector
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    integer, intent(in) :: iw_b(3)
    type(nlfff_physical_metrics), intent(out) :: metrics

    double precision :: bvec(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,3)
    double precision :: current(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,3)
    double precision :: divb(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
    double precision :: local(5),global(5),b(3),j(3),jxb(3)
    double precision :: b2,j2,jxb2,volume,hcell
    integer :: iigrid,igrid,ix1,ix2,ix3,idirmin

    local=0.d0
    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      block=>ps(igrid)
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      bvec=0.d0
      bvec(:,:,:,1:3)=ps(igrid)%w(:,:,:,iw_b(:))
      if(B0field) bvec(:,:,:,1:3)=bvec(:,:,:,1:3)+&
         block%B0(:,:,:,iw_b(:),0)
      current=0.d0
      divb=0.d0
      idirmin=1
      call curlvector(bvec,ixG^LL,ixM^LL,current,idirmin,1,3)
      call divvector(bvec,ixG^LL,ixM^LL,divb,1)
      do ix3=ixMlo3,ixMhi3
        do ix2=ixMlo2,ixMhi2
          do ix1=ixMlo1,ixMhi1
            b=bvec(ix1,ix2,ix3,:)
            j=current(ix1,ix2,ix3,:)
            jxb=nlfff_cross3(j,b)
            b2=dot_product(b,b)
            j2=dot_product(j,j)
            jxb2=dot_product(jxb,jxb)
            volume=block%dvolume(ix1,ix2,ix3)
            hcell=volume**(1.d0/3.d0)
            if(b2>0.d0) then
              local(1)=local(1)+dsqrt(jxb2/b2)*volume
              local(3)=local(3)+hcell**2*jxb2/b2*volume
            end if
            local(2)=local(2)+dsqrt(j2)*volume
            local(4)=local(4)+hcell**2*divb(ix1,ix2,ix3)**2*volume
            local(5)=local(5)+b2*volume
          end do
        end do
      end do
    end do
    call MPI_ALLREDUCE(local,global,5,MPI_DOUBLE_PRECISION,MPI_SUM,&
       icomm,ierrmpi)
    if(.not.all(ieee_is_finite(global))) &
       call mpistop('non-finite common NLFFF diagnostic')
    metrics=nlfff_physical_metrics()
    if(global(2)>0.d0) metrics%cw_sin_theta=global(1)/global(2)
    if(global(5)>0.d0) then
      metrics%epsilon_force=dsqrt(max(0.d0,global(3))/global(5))
      metrics%epsilon_div=dsqrt(max(0.d0,global(4))/global(5))
    end if
    metrics%magnetic_energy=0.5d0*global(5)
  end subroutine evaluate_nlfff_metrics_amrvac

  !> Evaluate the same metrics on a replicated uniform dense field. Plane 1
  !> is the lower boundary and planes 2:nz are the active-cell centres.
  subroutine evaluate_nlfff_metrics_dense(b,dx1,dx2,dx3,metrics)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    double precision, intent(in) :: b(:,:,:,:),dx1,dx2,dx3
    type(nlfff_physical_metrics), intent(out) :: metrics
    double precision :: db(3,3),j(3),bv(3),jxb(3),divb,b2,j2,jxb2
    double precision :: sum_cross,sum_j,sum_force,sum_div,sum_b2,h,volume
    integer :: i,jj,k,ic

    sum_cross=0.d0; sum_j=0.d0; sum_force=0.d0
    sum_div=0.d0; sum_b2=0.d0
    h=(dx1*dx2*dx3)**(1.d0/3.d0)
    volume=dx1*dx2*dx3
    do k=2,size(b,3)
      do jj=1,size(b,2)
        do i=1,size(b,1)
          do ic=1,3
            db(ic,1)=dense_derivative(b(:,jj,k,ic),i,dx1)
            db(ic,2)=dense_derivative(b(i,:,k,ic),jj,dx2)
            db(ic,3)=dense_derivative(b(i,jj,:,ic),k,dx3)
          end do
          j=(/db(3,2)-db(2,3),db(1,3)-db(3,1),db(2,1)-db(1,2)/)
          bv=b(i,jj,k,:)
          jxb=nlfff_cross3(j,bv)
          b2=dot_product(bv,bv)
          j2=dot_product(j,j)
          jxb2=dot_product(jxb,jxb)
          if(b2>0.d0) then
            sum_cross=sum_cross+dsqrt(jxb2/b2)*volume
            sum_force=sum_force+h**2*jxb2/b2*volume
          end if
          sum_j=sum_j+dsqrt(j2)*volume
          divb=db(1,1)+db(2,2)+db(3,3)
          sum_div=sum_div+h**2*divb**2*volume
          sum_b2=sum_b2+b2*volume
        end do
      end do
    end do
    if(.not.ieee_is_finite(sum_cross+sum_j+sum_force+sum_div+sum_b2)) &
       error stop 'non-finite common dense NLFFF diagnostic'
    metrics=nlfff_physical_metrics()
    if(sum_j>0.d0) metrics%cw_sin_theta=sum_cross/sum_j
    if(sum_b2>0.d0) then
      metrics%epsilon_force=dsqrt(max(0.d0,sum_force)/sum_b2)
      metrics%epsilon_div=dsqrt(max(0.d0,sum_div)/sum_b2)
    end if
    metrics%magnetic_energy=0.5d0*sum_b2
  end subroutine evaluate_nlfff_metrics_dense

  pure function nlfff_cross3(a,b) result(c)
    double precision, intent(in) :: a(3),b(3)
    double precision :: c(3)
    c=(/a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),&
       a(1)*b(2)-a(2)*b(1)/)
  end function nlfff_cross3

  pure double precision function dense_derivative(f,i,h)
    double precision, intent(in) :: f(:),h
    integer, intent(in) :: i
    integer :: n

    n=size(f)
    if(n<3) then
      dense_derivative=0.d0
    else if(i==1) then
      dense_derivative=(-3.d0*f(1)+4.d0*f(2)-f(3))/(2.d0*h)
    else if(i==n) then
      dense_derivative=(3.d0*f(n)-4.d0*f(n-1)+f(n-2))/(2.d0*h)
    else
      dense_derivative=(f(i+1)-f(i-1))/(2.d0*h)
    end if
  end function dense_derivative

  subroutine write_nlfff_metrics_header(unit)
    integer, intent(in) :: unit
    write(unit,'(a)') &
       'iteration,CW_sin_theta,epsilon_force,epsilon_div,magnetic_energy'
    flush(unit)
  end subroutine write_nlfff_metrics_header

  subroutine write_nlfff_metrics_row(unit,iteration,metrics)
    integer, intent(in) :: unit,iteration
    type(nlfff_physical_metrics), intent(in) :: metrics
    write(unit,'(i0,4(",",es24.16))') iteration,metrics%cw_sin_theta,&
       metrics%epsilon_force,metrics%epsilon_div,metrics%magnetic_energy
    flush(unit)
  end subroutine write_nlfff_metrics_row
}

end module mod_nlfff_diagnostics
