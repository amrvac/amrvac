!=============================================================================
subroutine write_analysis
  ! This is an example file how to use the analysis capability.  You can schedule 
  ! this routine using the slot 5 in itsave, dtsave and ditsave.  
  ! To use, just copy this file to your working directory and make your modifications.
  use constants
  use mod_forest, only: Morton_sub_start, Morton_sub_stop
  include 'amrvacdef.f'
  character(len=20):: userconvert_type
  !-----------------------------------------------------------------------------
  logical :: fileopen
  integer :: iigrid, igrid, level
  double precision :: volume(1:nlevelshi)
  double precision :: e2, b2, jmax, jmin, jbar, voltotal, val, vbl, dal, dbl, valph_a, valph_b, va_a, va_b
  double precision :: dvolume(ixG^T), current(ixG^T,1:ndir), jabs(ixG^T)
  double precision, dimension(1:8) :: dsum_send, dsum_recv
  double precision, dimension(1:2) :: dmax_send, dmax_recv
  character(len=80) :: filename
  logical, save :: opened=.false.
  logical, save :: file_exists=.false.
  integer :: amode, status(MPI_STATUS_SIZE), Njgrid, jgrid, ix2
  double precision,dimension(0:nw+nwauxio)       :: normconv 
  double precision, parameter                    :: xslicea=0.05, xsliceb=0.1, slicewidth=0.1
  logical                                        :: sliceascii_save, saveprim_save
  !-----------------------------------------------------------------------------


  volume(1:mxnest)=zero

  e2   = zero
  b2   = zero
  jmax = -bigdouble
  jmin = +bigdouble
  jbar = zero
  val  = zero
  vbl  = zero
  dal  = zero
  dbl  = zero
  valph_a = zero
  valph_b = zero
  va_a = zero
  va_b = zero


  do iigrid=1,igridstail; igrid=igrids(iigrid);
     level = node(plevel_,igrid)

     if (slab) then
        dvolume(ixM^T)={rnode(rpdx^D_,igrid)|*}
     else
        dvolume(ixM^T)=pgeo(igrid)%dvolume(ixM^T)
     end if
     volume(level)=volume(level)+sum(dvolume(ixM^T))

     ! Compute integral quantities:
     e2 = e2 + sum( ({^C& pw(igrid)%w(ixM^T,e^C_)**2|+})*dvolume(ixM^T) )
     b2 = b2 + sum( ({^C& pw(igrid)%w(ixM^T,b^C_)**2|+})*dvolume(ixM^T) )
     call getcurrent(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,.false.,current)
     jabs(ixM^T) = sqrt({^C& current(ixM^T,^C)**2|+})
     jbar = jbar + sum( jabs(ixM^T)*dvolume(ixM^T) )

     ! Extrema:
     jmax = max(jmax, maxval(jabs(ixM^T)))
     jmin = min(jmin, minval(jabs(ixM^T)))
  end do

  voltotal=sum(volume(levmin:levmax))

  dsum_send(1)=e2
  dsum_send(2)=b2
  dsum_send(3)=jbar
  dsum_send(4)=voltotal

  call MPI_REDUCE(dsum_send,dsum_recv,4,MPI_DOUBLE_PRECISION, &
       MPI_SUM,0,icomm,ierrmpi)

  dmax_send(1) = -jmin
  dmax_send(2) =  jmax
  call MPI_REDUCE(dmax_send,dmax_recv,2,MPI_DOUBLE_PRECISION, &
       MPI_MAX,0,icomm,ierrmpi)

  if (mype==0) then
     e2 = dsum_recv(1)
     b2 = dsum_recv(2)
     jbar = dsum_recv(3)
     voltotal = dsum_recv(4)
     e2 = e2/voltotal
     b2 = b2/voltotal
     jbar = jbar/voltotal
     jmin = - dmax_recv(1)
     jmax = dmax_recv(2)
  end if


! Compute reconnection rate through two slices:
! requires sliceascii and saveprim both true!

{#IFNDEF D2
call mpistop('write_analysis: Only implemented for 2D problem!')
}
  sliceascii_save = sliceascii
  sliceascii      = .true.
  saveprim_save   = saveprim
  sliceascii      = .true.
  
  !--------------------------------------
  ! Inner slice:
  !--------------------------------------

  call select_slice(1,-xslicea,.false.,unitslice,normconv)
  ! local number of sub-grids:
  Njgrid = Morton_sub_stop(mype) - Morton_sub_start(mype) + 1
  if (Njgrid>0) then 
     do jgrid=1,Njgrid
        do ix2=ixMlo2,ixMhi2
           if (abs(px_sub(jgrid)%x(ix2,2)) .le. slicewidth) then 
              val = val + abs(pw_sub(jgrid)%w(ix2,u1_)/pw_sub(jgrid)%w(ix2,lfac_)) * rnode_sub(rpdx1_,jgrid)
              valph_a = valph_a &
                   + sqrt(({^C& pw_sub(jgrid)%w(ix2,b^C_)**2|+})/(({^C& pw_sub(jgrid)%w(ix2,b^C_)**2|+})+pw_sub(jgrid)%w(ix2,xi_)/pw_sub(jgrid)%w(ix2,lfac_)**2)) & 
                   * abs(pw_sub(jgrid)%w(ix2,b1_)/sqrt(pw_sub(jgrid)%w(ix2,b1_)**2+pw_sub(jgrid)%w(ix2,b3_)**2)) &
                   * rnode_sub(rpdx1_,jgrid)
              va_a = va_a &
                   + sqrt(({^C& pw_sub(jgrid)%w(ix2,b^C_)**2|+})/(({^C& pw_sub(jgrid)%w(ix2,b^C_)**2|+})+pw_sub(jgrid)%w(ix2,xi_)/pw_sub(jgrid)%w(ix2,lfac_)**2)) & 
                   * rnode_sub(rpdx1_,jgrid)
              dal  = dal + rnode_sub(rpdx1_,jgrid)
           end if
        end do
     end do
     do jgrid=1,Njgrid
        call dealloc_subnode(jgrid)
     end do
  end if ! Njgrid>0

  call select_slice(1,xslicea,.false.,unitslice,normconv)
  ! local number of sub-grids:
  Njgrid = Morton_sub_stop(mype) - Morton_sub_start(mype) + 1
  if (Njgrid>0) then 
     do jgrid=1,Njgrid
        do ix2=ixMlo2,ixMhi2
           if (abs(px_sub(jgrid)%x(ix2,2)) .le. slicewidth) then 
              val = val + abs(pw_sub(jgrid)%w(ix2,u1_)/pw_sub(jgrid)%w(ix2,lfac_)) * rnode_sub(rpdx1_,jgrid)
              valph_a = valph_a &
                   + sqrt(({^C& pw_sub(jgrid)%w(ix2,b^C_)**2|+})/(({^C& pw_sub(jgrid)%w(ix2,b^C_)**2|+})+pw_sub(jgrid)%w(ix2,xi_)/pw_sub(jgrid)%w(ix2,lfac_)**2)) & 
                   * abs(pw_sub(jgrid)%w(ix2,b1_)/sqrt(pw_sub(jgrid)%w(ix2,b1_)**2+pw_sub(jgrid)%w(ix2,b3_)**2)) &
                   * rnode_sub(rpdx1_,jgrid)
              va_a = va_a &
                   + sqrt(({^C& pw_sub(jgrid)%w(ix2,b^C_)**2|+})/(({^C& pw_sub(jgrid)%w(ix2,b^C_)**2|+})+pw_sub(jgrid)%w(ix2,xi_)/pw_sub(jgrid)%w(ix2,lfac_)**2)) & 
                   * rnode_sub(rpdx1_,jgrid)
              dal  = dal + rnode_sub(rpdx1_,jgrid)
           end if
        end do
     end do
     do jgrid=1,Njgrid
        call dealloc_subnode(jgrid)
     end do
  end if ! Njgrid>0


  !--------------------------------------
  ! Outer slice:
  !--------------------------------------
  call select_slice(1,-xsliceb,.false.,unitslice,normconv)
  ! local number of sub-grids:
  Njgrid = Morton_sub_stop(mype) - Morton_sub_start(mype) + 1
  if (Njgrid>0) then 
     do jgrid=1,Njgrid
        do ix2=ixMlo2,ixMhi2
           if (abs(px_sub(jgrid)%x(ix2,2)) .le. slicewidth) then 
              vbl = vbl + abs(pw_sub(jgrid)%w(ix2,u1_)/pw_sub(jgrid)%w(ix2,lfac_)) * rnode_sub(rpdx1_,jgrid)
              valph_b = valph_b &
                   + sqrt(({^C& pw_sub(jgrid)%w(ix2,b^C_)**2|+})/(({^C& pw_sub(jgrid)%w(ix2,b^C_)**2|+})+pw_sub(jgrid)%w(ix2,xi_)/pw_sub(jgrid)%w(ix2,lfac_)**2)) & 
                   * abs(pw_sub(jgrid)%w(ix2,b1_)/sqrt(pw_sub(jgrid)%w(ix2,b1_)**2+pw_sub(jgrid)%w(ix2,b3_)**2)) &
                   * rnode_sub(rpdx1_,jgrid)
              va_b = va_b &
                   + sqrt(({^C& pw_sub(jgrid)%w(ix2,b^C_)**2|+})/(({^C& pw_sub(jgrid)%w(ix2,b^C_)**2|+})+pw_sub(jgrid)%w(ix2,xi_)/pw_sub(jgrid)%w(ix2,lfac_)**2)) & 
                   * rnode_sub(rpdx1_,jgrid)
              dbl  = dbl + rnode_sub(rpdx1_,jgrid)
           end if
        end do
     end do
     do jgrid=1,Njgrid
        call dealloc_subnode(jgrid)
     end do
  end if ! Njgrid>0

  call select_slice(1,xsliceb,.false.,unitslice,normconv)
  ! local number of sub-grids:
  Njgrid = Morton_sub_stop(mype) - Morton_sub_start(mype) + 1
  if (Njgrid>0) then 
     do jgrid=1,Njgrid
        do ix2=ixMlo2,ixMhi2
           if (abs(px_sub(jgrid)%x(ix2,2)) .le. slicewidth) then 
              vbl = vbl + abs(pw_sub(jgrid)%w(ix2,u1_)/pw_sub(jgrid)%w(ix2,lfac_)) * rnode_sub(rpdx1_,jgrid)
              valph_b = valph_b &
                   + sqrt(({^C& pw_sub(jgrid)%w(ix2,b^C_)**2|+})/(({^C& pw_sub(jgrid)%w(ix2,b^C_)**2|+})+pw_sub(jgrid)%w(ix2,xi_)/pw_sub(jgrid)%w(ix2,lfac_)**2)) & 
                   * abs(pw_sub(jgrid)%w(ix2,b1_)/sqrt(pw_sub(jgrid)%w(ix2,b1_)**2+pw_sub(jgrid)%w(ix2,b3_)**2)) &
                   * rnode_sub(rpdx1_,jgrid)
              va_b = va_b &
                   + sqrt(({^C& pw_sub(jgrid)%w(ix2,b^C_)**2|+})/(({^C& pw_sub(jgrid)%w(ix2,b^C_)**2|+})+pw_sub(jgrid)%w(ix2,xi_)/pw_sub(jgrid)%w(ix2,lfac_)**2)) & 
                   * rnode_sub(rpdx1_,jgrid)
              dbl  = dbl + rnode_sub(rpdx1_,jgrid)
           end if
        end do
     end do
     do jgrid=1,Njgrid
        call dealloc_subnode(jgrid)
     end do
  end if ! Njgrid>0



  dsum_send(1) = val
  dsum_send(2) = vbl
  dsum_send(3) = dal
  dsum_send(4) = dbl
  dsum_send(5) = valph_a
  dsum_send(6) = valph_b
  dsum_send(7) = va_a
  dsum_send(8) = va_b

  call MPI_REDUCE(dsum_send,dsum_recv,8,MPI_DOUBLE_PRECISION, &
       MPI_SUM,0,icomm,ierrmpi)


  if (mype .eq. 0) then

     val = dsum_recv(1) / dsum_recv(3)
     vbl = dsum_recv(2) / dsum_recv(4)
     valph_a = dsum_recv(5) / dsum_recv(3)
     valph_b = dsum_recv(6) / dsum_recv(4)
     va_a = dsum_recv(7) / dsum_recv(3)
     va_b = dsum_recv(8) / dsum_recv(4)
     
     if (.not.opened) then
        ! generate filename
        write(filename,"(a,a,a)") trim(filenameout),TRIM("_analysis"),".csv"
        INQUIRE(FILE=filename, EXIST=file_exists)
        if (.not. file_exists) then
           open(unit=unitanalysis,file=filename,status='unknown',access='append')
           write(unitanalysis,"(a)") trim('# t e2 b2 jbar jmin jmax vra vrb vaacostheta vabcostheta vaa vab')
        else
           open(unit=unitanalysis,file=filename,status='unknown',access='append')
        end if
        opened=.true.
     end if
     write(unitanalysis,'(12(es14.6))')t, e2, b2, jbar, jmin, jmax, val, vbl, valph_a, valph_b, va_a, va_b
     
  end if


  sliceascii = sliceascii_save
  saveprim   = saveprim_save

end subroutine write_analysis
!=============================================================================
