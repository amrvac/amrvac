{#IFDEF PARTICLES
!=============================================================================
subroutine init_tracerparticles

use mod_particles
!-----------------------------------------------------------------------------

call init_particles_vars
call init_particles_com
call init_particle_integrator

end subroutine init_tracerparticles
!=============================================================================
subroutine handle_particles()

use mod_gridvars
use mod_particles
use mod_timing

include 'mpif.h'
!-----------------------------------------------------------------------------

call set_neighbor_ipe

tpartc_grid_0=MPI_WTIME()
call init_gridvars
tpartc_grid = tpartc_grid + (MPI_WTIME()-tpartc_grid_0)

tpartc_com0=MPI_WTIME()
call comm_particles_global
tpartc_com=tpartc_com + (MPI_WTIME()-tpartc_com0)


! main integration loop:

particle_evol: do


   call select_active_particles


   call set_particles_dt


   tpartc_io_0=MPI_WTIME()
   call check_particles_output
   timeio_tot=timeio_tot+(MPI_WTIME()-tpartc_io_0)
   tpartc_io=tpartc_io+(MPI_WTIME()-tpartc_io_0)


   if (exit_condition() .eqv. .true.) exit particle_evol 


   tpartc_int_0=MPI_WTIME()
   call integrate_particles
   tpartc_int=tpartc_int+(MPI_WTIME()-tpartc_int_0)


   tpartc_com0=MPI_WTIME()
   call comm_particles
   tpartc_com=tpartc_com + (MPI_WTIME()-tpartc_com0)


   it_particles = it_particles + 1


end do particle_evol



tpartc_com0=MPI_WTIME()
call comm_particles
tpartc_com=tpartc_com + (MPI_WTIME()-tpartc_com0)

call finish_gridvars

end subroutine handle_particles
!=============================================================================
subroutine finish_tracerparticles

use mod_particles
!-----------------------------------------------------------------------------

call finish_particles
call finish_particles_com
call finish_particles_vars

end subroutine finish_tracerparticles
!=============================================================================
subroutine time_spent_on_particles()

use mod_timing
include 'amrvacdef.f'

double precision         :: tpartc_avg, tpartc_int_avg, tpartc_io_avg, tpartc_com_avg, tpartc_grid_avg

!-----------------------------------------------------------------------------
call MPI_REDUCE(tpartc,tpartc_avg,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
call MPI_REDUCE(tpartc_int,tpartc_int_avg,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
call MPI_REDUCE(tpartc_io,tpartc_io_avg,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
call MPI_REDUCE(tpartc_com,tpartc_com_avg,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
call MPI_REDUCE(tpartc_grid,tpartc_grid_avg,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)

if( mype ==0 ) then
   tpartc_avg     = tpartc_avg/dble(npe)
   tpartc_int_avg = tpartc_int_avg/dble(npe)
   tpartc_io_avg  = tpartc_io_avg/dble(npe)
   tpartc_com_avg  = tpartc_com_avg/dble(npe)
   tpartc_grid_avg  = tpartc_grid_avg/dble(npe)
   write(*,'(a,f12.3,a)')' Particle handling took     : ',tpartc,' sec'
   write(*,'(a,f12.2,a)')'                  Percentage: ',100.0d0*tpartc/timeloop,' %'
   write(*,'(a,f12.3,a)')' Particle IO took           : ',tpartc_io_avg,' sec'
   write(*,'(a,f12.2,a)')'                  Percentage: ',100.0d0*tpartc_io_avg/timeloop,' %'
   write(*,'(a,f12.3,a)')' Particle COM took          : ',tpartc_com_avg,' sec'
   write(*,'(a,f12.2,a)')'                  Percentage: ',100.0d0*tpartc_com_avg/timeloop,' %'
   write(*,'(a,f12.3,a)')' Particle integration took  : ',tpartc_int_avg,' sec'
   write(*,'(a,f12.2,a)')'                  Percentage: ',100.0d0*tpartc_int_avg/timeloop,' %'
   write(*,'(a,f12.3,a)')' Particle init grid took    : ',tpartc_grid_avg,' sec'
   write(*,'(a,f12.2,a)')'                  Percentage: ',100.0d0*tpartc_grid_avg/timeloop,' %'
end if


end subroutine time_spent_on_particles
!=============================================================================
subroutine read_particles_snapshot()

use mod_particles
include 'amrvacdef.f'
logical,save                    :: file_exists=.false.
character(len=128)              :: filename
integer                         :: mynpayload, mynparticles
integer, dimension(0:1)         :: buff
!-----------------------------------------------------------------------------
! some initialisations:
nparticles_on_mype = 0
mynparticles       = 0
nparticles         = 0
it_particles       = 0

! open the snapshot file on the headnode
if (mype .eq. 0) then
   write(filename,"(a,a,i4.4,a)") trim(filenameini),'_particles',snapshotini,'.dat'
   INQUIRE(FILE=filename, EXIST=file_exists)

   if (.not. file_exists) then
      write(*,*) 'File '//trim(filename)//' with particle data does not exist, calling init_particles instead'
      buff(0) = -1
      buff(1) = -1

   else
      open(unit=unitparticles,file=filename,form='unformatted',action='read',status='unknown',access='stream')
      read(unitparticles) nparticles,it_particles,mynpayload
      if (mynpayload .ne. npayload) &
           call mpistop('npayload in restart file does not match npayload in mod_particles')
      buff(0) = nparticles
      buff(1) = it_particles
   end if

end if
 
if (npe>0) call MPI_BCAST(buff,2,MPI_INTEGER,0,icomm,ierrmpi)

! check if the particle data was found:
if (buff(1) .eq. -1) then 

   call init_particles
   return

end if

! particle data is there, fill variables:
nparticles   = buff(0)
it_particles = buff(1)


do while (mynparticles .lt. nparticles)

   if (mype .eq. 0) then

      do while (nparticles_on_mype .lt. nparticles_per_cpu_hi &
           .and. mynparticles .lt. nparticles)
         call read_from_snapshot
         mynparticles = mynparticles + 1
      end do

   end if ! mype==0

   if (npe>0) call MPI_BCAST(mynparticles,1,MPI_INTEGER,0,icomm,ierrmpi)
      
   call comm_particles_global
end do


if (mype .eq. 0) close(unit=unitparticles)

itsavelast_particles = it_particles
end subroutine read_particles_snapshot
!=============================================================================
subroutine write_particles_snapshot()

use mod_particles
include 'amrvacdef.f'
character(len=128)              :: filename
type(particle_node), dimension(0:nparticles_per_cpu_hi-1)  :: send_particles
type(particle_node), dimension(0:nparticles_per_cpu_hi-1)  :: receive_particles
integer                         :: status(MPI_STATUS_SIZE)
integer,dimension(0:npe-1)      :: receive_n_particles_for_output_from_ipe
integer                         :: ipe, ipart, iipart, send_n_particles_for_output
logical,save                    :: file_exists=.false.
!-----------------------------------------------------------------------------
receive_n_particles_for_output_from_ipe(:) = 0

! open the snapshot file on the headnode
if (mype .eq. 0) then
   write(filename,"(a,a,i4.4,a)") trim(filenameout),'_particles',snapshot-1,'.dat'
   INQUIRE(FILE=filename, EXIST=file_exists)
   if (.not. file_exists) then
      open(unit=unitparticles,file=filename,form='unformatted',status='new',access='stream')
   else
      open(unit=unitparticles,file=filename,form='unformatted',status='replace',access='stream')
   end if
   write(unitparticles) nparticles,it_particles,npayload
end if


if (npe==1) then 
   do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);
      call append_to_snapshot(particle(ipart)%self)
   end do
   return
end if


if (mype .ne. 0) then
   call MPI_SEND(nparticles_on_mype,1,MPI_INTEGER,0,mype,icomm,ierrmpi)
! fill the send_buffer
   send_n_particles_for_output = nparticles_on_mype
   do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);
      send_particles(iipart-1) = particle(ipart)%self
   end do
else
! directly fill the receive buffer
   receive_n_particles_for_output_from_ipe(0) = nparticles_on_mype
   do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);
      receive_particles(iipart-1) = particle(ipart)%self
   end do
! get number of particles on other ipes
   do ipe=1,npe-1
      call MPI_RECV(receive_n_particles_for_output_from_ipe(ipe),1,MPI_INTEGER,ipe,ipe,icomm,status,ierrmpi)
   end do
end if



if (mype .ne. 0) &
     call MPI_SEND(send_particles,send_n_particles_for_output,type_particle,0,mype,icomm,ierrmpi)

if (mype==0) then

   ! first output own particles (already in receive buffer)
   do ipart=1,receive_n_particles_for_output_from_ipe(0)
      call append_to_snapshot(receive_particles(ipart-1))
   end do

   ! now output the particles sent from the other ipes
   do ipe=1,npe-1
      call MPI_RECV(receive_particles,receive_n_particles_for_output_from_ipe(ipe),type_particle,ipe,ipe,icomm,status,ierrmpi)
      do ipart=1,receive_n_particles_for_output_from_ipe(ipe)
         call append_to_snapshot(receive_particles(ipart-1))
      end do ! ipart
   end do ! ipe


   close(unit=unitparticles)
end if ! mype == 0

end subroutine write_particles_snapshot
!=============================================================================

}
