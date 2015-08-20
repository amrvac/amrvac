{#IFDEF PARTICLES
!=============================================================================
module mod_particles
!
! Contains global variables and routines for particles
!
implicit none
  
integer,parameter                      :: nparticleshi = 100000, nparticles_per_cpu_hi = 100000
integer, parameter                     :: npayload=1
double precision, save                 :: t_particles, dt_particles

double precision, save                 :: tmax_particles
integer, save                          :: itmax_particles
integer, save                          :: ditsave_particles
double precision, save                 :: dtsave_ensemble
double precision, save                 :: dtheta
logical,save                           :: losses

integer, save                          :: nparticles, it_particles
integer, save                          :: itsavelast_particles

! these two save the list of neighboring cpus:
integer, dimension(:), allocatable,save     :: ipe_neighbor
integer, save                               :: npe_neighbors

integer, save                               :: type_particle

integer, parameter                      :: unitparticles=15 

! list of particles on the processor:
integer, dimension(nparticles_per_cpu_hi), save    :: particles_on_mype, particles_active_on_mype
integer, save                                      :: nparticles_on_mype, nparticles_active_on_mype

type particle_ptr
   type(particle_node), pointer         :: self
   integer                              :: igrid, ipe
end type particle_ptr

type particle_node
   logical                               :: follow
   integer                               :: index
   double precision                      :: q, m
   double precision                      :: t, dt
   double precision, dimension(npayload) :: payload
   double precision, dimension(^NC)      :: x
   double precision, dimension(^NC)      :: u
end type particle_node


type(particle_ptr),save,dimension(:),allocatable               :: particle

!INTERFACES
interface 
   function user_destroy(myparticle)
     import particle_node
     logical                         :: user_destroy
     type(particle_node), intent(in) :: myparticle
   end function user_destroy
end interface
!=============================================================================
contains
!=============================================================================
logical function exit_condition()
! check if we should go out of the integration loop

!-----------------------------------------------------------------------------

exit_condition = (&
     t_particles >= tmax_particles &
     .or. nparticles == 0 &
     .or. it_particles >= itmax_particles &
     )

end function exit_condition
!=============================================================================
subroutine init_particles_output()
! overwrite the output files

include 'amrvacdef.f'
character(128)                    :: filename
integer                           :: iipart, ipart
!-----------------------------------------------------------------------------

do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);
   if (particle(ipart)%self%follow) then
      filename = make_particle_filename(particle(ipart)%self%t,particle(ipart)%self%index,'individual')
      open(unit=unitparticles,file=filename,status='replace')
      write(unitparticles,"(a)") trim('# t, dt, x^C, u^C, payload(1:npayload), ipe, iteration, index')
      close(unit=unitparticles)
   end if
end do

end subroutine init_particles_output
!=============================================================================
subroutine init_particles_vars()
! allocate some dynamic arrays

include 'amrvacdef.f'
!-----------------------------------------------------------------------------

allocate(particle(1:nparticleshi))
allocate(ipe_neighbor(0:npe-1))

! initialise 

dt_particles      = bigdouble
t_particles       = 0.0d0
tmax_particles    = bigdouble
itmax_particles   = biginteger
ditsave_particles = 8 
dtsave_ensemble   = bigdouble
dtheta            = 2.0d0*dpi / 60.0d0

losses            = .false.

nparticles = 0
it_particles = 0
itsavelast_particles = 0

nparticles_on_mype   = 0
particles_on_mype(:) = 0

nparticles_active_on_mype   = 0
particles_active_on_mype(:) = 0

end subroutine init_particles_vars
!=============================================================================
subroutine finish_particles_vars()
! de-allocate some dynamic arrays

include 'amrvacdef.f'
deallocate(particle)
deallocate(ipe_neighbor)


end subroutine finish_particles_vars
!=============================================================================
subroutine select_active_particles

include 'amrvacdef.f'
integer                                         :: ipart, iipart
logical                                         :: activate
!-----------------------------------------------------------------------------

nparticles_active_on_mype = 0
do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);

   if (time_advance) then
      activate = particle(ipart)%self%t .lt. tmax_particles
   else
      activate = particle(ipart)%self%t .le. tmax_particles
   end if
   
   if (activate) then 
      nparticles_active_on_mype = nparticles_active_on_mype + 1
      particles_active_on_mype(nparticles_active_on_mype) = ipart
   end if
   
end do


end subroutine select_active_particles
!=============================================================================
subroutine locate_particle(index,igrid_particle,ipe_particle)
! given the particles unique index, tell me on which cpu and igrid it is
! returns -1,-1 if particle was not found
include 'amrvacdef.f'
integer, intent(in)                            :: index
integer, intent(out)                           :: igrid_particle, ipe_particle
integer                                        :: iipart,ipart,ipe_has_particle,ipe
logical                                        :: has_particle(0:npe-1)
integer,dimension(0:1)                         :: buff
!-----------------------------------------------------------------------------

has_particle(:) = .false.
do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);
   if (particle(ipart)%self%index == index) then
      has_particle(mype) = .true.
      exit
   end if
end do

if (has_particle(mype)) then
   buff(0) = particle(ipart)%igrid
   buff(1) = mype
end if

if (npe>0) call MPI_ALLGATHER(MPI_IN_PLACE,1,MPI_LOGICAL,has_particle,1,MPI_LOGICAL,icomm,ierrmpi)

ipe_has_particle = -1
do ipe=0,npe-1
   if (has_particle(ipe) .eqv. .true.) then 
      ipe_has_particle = ipe
      exit
   end if
end do

if (ipe_has_particle .ne. -1) then
   if (npe>0) call MPI_BCAST(buff,2,MPI_INTEGER,ipe_has_particle,icomm,ierrmpi)
   igrid_particle=buff(0)
   ipe_particle = buff(1)
else
   igrid_particle=-1
   ipe_particle = -1
end if

end subroutine locate_particle
!=============================================================================
subroutine find_particle_ipe(x,igrid_particle,ipe_particle)

use mod_forest, only: tree_node_ptr, tree_root
include 'amrvacdef.f'

double precision, dimension(ndir), intent(in)   :: x
integer, intent(out)                            :: igrid_particle, ipe_particle

integer, dimension(ndir,nlevelshi)              :: ig
integer                                         :: idim, ic(ndim)
type(tree_node_ptr)                             :: branch
!-----------------------------------------------------------------------------

! first check if the particle is in the domain
if (.not. particle_in_domain(x)) then
   igrid_particle = -1
   ipe_particle   = -1
   return
end if

! get the index on each level
do idim = 1, ndim
   call get_igslice(idim,x(idim),ig(idim,:))
end do

! traverse the tree until leaf is found
branch=tree_root({ig(^D,1)})
do while (.not.branch%node%leaf)
   {ic(^D)=ig(^D,branch%node%level+1) - 2 * branch%node%ig^D +2\}
   branch%node => branch%node%child({ic(^D)})%node
end do

igrid_particle = branch%node%igrid
ipe_particle   = branch%node%ipe

end subroutine find_particle_ipe
!=============================================================================
logical function particle_in_domain(x)

include 'amrvacdef.f'
double precision, dimension(ndim), intent(in)  :: x
integer                                        :: idim
!-----------------------------------------------------------------------------

particle_in_domain = .true.

do idim=1,ndim
   select case(idim)
      {case (^D)
      if (x(^D) .lt. xprobmin^D) then
         particle_in_domain = .false.
         exit
      end if
      if (x(^D) .ge. xprobmax^D) then
         particle_in_domain = .false.
         exit
      end if
      \}
   end select
end do

end function particle_in_domain
!=============================================================================
logical function particle_in_igrid(ipart,igrid)
! quick check if particle is still in igrid

include 'amrvacdef.f'
integer, intent(in)                            :: igrid,ipart
integer                                        :: idim
!-----------------------------------------------------------------------------

! first check if the igrid is still there:
if (.not. associated(pw(igrid)%w)) then
   particle_in_igrid = .false.
   return
end if

particle_in_igrid = .true.

do idim=1,ndim
   select case(idim)
      {case (^D)
      if (particle(ipart)%self%x(^D) .lt. rnode(rpxmin^D_,igrid)) then
         particle_in_igrid = .false.
         exit
      end if
      if (particle(ipart)%self%x(^D) .ge. rnode(rpxmax^D_,igrid)) then
         particle_in_igrid = .false.
         exit
      end if
      \}
   end select
end do

return
end function particle_in_igrid
!=============================================================================
subroutine init_particles_com()
! initialise communicators for particles
include 'amrvacdef.f'

integer(kind=MPI_ADDRESS_KIND)         :: size_int, size_double, size_logical, lb
integer                                :: oldtypes(0:8), offsets(0:8), &
     blocklengths(0:8)
!-----------------------------------------------------------------------------

! create the MPI-datatype for particle

call MPI_TYPE_GET_EXTENT(MPI_INTEGER, lb, size_int, ierrmpi)
call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION, lb, size_double, ierrmpi)
call MPI_TYPE_GET_EXTENT(MPI_LOGICAL, lb, size_logical, ierrmpi)

oldtypes(0) = MPI_LOGICAL
oldtypes(1) = MPI_INTEGER
oldtypes(2:8) = MPI_DOUBLE_PRECISION

data blocklengths/1,1,1,1,1,1,npayload,ndir,ndir/

offsets(0) = 0
offsets(1) = size_logical * blocklengths(0)
offsets(2) = offsets(1) + size_int * blocklengths(1)
offsets(3) = offsets(2) + size_double * blocklengths(2)
offsets(4) = offsets(3) + size_double * blocklengths(3)
offsets(5) = offsets(4) + size_double * blocklengths(4)
offsets(6) = offsets(5) + size_double * blocklengths(5)
offsets(7) = offsets(6) + size_double * blocklengths(6)
offsets(8) = offsets(7) + size_double * blocklengths(7)

call MPI_TYPE_STRUCT(9,blocklengths,offsets,oldtypes,type_particle,ierrmpi)
call MPI_TYPE_COMMIT(type_particle,ierrmpi)

end subroutine init_particles_com
!=============================================================================
subroutine finish_particles_com()

include 'amrvacdef.f'

call MPI_TYPE_FREE(type_particle,ierrmpi)

end subroutine finish_particles_com
!=============================================================================
subroutine finish_particles()

integer             :: ipart, iipart
!-----------------------------------------------------------------------------

call destroy_particles(nparticles_on_mype,particles_on_mype(1:nparticles_on_mype))

end subroutine finish_particles
!=============================================================================
subroutine set_neighbor_ipe()

include 'amrvacdef.f'
integer              :: igrid, iigrid,ipe
logical              :: ipe_is_neighbor(0:npe-1)
integer              :: my_neighbor_type, i^D
!-----------------------------------------------------------------------------

ipe_is_neighbor(:) = .false.
do iigrid=1,igridstail; igrid=igrids(iigrid);

   {do i^DB=-1,1\}
   if (i^D==0|.and.) then
      cycle
   end if
   my_neighbor_type=neighbor_type(i^D,igrid)
   
   select case (my_neighbor_type)
    case (1) ! boundary
       ! do nothing
    case (2) ! fine-coarse
       call ipe_fc(i^D,igrid,ipe_is_neighbor)
    case (3) ! same level
       call ipe_srl(i^D,igrid,ipe_is_neighbor)
    case (4) ! coarse-fine
       call ipe_cf(i^D,igrid,ipe_is_neighbor)
    end select
    
    {end do\}

end do

! remove self as neighbor
ipe_is_neighbor(mype) = .false.

! Now make the list of neighbors
npe_neighbors = 0
do ipe=0,npe-1
   if (ipe_is_neighbor(ipe)) then
      npe_neighbors = npe_neighbors + 1
      ipe_neighbor(npe_neighbors) = ipe
   end if
end do ! ipe

end subroutine set_neighbor_ipe
!=============================================================================
subroutine ipe_fc(i^D,igrid,ipe_is_neighbor)

include 'amrvacdef.f'
integer, intent(in) :: i^D, igrid
logical, intent(inout) :: ipe_is_neighbor(0:npe-1)
!-----------------------------------------------------------------------------

ipe_is_neighbor(neighbor(2,i^D,igrid)) = .true.

end subroutine ipe_fc
!=============================================================================
subroutine ipe_srl(i^D,igrid,ipe_is_neighbor)

include 'amrvacdef.f'
integer, intent(in) :: i^D, igrid
logical, intent(inout) :: ipe_is_neighbor(0:npe-1)
!-----------------------------------------------------------------------------

ipe_is_neighbor(neighbor(2,i^D,igrid)) = .true.

end subroutine ipe_srl
!=============================================================================
subroutine ipe_cf(i^D,igrid,ipe_is_neighbor)

include 'amrvacdef.f'
integer, intent(in)    :: i^D, igrid
logical, intent(inout) :: ipe_is_neighbor(0:npe-1)
integer                :: ic^D, inc^D
!-----------------------------------------------------------------------------

{do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
    inc^DB=2*i^DB+ic^DB
\}
ipe_is_neighbor( neighbor_child(2,inc^D,igrid) ) = .true.

{end do\}

end subroutine ipe_cf
!=============================================================================
subroutine check_particles_output()

include 'amrvacdef.f'

integer                         :: ipart,iipart
type(particle_node), dimension(0:nparticles_per_cpu_hi-1)  :: send_particles
integer                         :: send_n_particles_for_output
integer                         :: nout
double precision                :: tout
!-----------------------------------------------------------------------------

! is it time for output of individual particle?
if ((it_particles == itsavelast_particles + ditsave_particles) &
     .or.  (it_particles == 0)) call output_individual()

! append to ensemble files if its time to do so
send_n_particles_for_output = 0

do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);

   ! corresponding output slot:
   nout = nint(particle(ipart)%self%t/dtsave_ensemble)
   tout = dble(nout) * dtsave_ensemble
   ! should the particle be dumped?
   if (particle(ipart)%self%t .le. tout &
        .and. particle(ipart)%self%t+particle(ipart)%self%dt .gt. tout) then
      ! have to send particle to rank zero for output
      send_n_particles_for_output = send_n_particles_for_output + 1
      send_particles(send_n_particles_for_output-1) = particle(ipart)%self
   end if ! time to output?

end do

call output_ensemble(send_n_particles_for_output,send_particles(0:send_n_particles_for_output-1),'ensemble')

end subroutine check_particles_output
!=============================================================================
character(len=128) function make_particle_filename(tout,index,type)

include 'amrvacdef.f'

character(len=*), intent(in)    :: type
double precision, intent(in)    :: tout
integer, intent(in)             :: index
integer                         :: nout, mysnapshot
!-----------------------------------------------------------------------------

if (snapshotini .ne. -1) then
   mysnapshot = snapshotini
else
   mysnapshot = 0
end if

if (.not. time_advance) then
  select case(type)
  case ('ensemble')
     nout = nint(tout / dtsave_ensemble)
     write(make_particle_filename,"(a,i4.4,a,i6.6,a)") trim(filenameout),mysnapshot,'_ensemble',nout,'.csv'
  case ('destroy')
     write(make_particle_filename,"(a,i4.4,a)") trim(filenameout),mysnapshot,'_destroyed.csv'
  case ('individual')
     write(make_particle_filename,"(a,i4.4,a,i6.6,a)") trim(filenameout),mysnapshot,'_particle',index,'.csv'
  case ('followed')
     nout = nint(tout / dtsave_ensemble)
     write(make_particle_filename,"(a,i4.4,a,i6.6,a)") trim(filenameout),mysnapshot,'_followed',nout,'.csv'
  end select
else
  select case(type)
  case ('ensemble')
     nout = nint(tout / dtsave_ensemble)
     write(make_particle_filename,"(a,a,i6.6,a)") trim(filenameout),'_ensemble',nout,'.csv'
  case ('destroy')
     write(make_particle_filename,"(a,a)") trim(filenameout),'_destroyed.csv'
  case ('individual')
     write(make_particle_filename,"(a,a,i6.6,a)") trim(filenameout),'_particle',index,'.csv'
  case ('followed')
     nout = nint(tout / dtsave_ensemble)
     write(make_particle_filename,"(a,a,i6.6,a)") trim(filenameout),'_followed',nout,'.csv'
  end select
end if

end function make_particle_filename
!=============================================================================
subroutine output_ensemble(send_n_particles_for_output,send_particles,type)

include 'amrvacdef.f'

integer, intent(in)             :: send_n_particles_for_output
type(particle_node), dimension(0:send_n_particles_for_output-1), intent(in)  :: send_particles
character(len=*), intent(in)    :: type
character(len=128)              :: filename, filename2
type(particle_node), dimension(0:nparticles_per_cpu_hi-1)  :: receive_particles
integer                         :: status(MPI_STATUS_SIZE)
integer,dimension(0:npe-1)      :: receive_n_particles_for_output_from_ipe
integer                         :: ipe, ipart
!-----------------------------------------------------------------------------
receive_n_particles_for_output_from_ipe(:) = 0

if (npe==1) then 
   do ipart=1,send_n_particles_for_output
      filename = make_particle_filename(send_particles(ipart-1)%t,send_particles(ipart-1)%index,type)
      call output_particle(send_particles(ipart-1),mype,filename)
      if(type=='ensemble' .and. send_particles(ipart-1)%follow) then
        filename2 = make_particle_filename(send_particles(ipart-1)%t,send_particles(ipart-1)%index,'followed')
        call output_particle(send_particles(ipart-1),mype,filename2)
      end if
   end do
   return
end if


if (mype .ne. 0) then
   call MPI_SEND(send_n_particles_for_output,1,MPI_INTEGER,0,mype,icomm,ierrmpi)
else
   receive_n_particles_for_output_from_ipe(0) = send_n_particles_for_output
   if (send_n_particles_for_output .gt. 0) &
   receive_particles(0:send_n_particles_for_output-1) = &
        send_particles(0:send_n_particles_for_output-1)

   do ipe=1,npe-1
      call MPI_RECV(receive_n_particles_for_output_from_ipe(ipe),1,MPI_INTEGER,ipe,ipe,icomm,status,ierrmpi)
   end do

end if



if (mype .ne. 0) &
     call MPI_SEND(send_particles,send_n_particles_for_output,type_particle,0,mype,icomm,ierrmpi)

if (mype==0) then
   do ipart=1,receive_n_particles_for_output_from_ipe(0)
      filename = make_particle_filename(receive_particles(ipart-1)%t,receive_particles(ipart-1)%index,type)
      call output_particle(receive_particles(ipart-1),0,filename)
      if(type=='ensemble' .and. receive_particles(ipart-1)%follow) then
        filename2 = make_particle_filename(receive_particles(ipart-1)%t,receive_particles(ipart-1)%index,'followed')
        call output_particle(receive_particles(ipart-1),0,filename2)
      end if
   end do

   do ipe=1,npe-1
      call MPI_RECV(receive_particles,receive_n_particles_for_output_from_ipe(ipe),type_particle,ipe,ipe,icomm,status,ierrmpi)
      do ipart=1,receive_n_particles_for_output_from_ipe(ipe)
         filename = make_particle_filename(receive_particles(ipart-1)%t,receive_particles(ipart-1)%index,type)
         call output_particle(receive_particles(ipart-1),ipe,filename)
         if(type=='ensemble' .and. receive_particles(ipart-1)%follow) then
           filename2 = make_particle_filename(receive_particles(ipart-1)%t,receive_particles(ipart-1)%index,'followed')
           call output_particle(receive_particles(ipart-1),ipe,filename2)
         end if
      end do ! ipart
   end do ! ipe
end if ! mype == 0

end subroutine output_ensemble
!=============================================================================
subroutine output_individual()

include 'amrvacdef.f'
logical,parameter               :: output_from_root=.false.

character(len=128)              :: filename
integer                         :: ipart,iipart,ipe
integer                         :: send_n_particles_for_output
type(particle_node), dimension(0:nparticles_per_cpu_hi-1)  :: send_particles
type(particle_node), dimension(0:nparticles_per_cpu_hi-1)  :: receive_particles
integer                         :: status(MPI_STATUS_SIZE)
integer,dimension(0:npe-1)      :: receive_n_particles_for_output_from_ipe
!-----------------------------------------------------------------------------

send_n_particles_for_output = 0
receive_n_particles_for_output_from_ipe(:) = 0

do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);

   ! should the particle be dumped?
   if (particle(ipart)%self%follow) then
      ! have to send particle to rank zero for output
      send_n_particles_for_output = send_n_particles_for_output + 1
      send_particles(send_n_particles_for_output-1) = particle(ipart)%self
   end if ! follow
   
end do ! ipart


if (npe==1) then 
   do ipart=1,send_n_particles_for_output
      filename = make_particle_filename(send_particles(ipart-1)%t,send_particles(ipart-1)%index,'individual')
      call output_particle(send_particles(ipart-1),mype,filename)
   end do
      itsavelast_particles = it_particles
   return
end if

if (output_from_root) then

   if (mype .ne. 0) then
      call MPI_SEND(send_n_particles_for_output,1,MPI_INTEGER,0,mype,icomm,ierrmpi)
   else
      receive_n_particles_for_output_from_ipe(0) = send_n_particles_for_output
      receive_particles(0:send_n_particles_for_output) = &
           send_particles(0:send_n_particles_for_output)
      
      do ipe=1,npe-1
         call MPI_RECV(receive_n_particles_for_output_from_ipe(ipe),1,MPI_INTEGER,ipe,ipe,icomm,status,ierrmpi)
      end do

   end if


   if (mype .ne. 0) &
        call MPI_SEND(send_particles,send_n_particles_for_output,type_particle,0,mype,icomm,ierrmpi)
   
   if (mype==0) then
      do ipart=1,receive_n_particles_for_output_from_ipe(0)
      filename = make_particle_filename(receive_particles(ipart-1)%t,receive_particles(ipart-1)%index,'individual')
         call output_particle(receive_particles(ipart-1),0,filename)
      end do
      
      do ipe=1,npe-1
         call MPI_RECV(receive_particles,receive_n_particles_for_output_from_ipe(ipe),type_particle,ipe,ipe,icomm,status,ierrmpi)
         do ipart=1,receive_n_particles_for_output_from_ipe(ipe)
            filename = make_particle_filename(receive_particles(ipart-1)%t,receive_particles(ipart-1)%index,'individual')
            call output_particle(receive_particles(ipart-1),ipe,filename)
         end do
      end do
   end if

else ! output_from_root

   do ipart=1,send_n_particles_for_output
      ! generate filename
      filename = make_particle_filename(send_particles(ipart-1)%t,send_particles(ipart-1)%index,'individual')
      call output_particle(send_particles(ipart-1),mype,filename)
   end do

end if

itsavelast_particles = it_particles
end subroutine output_individual
!=============================================================================
subroutine output_particle(myparticle,ipe,filename)

include 'amrvacdef.f'
type(particle_node),intent(in)                 :: myparticle
integer, intent(in)                            :: ipe
character(len=128),intent(in)                  :: filename
logical,save                                   :: file_exists=.false.
integer,parameter                              :: ndoubles=2+^NC+^NC+npayload
character(len=20)                              :: formatstring
double precision                               :: x(^NC)
character(len=1024)                            :: line, data
integer                                        :: icomp
double precision, parameter                    :: minvalue = 1.0d-99
double precision                               :: roundoff
!-----------------------------------------------------------------------------

! normalize the position
if (typeaxial == 'slab') x = myparticle%x*UNIT_LENGTH
if (typeaxial == 'cylindrical') then
   x(r_)   = myparticle%x(r_)*UNIT_LENGTH
   x(zz_)   = myparticle%x(zz_)*UNIT_LENGTH
   x(pphi_) = myparticle%x(pphi_) 
end if
if (typeaxial == 'spherical') then
   {x(^C) = myparticle%x(^C)\}
   x(1) = x(1)*UNIT_LENGTH
end if


INQUIRE(FILE=filename, EXIST=file_exists)

if (.not. file_exists) then
   open(unit=unitparticles,file=filename,status='unknown',access='append')
   line=''
   do icomp=1, npayload
      write(data,"(a,i2.2,a)") 'payload',icomp,','
      line = trim(line)//trim(data)
   end do
   write(unitparticles,"(a,a,a)") 't,dt,x^C,u^C,',trim(line),'ipe,iteration,index'
else
   open(unit=unitparticles,file=filename,status='unknown',access='append')
end if

! create the formatstring:
!write(formatstring,"(a,i1,a)") '(',ndoubles,'(es14.6),3(i12.10))'

!write(unitparticles,formatstring) myparticle%t, myparticle%dt, x, &
!           myparticle%u, myparticle%payload(1:npayload), ipe, it_particles, myparticle%index

line = ''
write(data,"(es13.6, a)")roundoff(myparticle%t,minvalue), ','
line = trim(line)//trim(data)
write(data,"(es13.6, a)")roundoff(myparticle%dt,minvalue), ','
line = trim(line)//trim(data)
do icomp = 1, ^NC
   write(data,"(es13.6, a)")roundoff(x(icomp),minvalue), ','
   line = trim(line)//trim(data)
end do
do icomp = 1, ^NC
   write(data,"(es13.6, a)")roundoff(myparticle%u(icomp),minvalue), ','
   line = trim(line)//trim(data)
end do
do icomp = 1, npayload
   write(data,"(es13.6, a)")roundoff(myparticle%payload(icomp),minvalue), ','
   line = trim(line)//trim(data)
end do
write(data,"(i8.7, a)")ipe, ','
line = trim(line)//trim(data)
write(data,"(i11.10, a)")it_particles, ','
line = trim(line)//trim(data)
write(data,"(i8.7)")myparticle%index
line = trim(line)//trim(data)

write(unitparticles,"(a)") trim(line)

close(unit=unitparticles)
end subroutine output_particle
!=============================================================================
subroutine append_to_snapshot(myparticle)

type(particle_node),intent(in)                 :: myparticle
!-----------------------------------------------------------------------------


write(unitparticles) myparticle%index
write(unitparticles) myparticle%follow
write(unitparticles) myparticle%q
write(unitparticles) myparticle%m
write(unitparticles) myparticle%t
write(unitparticles) myparticle%dt
write(unitparticles) myparticle%x
write(unitparticles) myparticle%u 
write(unitparticles) myparticle%payload(1:npayload)


end subroutine append_to_snapshot
!=============================================================================
subroutine read_from_snapshot()

integer                         :: index,igrid_particle,ipe_particle
!-----------------------------------------------------------------------------


read(unitparticles) index
allocate(particle(index)%self)
particle(index)%self%index = index

read(unitparticles) particle(index)%self%follow
read(unitparticles) particle(index)%self%q
read(unitparticles) particle(index)%self%m
read(unitparticles) particle(index)%self%t
read(unitparticles) particle(index)%self%dt
read(unitparticles) particle(index)%self%x
read(unitparticles) particle(index)%self%u
read(unitparticles) particle(index)%self%payload(1:npayload)

if (particle(index)%self%follow) print*, 'follow index:', index

call find_particle_ipe(particle(index)%self%x,igrid_particle,ipe_particle)
particle(index)%igrid = igrid_particle
particle(index)%ipe   = ipe_particle

call push_particle_into_particles_on_mype(index)

end subroutine read_from_snapshot
!=============================================================================
subroutine comm_particles()

include 'amrvacdef.f'

integer                         :: ipart, iipart, igrid_particle, ipe_particle, ipe, iipe
integer                         :: index
integer                         :: tag_send, tag_receive, send_buff, rcv_buff
integer                         :: status(MPI_STATUS_SIZE)
integer, dimension(0:npe-1)     :: send_n_particles_to_ipe
integer, dimension(0:npe-1)    :: receive_n_particles_from_ipe
type(particle_node), dimension(0:nparticles_per_cpu_hi-1)  :: send_particles
type(particle_node), dimension(0:nparticles_per_cpu_hi-1)  :: receive_particles
integer, dimension(nparticles_per_cpu_hi,0:npe-1)  :: particle_index_to_be_sent_to_ipe
integer, dimension(nparticles_per_cpu_hi) :: particle_index_to_be_destroyed
integer                                   :: destroy_n_particles_mype
logical                                   :: BC_applied
!-----------------------------------------------------------------------------

send_n_particles_to_ipe(:)      = 0
receive_n_particles_from_ipe(:) = 0
destroy_n_particles_mype        = 0

! check if and where to send each particle, destroy if necessary
do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);

   ! first check if the particle should be destroyed (user defined criterion)
   if ( user_destroy(particle(ipart)%self) .or. &
        (.not.time_advance .and. particle(ipart)%self%t .gt. tmax_particles) ) then
      destroy_n_particles_mype  = destroy_n_particles_mype + 1
      particle_index_to_be_destroyed(destroy_n_particles_mype) = ipart
      cycle   
   end if

   ! is my particle still in the same igrid?
   if (.not.particle_in_igrid(ipart,particle(ipart)%igrid)) then
      call find_particle_ipe(particle(ipart)%self%x,igrid_particle,ipe_particle)
      
      ! destroy particle if out of domain (signalled by return value of -1)
      if (igrid_particle == -1 )then 
         call apply_periodB(particle(ipart)%self,igrid_particle,ipe_particle,BC_applied)
         if (.not. BC_applied .or. igrid_particle == -1) then
            destroy_n_particles_mype  = destroy_n_particles_mype + 1
            particle_index_to_be_destroyed(destroy_n_particles_mype) = ipart
            cycle
         end if
      end if
      
      ! particle still present
      particle(ipart)%igrid = igrid_particle
      particle(ipart)%ipe = ipe_particle
      
      ! if we have more than one core, is it on another cpu?
      if (npe .gt. 1 .and. particle(ipart)%ipe .ne. mype) then
         send_n_particles_to_ipe(ipe_particle) = &
              send_n_particles_to_ipe(ipe_particle) + 1
         particle_index_to_be_sent_to_ipe(send_n_particles_to_ipe(ipe_particle),ipe_particle) = ipart
      end if ! ipe_particle
      
   end if ! particle_in_grid

end do ! ipart

call destroy_particles(destroy_n_particles_mype,particle_index_to_be_destroyed(1:destroy_n_particles_mype))

! get out when only one core:
if (npe == 1) return


! communicate amount of particles to be sent/received
do iipe=1,npe_neighbors;ipe=ipe_neighbor(iipe);
   tag_send    = mype * npe + ipe
   tag_receive = ipe * npe + mype
   call MPI_SEND(send_n_particles_to_ipe(ipe),1,MPI_INTEGER,ipe,tag_send,icomm,ierrmpi)
   call MPI_RECV(receive_n_particles_from_ipe(ipe),1,MPI_INTEGER,ipe,tag_receive,icomm,status,ierrmpi)
end do


! send and receive the data of the particles
do iipe=1,npe_neighbors;ipe=ipe_neighbor(iipe);
   tag_send    = mype * npe + ipe
   tag_receive = ipe * npe + mype

   ! should i send some particles to ipe?
   if (send_n_particles_to_ipe(ipe) .gt. 0) then

      ! create the send buffer
      do ipart = 1, send_n_particles_to_ipe(ipe)
         send_particles(ipart-1) = particle(particle_index_to_be_sent_to_ipe(ipart,ipe))%self
      end do ! ipart

      call MPI_SEND(send_particles,send_n_particles_to_ipe(ipe),type_particle,ipe,tag_send,icomm,ierrmpi)

      do ipart = 1, send_n_particles_to_ipe(ipe)
         deallocate(particle(particle_index_to_be_sent_to_ipe(ipart,ipe))%self)
         call pull_particle_from_particles_on_mype(particle_index_to_be_sent_to_ipe(ipart,ipe))
      end do ! ipart

   end if ! send .gt. 0
  
   ! should i receive some particles from ipe?
   if (receive_n_particles_from_ipe(ipe) .gt. 0) then

      call MPI_RECV(receive_particles,receive_n_particles_from_ipe(ipe),type_particle,ipe,tag_receive,icomm,status,ierrmpi)

      do ipart = 1, receive_n_particles_from_ipe(ipe)

         index = receive_particles(ipart-1)%index
         allocate(particle(index)%self)
         particle(index)%self = receive_particles(ipart-1)
         call push_particle_into_particles_on_mype(index)

         ! since we don't send the igrid, need to re-locate it
         call find_particle_ipe(particle(index)%self%x,igrid_particle,ipe_particle)
         particle(index)%igrid = igrid_particle
         particle(index)%ipe = ipe_particle

      end do ! ipart

   end if ! receive .gt. 0
end do ! ipe


end subroutine comm_particles
!=============================================================================
subroutine comm_particles_global()

include 'amrvacdef.f'

integer                         :: ipart, iipart, igrid_particle, ipe_particle, ipe, iipe
integer                         :: index
integer                         :: tag_send, tag_receive, send_buff, rcv_buff
integer                         :: status(MPI_STATUS_SIZE)
integer, dimension(0:npe-1)     :: send_n_particles_to_ipe
integer, dimension(0:npe-1)    :: receive_n_particles_from_ipe
type(particle_node), dimension(0:nparticles_per_cpu_hi-1)  :: send_particles
type(particle_node), dimension(0:nparticles_per_cpu_hi-1)  :: receive_particles
integer, dimension(nparticles_per_cpu_hi,0:npe-1)  :: particle_index_to_be_sent_to_ipe
integer, dimension(nparticles_per_cpu_hi) :: particle_index_to_be_destroyed
integer                                   :: destroy_n_particles_mype
logical                                   :: BC_applied
!-----------------------------------------------------------------------------
send_n_particles_to_ipe(:)      = 0
receive_n_particles_from_ipe(:) = 0
destroy_n_particles_mype        = 0

! check if and where to send each particle, relocate all of them (in case grid has changed)
do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);

   call find_particle_ipe(particle(ipart)%self%x,igrid_particle,ipe_particle)
      
   particle(ipart)%igrid = igrid_particle
   particle(ipart)%ipe = ipe_particle
   
      ! if we have more than one core, is it on another cpu?
   if (particle(ipart)%ipe .ne. mype) then
      send_n_particles_to_ipe(ipe_particle) = &
           send_n_particles_to_ipe(ipe_particle) + 1
      particle_index_to_be_sent_to_ipe(send_n_particles_to_ipe(ipe_particle),ipe_particle) = ipart
   end if ! ipe_particle
      
end do ! ipart

! get out when only one core:
if (npe == 1) return

! communicate amount of particles to be sent/received
do ipe=0,npe-1;if (ipe .eq. mype) cycle;
   tag_send    = mype * npe + ipe
   tag_receive = ipe * npe + mype
   call MPI_SEND(send_n_particles_to_ipe(ipe),1,MPI_INTEGER,ipe,tag_send,icomm,ierrmpi)
   call MPI_RECV(receive_n_particles_from_ipe(ipe),1,MPI_INTEGER,ipe,tag_receive,icomm,status,ierrmpi)
end do


! send and receive the data of the particles
do ipe=0,npe-1;if (ipe .eq. mype) cycle;
   tag_send    = mype * npe + ipe
   tag_receive = ipe * npe + mype

   ! should i send some particles to ipe?
   if (send_n_particles_to_ipe(ipe) .gt. 0) then

      ! create the send buffer
      do ipart = 1, send_n_particles_to_ipe(ipe)
         send_particles(ipart-1) = particle(particle_index_to_be_sent_to_ipe(ipart,ipe))%self
      end do ! ipart

      call MPI_SEND(send_particles,send_n_particles_to_ipe(ipe),type_particle,ipe,tag_send,icomm,ierrmpi)

      do ipart = 1, send_n_particles_to_ipe(ipe)
         deallocate(particle(particle_index_to_be_sent_to_ipe(ipart,ipe))%self)
         call pull_particle_from_particles_on_mype(particle_index_to_be_sent_to_ipe(ipart,ipe))
      end do ! ipart

   end if ! send .gt. 0
  
   ! should i receive some particles from ipe?
   if (receive_n_particles_from_ipe(ipe) .gt. 0) then

      call MPI_RECV(receive_particles,receive_n_particles_from_ipe(ipe),type_particle,ipe,tag_receive,icomm,status,ierrmpi)

      do ipart = 1, receive_n_particles_from_ipe(ipe)

         index = receive_particles(ipart-1)%index
         allocate(particle(index)%self)
         particle(index)%self = receive_particles(ipart-1)
         call push_particle_into_particles_on_mype(index)

         ! since we don't send the igrid, need to re-locate it
         call find_particle_ipe(particle(index)%self%x,igrid_particle,ipe_particle)
         particle(index)%igrid = igrid_particle
         particle(index)%ipe = ipe_particle

      end do ! ipart

   end if ! receive .gt. 0
end do ! ipe


end subroutine comm_particles_global
!=============================================================================
subroutine apply_periodB(particle,igrid_particle,ipe_particle,BC_applied)

include 'amrvacdef.f'
type(particle_node), intent(inout)        :: particle
integer, intent(inout)                    :: igrid_particle, ipe_particle
logical,intent(out)                       :: BC_applied
integer                                   :: idim
!-----------------------------------------------------------------------------
BC_applied = .false.
! get out if we don't have any periodic BC:
if (.not. any(periodB(1:ndim))) return

! go through dimensions and try re-inject the particle at the other side
do idim=1,ndim

   if (.not. periodB(idim)) cycle

   select case(idim)
      {case (^D)
      if (particle%x(^D) .lt. xprobmin^D) then
         particle%x(^D) = particle%x(^D) + (xprobmax^D - xprobmin^D)
         BC_applied = .true.
      end if
      if (particle%x(^D) .ge. xprobmax^D) then
         particle%x(^D) = particle%x(^D) - (xprobmax^D - xprobmin^D)
         BC_applied = .true.
      end if
      \}
   end select

end do

call find_particle_ipe(particle%x,igrid_particle,ipe_particle)

end subroutine apply_periodB
!=============================================================================
subroutine destroy_particles(destroy_n_particles_mype,particle_index_to_be_destroyed)
! clean up destroyed particles on all cores

include 'amrvacdef.f'

integer, intent(in)                                   :: destroy_n_particles_mype
integer, dimension(1:destroy_n_particles_mype), intent(in) :: particle_index_to_be_destroyed
type(particle_node), dimension(0:destroy_n_particles_mype-1):: destroy_particles_mype
integer                                               :: iipart,ipart,destroy_n_particles
!-----------------------------------------------------------------------------
destroy_n_particles             = 0

! append the particle to list of destroyed particles
do iipart=1,destroy_n_particles_mype;ipart=particle_index_to_be_destroyed(iipart);
   destroy_particles_mype(iipart-1) = particle(ipart)%self
end do

call output_ensemble(destroy_n_particles_mype,destroy_particles_mype,'destroy')


if (npe > 1) then
   call MPI_ALLREDUCE(destroy_n_particles_mype,destroy_n_particles,1,MPI_INTEGER,MPI_SUM,icomm,ierrmpi)
else
   destroy_n_particles = destroy_n_particles_mype
end if

nparticles = nparticles - destroy_n_particles

do iipart=1,destroy_n_particles_mype;ipart=particle_index_to_be_destroyed(iipart);

   particle(ipart)%igrid = -1
   particle(ipart)%ipe   = -1
   if(time_advance) then
     write(*,*), particle(ipart)%self%t, ': particle',ipart,' has left at it ',it_particles,' on pe', mype
     write(*,*), particle(ipart)%self%t, ': particles remaining:',nparticles, '(total)', nparticles_on_mype-1, 'on pe', mype
   endif
   deallocate(particle(ipart)%self)
   call pull_particle_from_particles_on_mype(ipart)
end do

end subroutine destroy_particles
!=============================================================================
subroutine push_particle_into_particles_on_mype(ipart)

implicit none

integer, intent(in)            :: ipart
!-----------------------------------------------------------------------------

nparticles_on_mype = nparticles_on_mype + 1
particles_on_mype(nparticles_on_mype) = ipart

end subroutine push_particle_into_particles_on_mype
!=============================================================================
subroutine pull_particle_from_particles_on_mype(ipart)

implicit none

integer, intent(in)            :: ipart
integer                        :: i
!-----------------------------------------------------------------------------


do i=1,nparticles_on_mype
   if (particles_on_mype(i) == ipart) then
      particles_on_mype(i) = particles_on_mype(nparticles_on_mype)
      exit
   end if
end do

nparticles_on_mype = nparticles_on_mype - 1

end subroutine pull_particle_from_particles_on_mype
!=============================================================================
end module mod_particles
!=============================================================================
!=============================================================================
}
