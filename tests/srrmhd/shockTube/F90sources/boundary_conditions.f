!=============================================================================
subroutine bc_phys(iside,idims,time,w,x,ixGmin1,ixGmax1,ixBmin1,ixBmax1)

include 'amrvacdef.f'

integer, intent(in) :: iside, idims, ixGmin1,ixGmax1,ixBmin1,ixBmax1
double precision, intent(in) :: time
double precision, intent(inout) :: w(ixGmin1:ixGmax1,1:nw)
double precision, intent(in) :: x(ixGmin1:ixGmax1,1:ndim)

integer :: iw, iB, ix1, ixImin1,ixImax1
!-----------------------------------------------------------------------------
select case (idims)
case (1)
   if (iside==2) then
      ! maximal boundary
      iB=ismax1
      ixImin1=ixBmax1+1-dixB;
      ixImax1=ixBmax1;
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux
         select case (typeB(iw,iB))
         case ("symm")
            w(ixImin1:ixImax1,iw) = w(ixImin1-1:ixImin1-dixB:-1,iw)
         case ("asymm")
            w(ixImin1:ixImax1,iw) =-w(ixImin1-1:ixImin1-dixB:-1,iw)
         case ("cont")
            do ix1=ixImin1,ixImax1
               w(ix1,iw) = w(ixImin1-1,iw)
            end do
         case("noinflow")
            if (iw==1+1)then
              do ix1=ixImin1,ixImax1
                  w(ix1,iw) = max(w(ixImin1-1,iw),zero)
              end do
            else
              do ix1=ixImin1,ixImax1
                  w(ix1,iw) = w(ixImin1-1,iw)
              end do
            end if
         case("limitinflow")
            if (iw==1+1)then
              do ix1=ixImin1,ixImax1
                  w(ix1,iw) = max(w(ixImin1-1,iw), &
                                           w(ixImin1-1,iw)*ratebdflux)
              end do
            else
              do ix1=ixImin1,ixImax1
                  w(ix1,iw) = w(ixImin1-1,iw)
              end do
            end if
         case ("special")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
            w(ixImin1:ixImax1,iw) = - w(ixImin1:ixImax1,iw)
         case ("periodic")
!            call mpistop("periodic bc info should come from neighbors")
         case default
            write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
               "for variable iw=",iw," and side iB=",iB
         end select
      end do
   else
      ! minimal boundary
      iB=ismin1
      ixImin1=ixBmin1;
      ixImax1=ixBmin1-1+dixB;
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux
         select case (typeB(iw,iB))
         case ("symm")
            w(ixImin1:ixImax1,iw) = w(ixImax1+dixB:ixImax1+1:-1,iw)
         case ("asymm")
            w(ixImin1:ixImax1,iw) =-w(ixImax1+dixB:ixImax1+1:-1,iw)
         case ("cont")
            do ix1=ixImin1,ixImax1
               w(ix1,iw) = w(ixImax1+1,iw)
            end do
         case("noinflow")
            if (iw==1+1)then
               do ix1=ixImin1,ixImax1
                 w(ix1,iw) = min(w(ixImax1+1,iw),zero)
               end do
            else
               do ix1=ixImin1,ixImax1
                 w(ix1,iw) = w(ixImax1+1,iw)
               end do
            end if
         case("limitinflow")
            if (iw==1+1)then
               do ix1=ixImin1,ixImax1
                 w(ix1,iw) = min(w(ixImax1+1,iw), &
                                          w(ixImax1+1,iw)*ratebdflux)
               end do
            else
               do ix1=ixImin1,ixImax1
                 w(ix1,iw) = w(ixImax1+1,iw)
               end do
            end if
         case ("special")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
            w(ixImin1:ixImax1,iw) = - w(ixImin1:ixImax1,iw)
         case ("periodic")
!            call mpistop("periodic bc info should come from neighbors")
         case default
            write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
               "for variable iw=",iw," and side iB=",iB
         end select
      end do
   end if 
end select

! do special case AFTER all normal cases are set
!do iw=1,nwflux+nwaux
! opedit: iw==0 since this breaks fewest of setups.
   if (any(typeB(1:nwflux+nwaux,iB)=="special")) then
      call specialbound_usr(time,ixGmin1,ixGmax1,ixImin1,ixImax1,0,iB,w,x)
   end if
!end do

end subroutine bc_phys
!=============================================================================
subroutine getintbc(time,ixGmin1,ixGmax1,pwuse)

include 'amrvacdef.f'

double precision, intent(in)   :: time
integer, intent(in)            :: ixGmin1,ixGmax1
type(walloc), dimension(ngridshi)          :: pwuse

! .. local ..
integer :: iigrid, igrid, ixOmin1,ixOmax1,level
!----------------------------------------------------------------------------
ixOmin1=ixGmin1+dixB;ixOmax1=ixGmax1-dixB;

do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
!do iigrid=1,igridstail; igrid=igrids(iigrid);
   dxlevel(1)=rnode(rpdx1_,igrid);
   if (.not.slab) mygeo => pgeo(igrid)
   if (B0field) then
      myB0_cell => pB0_cell(igrid)
      myB0_face1 => pB0_face1(igrid)
   end if
   typelimiter=typelimiter1(node(plevel_,igrid))
   typegradlimiter=typegradlimiter1(node(plevel_,igrid))
   level=node(plevel_,igrid)
   saveigrid=igrid
   call bc_int(level,time,ixGmin1,ixGmax1,ixOmin1,ixOmax1,pwuse(igrid)%w,&
      px(igrid)%x)
end do

      
end subroutine getintbc
!=============================================================================
