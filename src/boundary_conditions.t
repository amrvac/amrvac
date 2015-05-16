!=============================================================================
subroutine bc_phys(iside,idims,time,w,x,ixG^L,ixB^L)

include 'amrvacdef.f'

integer, intent(in) :: iside, idims, ixG^L,ixB^L
double precision, intent(in) :: time
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision, intent(in) :: x(ixG^S,1:ndim)

integer :: iw, iB, ix^D, ixI^L
!-----------------------------------------------------------------------------
select case (idims)
{case (^D)
   if (iside==2) then
      ! maximal boundary
      iB=ismax^D
      ixImin^DD=ixBmax^D+1-dixB^D%ixImin^DD=ixBmin^DD;
      ixImax^DD=ixBmax^DD;
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux
         select case (typeB(iw,iB))
         case ("symm")
            w(ixI^S,iw) = w(ixImin^D-1:ixImin^D-dixB:-1^D%ixI^S,iw)
         case ("asymm")
            w(ixI^S,iw) =-w(ixImin^D-1:ixImin^D-dixB:-1^D%ixI^S,iw)
         case ("cont")
            do ix^D=ixImin^D,ixImax^D
               w(ix^D^D%ixI^S,iw) = w(ixImin^D-1^D%ixI^S,iw)
            end do
         case("noinflow")
            if (iw==1+^D)then
              do ix^D=ixImin^D,ixImax^D
                  w(ix^D^D%ixI^S,iw) = max(w(ixImin^D-1^D%ixI^S,iw),zero)
              end do
            else
              do ix^D=ixImin^D,ixImax^D
                  w(ix^D^D%ixI^S,iw) = w(ixImin^D-1^D%ixI^S,iw)
              end do
            end if
         case("limitinflow")
            if (iw==1+^D)then
              do ix^D=ixImin^D,ixImax^D
                  w(ix^D^D%ixI^S,iw) = max(w(ixImin^D-1^D%ixI^S,iw), &
                                           w(ixImin^D-1^D%ixI^S,iw)*ratebdflux)
              end do
            else
              do ix^D=ixImin^D,ixImax^D
                  w(ix^D^D%ixI^S,iw) = w(ixImin^D-1^D%ixI^S,iw)
              end do
            end if
         case ("special")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("aperiodic")
            !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
            w(ixI^S,iw) = - w(ixI^S,iw)
         case ("periodic")
!            call mpistop("periodic bc info should come from neighbors")
         case default
            write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
               "for variable iw=",iw," and side iB=",iB
         end select
      end do
   else
      ! minimal boundary
      iB=ismin^D
      ixImin^DD=ixBmin^DD;
      ixImax^DD=ixBmin^D-1+dixB^D%ixImax^DD=ixBmax^DD;
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux
         select case (typeB(iw,iB))
         case ("symm")
            w(ixI^S,iw) = w(ixImax^D+dixB:ixImax^D+1:-1^D%ixI^S,iw)
         case ("asymm")
            w(ixI^S,iw) =-w(ixImax^D+dixB:ixImax^D+1:-1^D%ixI^S,iw)
         case ("cont")
            do ix^D=ixImin^D,ixImax^D
               w(ix^D^D%ixI^S,iw) = w(ixImax^D+1^D%ixI^S,iw)
            end do
         case("noinflow")
            if (iw==1+^D)then
               do ix^D=ixImin^D,ixImax^D
                 w(ix^D^D%ixI^S,iw) = min(w(ixImax^D+1^D%ixI^S,iw),zero)
               end do
            else
               do ix^D=ixImin^D,ixImax^D
                 w(ix^D^D%ixI^S,iw) = w(ixImax^D+1^D%ixI^S,iw)
               end do
            end if
         case("limitinflow")
            if (iw==1+^D)then
               do ix^D=ixImin^D,ixImax^D
                 w(ix^D^D%ixI^S,iw) = min(w(ixImax^D+1^D%ixI^S,iw), &
                                          w(ixImax^D+1^D%ixI^S,iw)*ratebdflux)
               end do
            else
               do ix^D=ixImin^D,ixImax^D
                 w(ix^D^D%ixI^S,iw) = w(ixImax^D+1^D%ixI^S,iw)
               end do
            end if
         case ("special")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("aperiodic")
            !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
            w(ixI^S,iw) = - w(ixI^S,iw)
         case ("periodic")
!            call mpistop("periodic bc info should come from neighbors")
         case default
            write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
               "for variable iw=",iw," and side iB=",iB
         end select
      end do
   end if \}
end select

! do special case AFTER all normal cases are set
!do iw=1,nwflux+nwaux
! opedit: iw==0 since this breaks fewest of setups.
   if (any(typeB(1:nwflux+nwaux,iB)=="special")) then
      call specialbound_usr(time,ixG^L,ixI^L,0,iB,w,x)
   end if
!end do

end subroutine bc_phys
!=============================================================================
subroutine getintbc(time,ixG^L,pwuse)

include 'amrvacdef.f'

double precision, intent(in)   :: time
integer, intent(in)            :: ixG^L
type(walloc), dimension(ngridshi)          :: pwuse

! .. local ..
integer :: iigrid, igrid, ixO^L,level
!----------------------------------------------------------------------------
ixO^L=ixG^L^LSUBdixB;

do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
!do iigrid=1,igridstail; igrid=igrids(iigrid);
   ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
   if (.not.slab) mygeo => pgeo(igrid)
   if (B0field) then
      myB0_cell => pB0_cell(igrid)
      {^D&myB0_face^D => pB0_face^D(igrid)\}
   end if
   typelimiter=typelimiter1(node(plevel_,igrid))
   typegradlimiter=typegradlimiter1(node(plevel_,igrid))
   level=node(plevel_,igrid)
   saveigrid=igrid
   call bc_int(level,time,ixG^L,ixO^L,pwuse(igrid)%w,px(igrid)%x)
end do

      
end subroutine getintbc
!=============================================================================
