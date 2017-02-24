!=============================================================================
subroutine bc_phys(iside,idims,time,qdt,w,x,ixG^L,ixB^L)
use mod_usr_methods, only: usr_special_bc
use mod_global_parameters

integer, intent(in) :: iside, idims, ixG^L,ixB^L
double precision, intent(in) :: time,qdt
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision :: wtmp(ixG^S,1:nwflux)

integer :: iw, iB, ix^D, ixI^L, ixM^L, nghostcellsi,iib^D
logical  :: isphysbound
!-----------------------------------------------------------------------------
select case (idims)
{case (^D)
   if (iside==2) then
      ! maximal boundary
      iB=ismax^D
      ixImin^DD=ixBmax^D+1-nghostcells^D%ixImin^DD=ixBmin^DD;
      ixImax^DD=ixBmax^DD;
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux
         select case (typeboundary(iw,iB))
         case ("symm")
            w(ixI^S,iw) = w(ixImin^D-1:ixImin^D-nghostcells:-1^D%ixI^S,iw)
         case ("asymm")
            w(ixI^S,iw) =-w(ixImin^D-1:ixImin^D-nghostcells:-1^D%ixI^S,iw)
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
         case ("special")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("character")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("aperiodic")
            !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
            w(ixI^S,iw) = - w(ixI^S,iw)
         case ("periodic")
!            call mpistop("periodic bc info should come from neighbors")
         case default
            write (unitterm,*) "Undefined boundarytype ",typeboundary(iw,iB), &
               "for variable iw=",iw," and side iB=",iB
         end select
      end do
   else
      ! minimal boundary
      iB=ismin^D
      ixImin^DD=ixBmin^DD;
      ixImax^DD=ixBmin^D-1+nghostcells^D%ixImax^DD=ixBmax^DD;
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux
         select case (typeboundary(iw,iB))
         case ("symm")
            w(ixI^S,iw) = w(ixImax^D+nghostcells:ixImax^D+1:-1^D%ixI^S,iw)
         case ("asymm")
            w(ixI^S,iw) =-w(ixImax^D+nghostcells:ixImax^D+1:-1^D%ixI^S,iw)
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
         case ("special")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("character")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("aperiodic")
            !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
            w(ixI^S,iw) = - w(ixI^S,iw)
         case ("periodic")
!            call mpistop("periodic bc info should come from neighbors")
         case default
            write (unitterm,*) "Undefined boundarytype ",typeboundary(iw,iB), &
               "for variable iw=",iw," and side iB=",iB
         end select
      end do
   end if \}
end select

! do special case AFTER all normal cases are set
!do iw=1,nwflux+nwaux
! opedit: iw==0 since this breaks fewest of setups.
if (any(typeboundary(1:nwflux+nwaux,iB)=="special")) then
   if (.not. associated(usr_special_bc)) &
        call mpistop("usr_special_bc not defined")
   call usr_special_bc(time,ixG^L,ixI^L,iB,w,x)
end if

{#IFDEF EVOLVINGBOUNDARY
if (any(typeboundary(1:nwflux,iB)=="character")) then
  ixM^L=ixM^LL;
  if(ixGmax1==ixGhi1) then
    nghostcellsi=nghostcells
  else
    nghostcellsi=ceiling(nghostcells*0.5d0)
  end if
  select case (idims)
  {case (^D)
     if (iside==2) then
        ! maximal boundary
        ixImin^DD=ixGmax^D+1-nghostcellsi^D%ixImin^DD=ixBmin^DD;
        ixImax^DD=ixBmax^DD;
        if(all(w(ixI^S,1:nwflux)==0.d0)) then
          do ix^D=ixImin^D,ixImax^D
             w(ix^D^D%ixI^S,1:nwflux) = w(ixImin^D-1^D%ixI^S,1:nwflux)
          end do
        end if
        if(qdt>0.d0.and.ixGmax^D==ixGhi^D) then
          ixImin^DD=ixImin^D^D%ixImin^DD=ixMmin^DD;
          ixImax^DD=ixImax^D^D%ixImax^DD=ixMmax^DD;
          wtmp(ixG^S,1:nw)=pw(saveigrid)%wold(ixG^S,1:nw)
          call characteristic_project(idims,iside,ixG^L,ixI^L,wtmp,x,dxlevel,qdt)
          w(ixI^S,1:nwflux)=wtmp(ixI^S,1:nwflux)
        end if
     else
        ! minimal boundary
        ixImin^DD=ixBmin^DD;
        ixImax^DD=ixGmin^D-1+nghostcellsi^D%ixImax^DD=ixBmax^DD;
        if(all(w(ixI^S,1:nwflux)==0.d0)) then
          do ix^D=ixImin^D,ixImax^D
             w(ix^D^D%ixI^S,1:nwflux) = w(ixImax^D+1^D%ixI^S,1:nwflux)
          end do
        end if
        if(qdt>0.d0.and.ixGmax^D==ixGhi^D) then
          ixImin^DD=ixImin^D^D%ixImin^DD=ixMmin^DD;
          ixImax^DD=ixImax^D^D%ixImax^DD=ixMmax^DD;
          wtmp(ixG^S,1:nw)=pw(saveigrid)%wold(ixG^S,1:nw)
          call characteristic_project(idims,iside,ixG^L,ixI^L,wtmp,x,dxlevel,qdt)
          w(ixI^S,1:nwflux)=wtmp(ixI^S,1:nwflux)
        end if
     end if \}
  end select
  if(ixGmax1==ixGhi1) then
    call identifyphysbound(saveigrid,isphysbound,iib^D)   
    if(iib1==-1.and.iib2==-1) then
      do ix2=nghostcells,1,-1 
        do ix1=nghostcells,1,-1 
          w(ix^D,1:nwflux)=(w(ix1+1,ix2+1,1:nwflux)+w(ix1+1,ix2,1:nwflux)+w(ix1,ix2+1,1:nwflux))/3.d0
        end do
      end do
    end if
    if(iib1== 1.and.iib2==-1) then
      do ix2=nghostcells,1,-1 
        do ix1=ixMmax1+1,ixGmax1
          w(ix^D,1:nwflux)=(w(ix1-1,ix2+1,1:nwflux)+w(ix1-1,ix2,1:nwflux)+w(ix1,ix2+1,1:nwflux))/3.d0
        end do
      end do
    end if
    if(iib1==-1.and.iib2== 1) then
      do ix2=ixMmax2+1,ixGmax2
        do ix1=nghostcells,1,-1 
          w(ix^D,1:nwflux)=(w(ix1+1,ix2-1,1:nwflux)+w(ix1+1,ix2,1:nwflux)+w(ix1,ix2-1,1:nwflux))/3.d0
        end do
      end do
    end if
    if(iib1== 1.and.iib2== 1) then
      do ix2=ixMmax2+1,ixGmax2
        do ix1=ixMmax1+1,ixGmax1
          w(ix^D,1:nwflux)=(w(ix1-1,ix2-1,1:nwflux)+w(ix1-1,ix2,1:nwflux)+w(ix1,ix2-1,1:nwflux))/3.d0
        end do
      end do
    end if
  end if
end if
}
!end do

end subroutine bc_phys
!=============================================================================
subroutine getintbc(time,ixG^L)
use mod_usr_methods, only: usr_internal_bc
use mod_global_parameters

double precision, intent(in)   :: time
integer, intent(in)            :: ixG^L

! .. local ..
integer :: iigrid, igrid, ixO^L,level
!----------------------------------------------------------------------------
ixO^L=ixG^L^LSUBnghostcells;

do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
!do iigrid=1,igridstail; igrid=igrids(iigrid);
   ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
   block=>pw(igrid)
   typelimiter=limiter(node(plevel_,igrid))
   typegradlimiter=gradient_limiter(node(plevel_,igrid))
   level=node(plevel_,igrid)
   saveigrid=igrid

   if (associated(usr_internal_bc)) then
      call usr_internal_bc(level,time,ixG^L,ixO^L,pw(igrid)%wb,pw(igrid)%x)
   end if
end do

end subroutine getintbc
!=============================================================================
