subroutine write_shape_io(igrid,ixO^L)
  use mod_global_parameters
  use mod_ghostcells_update, only: identifyphysbound

  integer, intent(in) :: igrid
  integer, intent(out) :: ixO^L

  integer :: iib^D
  logical :: isphysbound 

  ixO^L=ixG^LL^LSUBnghostcells;
  if(save_physical_boundary) then 
    call identifyphysbound(igrid,isphysbound,iib^D)
    {
    if(iib^D ==1) then
      ixOmax^D=ixOmax^D+nghostcells
    else if(iib^D==-1) then
      ixOmin^D=ixOmin^D-nghostcells
    end if
    \}
  end if
end subroutine write_shape_io
