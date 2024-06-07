!> Module to set boundary conditions from user data
module mod_bc_data
#include "amrvac.h"
  use mod_lookup_table
  use mod_global_parameters, only: std_len

  implicit none
  private

  integer, parameter :: max_boundary_conds = 10

  type bc_data_t
     integer                        :: n_variables
     double precision               :: origin(3)
     double precision               :: dx(3)
     integer                        :: n_points(3)
     character(len=40), allocatable :: names(:)
     double precision, allocatable  :: values(:, :, :, :)
  end type bc_data_t

  type(LT_t)  :: lt_1d(max_boundary_conds)
  type(LT2_t) :: lt_2d(max_boundary_conds)
  type(LT3_t) :: lt_3d(max_boundary_conds)

  !> Whether boundary condition data is time varying
  logical, public, protected :: bc_data_time_varying = .false.

  !> Integer array for indexing lookup tables per variable per direction
  integer, public, protected, allocatable :: bc_data_ix(:, :)

  !> data file name
  character(len=std_len), public, protected :: boundary_data_file_name
  logical, public, protected :: interp_phy_first_row=.true.
  logical, public, protected :: interp_phy_first_row_same=.false.
  logical, public, protected :: bc_phy_first_row=.false.
  logical, public, protected :: boundary_data_primitive=.false.

  public :: bc_data_init
  public :: bc_data_set
  public :: bc_data_get_2d
  public :: bc_data_get_3d

contains

  subroutine bc_data_init()
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    integer                :: i, iw, ib, n_files, n_bc
    double precision       :: xmax(3)
    type(bc_data_t)        :: bc

    call bc_read_params(par_files)

    allocate(bc_data_ix(nwfluxbc, 2*ndim))

    bc_data_ix(:, :) = -1
    n_bc             = 0

    if(any(typeboundary(1:nwfluxbc,1:2*ndim)==bc_data ) .or. any(typeboundary(1:nwfluxbc,1:2*ndim)==bc_icarus)) then
      call read_vtk_structured_points(trim(boundary_data_file_name), bc)
      xmax = bc%origin + (bc%n_points-1) * bc%dx
      bc_data_time_varying = (bc%n_points(ndim) > 1)

    endif

    do ib = 1, 2 * ndim
       do iw = 1, nwfluxbc
          if (typeboundary(iw, ib)==bc_data .or. typeboundary(iw, ib)==bc_icarus) then
             n_bc               = n_bc + 1
             bc_data_ix(iw, ib) = n_bc


             {^IFONED
             call mpistop("bc_data_init: 1D case not supported")
             }
             {^IFTWOD
             if (bc_data_time_varying) then
                lt_2d(n_bc) = LT2_create_from_data(bc%origin(1:ndim), &
                     xmax(1:ndim), bc%values(:, :, 1:1, n_bc))
             else
                ! Use first point in time
                lt_1d(n_bc) = LT_create_from_data(bc%origin(1), &
                     xmax(1), bc%values(:, 1, 1:1, n_bc))
             end if
             }
             {^IFTHREED
             if (bc_data_time_varying) then
                lt_3d(n_bc) = LT3_create_from_data(bc%origin(1:ndim), &
                     xmax(1:ndim), bc%values(:, :, :, n_bc:n_bc))
             else
                ! Use first point in time
                lt_2d(n_bc) = LT2_create_from_data(bc%origin(1:ndim-1), &
                     xmax(1:ndim-1), bc%values(:, :, 1:1, n_bc))
             end if
             }
          end if
       end do
    end do

  end subroutine bc_data_init

  !> Read this module"s parameters from a file
  subroutine bc_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /bd_list/ boundary_data_file_name, interp_phy_first_row, bc_phy_first_row,&
                       boundary_data_primitive, interp_phy_first_row_same

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, bd_list, end=111)
111    close(unitpar)
    end do
 
  end subroutine bc_read_params

  elemental function bc_data_get_3d(n_bc, x1, x2, qt) result(val)
    integer, intent(in)          :: n_bc
    double precision, intent(in) :: x1, x2, qt
    double precision             :: val

    if (bc_data_time_varying) then
       val = LT3_get_col(lt_3d(n_bc), 1, x1, x2, qt)
    else
#if defined(ICARUS_INTERP) && ICARUS_INTERP==1
       val = LT2_get_col_icarus(lt_2d(n_bc), 1, x1, x2)
#else
       val = LT2_get_col(lt_2d(n_bc), 1, x1, x2)
#endif
    end if
  end function bc_data_get_3d

  elemental function bc_data_get_2d(n_bc, x1, qt) result(val)
    integer, intent(in)          :: n_bc
    double precision, intent(in) :: x1, qt
    double precision             :: val

    if (bc_data_time_varying) then
       val = LT2_get_col(lt_2d(n_bc), 1, x1, qt)
    else
       val = LT_get_col(lt_1d(n_bc), 1, x1)
    end if
  end function bc_data_get_2d

  !> Set boundary conditions according to user data
  subroutine bc_data_set(qt,ixI^L,ixO^L,iB,w,x)
    use mod_global_parameters
    use mod_physics, only: phys_to_conserved,phys_to_primitive
    use mod_comm_lib, only: mpistop

    integer, intent(in)             :: ixI^L, ixO^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: tmp(ixO^S)
    integer                         :: i, ix, iw, n_bc

    integer :: ixOO^L

    ixOO^L=ixO^L;

    {^IFTWOD
    select case (iB)
    case (1, 2)
       if(interp_phy_first_row) then
            if (iB == 1) then
               ix = ixOmax1+1
            else
               ix = ixOmin1-1
            end if
            if(boundary_data_primitive) then
              call phys_to_primitive(ixI^L,ix,ixOmin2,ix,ixOmax2,w,x)
            endif   
       endif  
      ! the reason for this is that not all bc for all the vars 
      ! might be set with bc_data  
      if(bc_phy_first_row)then 
        if (iB == 1) then
          ixOOmax1=ixOOmax1+1 
        else
          ixOOmin1=ixOOmin1-1 
        endif
        if(boundary_data_primitive) then
          if (iB == 1) then
            ixOOmin1=ixOOmax1 
          else
            ixOOmax1=ixOOmin1
          endif
          call phys_to_primitive(ixI^L,ixOO^L,w,x)
          if (iB == 1) then
            ixOOmin1=ixOmin1 ! to be used later in to_conserved
          else
            ixOOmax1=ixOmax1
          endif
        endif
      endif
       do iw = 1, nwfluxbc
          n_bc = bc_data_ix(iw, iB)
          if (n_bc == -1) cycle
          
          tmp(ixOmin1, ixOmin2:ixOmax2) = bc_data_get_2d(n_bc, &
                  x(ixOmin1, ixOmin2:ixOmax2, 2), qt)
          if(interp_phy_first_row) then
            ! Approximate boundary value by linear interpolation to first ghost
            ! cell, rest of ghost cells contains the same value
            if(interp_phy_first_row_same) then
              do i = 0, ixOOmax1-ixOOmin1
                 w(ixOOmin1+i, ixOmin2:ixOmax2, iw) = &
                      2 * tmp(ixOmin1, ixOmin2:ixOmax2) - &
                      w(ix, ixOmin2:ixOmax2, iw)
              end do
            else
              if (iB==1) then
                do i = 0, ixOOmax1-ixOOmin1
                  if (i==0) then
                    w(ixOOmin1+i, ixOmin2:ixOmax2, iw) = &
                         4 * tmp(ixOmin1, ixOmin2:ixOmax2) - &
                         3 * w(ix, ixOmin2:ixOmax2, iw)
                  elseif (i==1)then  
                    w(ixOOmin1+i, ixOmin2:ixOmax2, iw) = &
                         2 * tmp(ixOmin1, ixOmin2:ixOmax2) - &
                         w(ix, ixOmin2:ixOmax2, iw)
                   else 
                    w(ixOOmin1+i, ixOmin2:ixOmax2, iw) = &
                         tmp(ixOmin1, ixOmin2:ixOmax2)
                   endif 
                  
                end do
              else
                do i = 0, ixOOmax1-ixOOmin1
                  if (i==ixOOmax1-ixOOmin1) then
                    w(ixOOmin1+i, ixOmin2:ixOmax2, iw) = &
                         4 * tmp(ixOmin1, ixOmin2:ixOmax2) - &
                         3 * w(ix, ixOmin2:ixOmax2, iw)
                  elseif (i==ixOOmax1-ixOOmin1-1)then  
                    w(ixOOmin1+i, ixOmin2:ixOmax2, iw) = &
                         2 * tmp(ixOmin1, ixOmin2:ixOmax2) - &
                         w(ix, ixOmin2:ixOmax2, iw)
                  else 
                    w(ixOOmin1+i, ixOmin2:ixOmax2, iw) = &
                         tmp(ixOmin1, ixOmin2:ixOmax2)
                  endif 
                  
                end do
              endif

            endif
          else
            do i = ixOOmin1, ixOOmax1
               w(i, ixOmin2:ixOmax2, iw) = &
                   tmp(ixOmin1, ixOmin2:ixOmax2)
            end do

          endif
       end do
       if(interp_phy_first_row) then
            if(boundary_data_primitive) then
              call phys_to_conserved(ixI^L,ix,ixOmin2,ix,ixOmax2,w,x)
            endif   
       endif   
    case (3, 4)
       if(interp_phy_first_row) then
            if (iB == 3) then
               ix = ixOmax2+1
            else
               ix = ixOmin2-1
            end if
            if(boundary_data_primitive) then
              call phys_to_primitive(ixI^L,ixOmin1,ix,ixOmax1,ix,w,x)
            endif   
       endif   
       if(bc_phy_first_row)then 
         if (iB == 3) then
           ixOOmax2=ixOOmax2+1 
         else
           ixOOmin2=ixOOmin2-1 
         endif
         if(boundary_data_primitive) then
           if (iB == 3) then
             ixOOmin2=ixOOmax2 
           else
             ixOOmax2=ixOOmin2
           endif
           call phys_to_primitive(ixI^L,ixOO^L,w,x)
           if (iB == 3) then
             ixOOmin2=ixOmin2 ! to be used later in to_conserved
           else
             ixOOmax2=ixOmax2
           endif
         endif
       endif
       do iw = 1, nwfluxbc
          n_bc = bc_data_ix(iw, iB)
          if (n_bc == -1) cycle

          tmp(ixOmin1:ixOmax1, ixOmin2) = bc_data_get_2d(n_bc, &
                  x(ixOmin1:ixOmax1, ixOmin2, 1), qt)

          if(interp_phy_first_row) then
  
            ! Approximate boundary value by linear interpolation to first ghost
            ! cell, rest of ghost cells contains the same value
            if(interp_phy_first_row_same) then
            do i = 0, ixOOmax2-ixOOmin2
               w(ixOmin1:ixOmax1, ixOOmin2+i, iw) = &
                    2 * tmp(ixOmin1:ixOmax1, ixOmin2) - &
                    w(ixOmin1:ixOmax1, ix, iw)
            end do
            else
              if (iB==3) then
                do i = 0, ixOOmax2-ixOOmin2
                  if (i==0) then
                    w(ixOmin1:ixOmax1, ixOOmin2+i, iw) = &
                    4 * tmp(ixOmin1:ixOmax1, ixOmin2) - &
                    3 * w(ixOmin1:ixOmax1, ix, iw)
                  elseif (i==1)then  
                    w(ixOmin1:ixOmax1, ixOOmin2+i, iw) = &
                    2 * tmp(ixOmin1:ixOmax1, ixOmin2) - &
                    w(ixOmin1:ixOmax1, ix, iw)
                   else 
                    w(ixOmin1:ixOmax1, ixOOmin2+i, iw) = &
                    tmp(ixOmin1:ixOmax1, ixOmin2) 
                   endif 
                  
                end do
              else
                do i = 0, ixOOmax2-ixOOmin2
                  if (i==ixOOmax2-ixOOmin2) then
                    w(ixOmin1:ixOmax1, ixOOmin2+i, iw) = &
                    4 * tmp(ixOmin1:ixOmax1, ixOmin2) - &
                    3 * w(ixOmin1:ixOmax1, ix, iw)
                  elseif (i==ixOOmax2-ixOOmin2-1)then  
                    w(ixOmin1:ixOmax1, ixOOmin2+i, iw) = &
                    2 * tmp(ixOmin1:ixOmax1, ixOmin2) - &
                    w(ixOmin1:ixOmax1, ix, iw)
                   else 
                    w(ixOmin1:ixOmax1, ixOOmin2+i, iw) = &
                    tmp(ixOmin1:ixOmax1, ixOmin2) 
                   endif 
                  
                end do
              endif

            endif
          else
            do i = ixOOmin2, ixOOmax2
             w(ixOmin1:ixOmax1, i, iw) = &
                  tmp(ixOmin1:ixOmax1, ixOmin2)
            end do

          endif
       end do
       if(interp_phy_first_row) then
            if(boundary_data_primitive) then
              call phys_to_conserved(ixI^L,ixOmin1,ix,ixOmax1,ix,w,x)
            endif   
       endif   
    case default
       call mpistop("bc_data_set: unknown iB")
    end select
    }

    {^IFTHREED
    select case (iB)
    case (1, 2)
       if(interp_phy_first_row) then
            if (iB == 1) then
               ix = ixOmax1+1
            else
               ix = ixOmin1-1
            end if
            if(boundary_data_primitive) then
              call phys_to_primitive(ixI^L,ix,ixOmin2,ixOmin3,ix,ixOmax2,ixOmax3,w,x)
            endif   
       endif 
       if(bc_phy_first_row)then 
         if (iB == 1) then
           ixOOmax1=ixOOmax1+1 
         else
           ixOOmin1=ixOOmin1-1 
         endif
         if(boundary_data_primitive) then
           if (iB == 1) then
             ixOOmin1=ixOOmax1 
           else
             ixOOmax1=ixOOmin1
           endif
           call phys_to_primitive(ixI^L,ixOO^L,w,x)
           if (iB == 1) then
             ixOOmin1=ixOmin1 ! to be used later in to_conserved
           else
             ixOOmax1=ixOmax1
           endif
         endif
       endif
       do iw = 1, nwfluxbc
          n_bc = bc_data_ix(iw, iB)
          if (n_bc == -1) cycle

          tmp(ixOmin1, ixOmin2:ixOmax2, ixOmin3:ixOmax3) = bc_data_get_3d(n_bc, &
               x(ixOmin1, ixOmin2:ixOmax2, ixOmin3:ixOmax3, 2), &
               x(ixOmin1, ixOmin2:ixOmax2, ixOmin3:ixOmax3, 3), qt)

          !print*,"I ",  mype, it,n_bc, minval(tmp(ixOmin1, ixOmin2:ixOmax2, ixOmin3:ixOmax3)),&
          !    maxval(tmp(ixOmin1, ixOmin2:ixOmax2, ixOmin3:ixOmax3))
          if(interp_phy_first_row) then
  
            ! Approximate boundary value by linear interpolation to first ghost
            ! cell, rest of ghost cells contains the same value
            if(interp_phy_first_row_same) then
              do i = 0, ixOOmax1-ixOOmin1
                 w(ixOOmin1+i, ixOmin2:ixOmax2, ixOmin3:ixOmax3, iw) = &
                      2 * tmp(ixOmin1, ixOmin2:ixOmax2, ixOmin3:ixOmax3) - &
                      w(ix, ixOmin2:ixOmax2, ixOmin3:ixOmax3, iw)
              end do
            else
              if (iB==1) then
                do i = 0, ixOOmax1-ixOOmin1
                  if (i==0) then
                    w(ixOOmin1+i, ixOmin2:ixOmax2, ixOmin3:ixOmax3, iw) = &
                         4 * tmp(ixOmin1, ixOmin2:ixOmax2, ixOmin3:ixOmax3) - &
                         3*w(ix, ixOmin2:ixOmax2, ixOmin3:ixOmax3, iw)
                  elseif (i==1)then  
                    w(ixOOmin1+i, ixOmin2:ixOmax2, ixOmin3:ixOmax3, iw) = &
                         2 * tmp(ixOmin1, ixOmin2:ixOmax2, ixOmin3:ixOmax3) - &
                         w(ix, ixOmin2:ixOmax2, ixOmin3:ixOmax3, iw)
                   else 
                    w(ixOOmin1+i, ixOmin2:ixOmax2, ixOmin3:ixOmax3, iw) = &
                         tmp(ixOmin1, ixOmin2:ixOmax2, ixOmin3:ixOmax3) 
                   endif 
                  
                end do
              else
                do i = 0, ixOOmax1-ixOOmin1
                  if (i==ixOOmax1-ixOOmin1) then
                    w(ixOOmin1+i, ixOmin2:ixOmax2, ixOmin3:ixOmax3, iw) = &
                         4 * tmp(ixOmin1, ixOmin2:ixOmax2, ixOmin3:ixOmax3) - &
                         3*w(ix, ixOmin2:ixOmax2, ixOmin3:ixOmax3, iw)
                  elseif (i==ixOOmax1-ixOOmin1-1)then  
                    w(ixOOmin1+i, ixOmin2:ixOmax2, ixOmin3:ixOmax3, iw) = &
                         2 * tmp(ixOmin1, ixOmin2:ixOmax2, ixOmin3:ixOmax3) - &
                         w(ix, ixOmin2:ixOmax2, ixOmin3:ixOmax3, iw)
                   else 
                    w(ixOOmin1+i, ixOmin2:ixOmax2, ixOmin3:ixOmax3, iw) = &
                         tmp(ixOmin1, ixOmin2:ixOmax2, ixOmin3:ixOmax3) 
                   endif 
                  
                end do
              endif

            endif
          else
            do i = ixOOmin1,ixOOmax1
               w(i, ixOmin2:ixOmax2, ixOmin3:ixOmax3, iw) = &
                     tmp(ixOmin1, ixOmin2:ixOmax2, ixOmin3:ixOmax3) 
            end do
          end if

       end do
       if(interp_phy_first_row) then
            if(boundary_data_primitive) then
              call phys_to_conserved(ixI^L,ix,ixOmin2,ixOmin3,ix,ixOmax2,ixOmax3,w,x)
            endif  
       endif 
    case (3, 4)
        if(interp_phy_first_row) then
            if (iB == 3) then
               ix = ixOmax2+1
            else
               ix = ixOmin2-1
            end if
            if(boundary_data_primitive) then
              call phys_to_primitive(ixI^L,ixOmin1,ix,ixOmin3,ixOmax1,ix,ixOmax3,w,x)
            endif   
        endif 
        if(bc_phy_first_row)then 
          if (iB == 3) then
            ixOOmax2=ixOOmax2+1 
          else
            ixOOmin2=ixOOmin2-1 
          endif
          if(boundary_data_primitive) then
            if (iB == 3) then
              ixOOmin2=ixOOmax2 
            else
              ixOOmax2=ixOOmin2
            endif
            call phys_to_primitive(ixI^L,ixOO^L,w,x)
            if (iB == 3) then
              ixOOmin2=ixOmin2 ! to be used later in to_conserved
            else
              ixOOmax2=ixOmax2
            endif
          endif
       endif
       do iw = 1, nwfluxbc
          n_bc = bc_data_ix(iw, iB)
          if (n_bc == -1) cycle

          tmp(ixOmin1:ixOmax1, ixOmin2, ixOmin3:ixOmax3) = bc_data_get_3d(n_bc, &
                  x(ixOmin1:ixOmax1, ixOmin2, ixOmin3:ixOmax3, 1), &
                  x(ixOmin1:ixOmax1, ixOmin2, ixOmin3:ixOmax3, 3), qt)

          if(interp_phy_first_row) then
  
            if(interp_phy_first_row_same) then
            ! Approximate boundary value by linear interpolation to first ghost
            ! cell, rest of ghost cells contains the same value
              do i = 0, ixOOmax2-ixOOmin2
                 w(ixOmin1:ixOmax1, ixOOmin2+i, ixOmin3:ixOmax3, iw) = &
                      2 * tmp(ixOmin1:ixOmax1, ixOmin2, ixOmin3:ixOmax3) - &
                      w(ixOmin1:ixOmax1, ix, ixOmin3:ixOmax3, iw)
              end do
            else
              if (iB==3) then
                do i = 0, ixOOmax2-ixOOmin2
                  if (i==0) then
                     w(ixOmin1:ixOmax1, ixOOmin2+i, ixOmin3:ixOmax3, iw) = &
                          4 * tmp(ixOmin1:ixOmax1, ixOmin2, ixOmin3:ixOmax3) - &
                          3* w(ixOmin1:ixOmax1, ix, ixOmin3:ixOmax3, iw)
                  elseif (i==1)then  
                    w(ixOmin1:ixOmax1, ixOOmin2+i, ixOmin3:ixOmax3, iw) = &
                          2 * tmp(ixOmin1:ixOmax1, ixOmin2, ixOmin3:ixOmax3) - &
                          w(ixOmin1:ixOmax1, ix, ixOmin3:ixOmax3, iw)
                   else 
                     w(ixOmin1:ixOmax1, ixOOmin2+i, ixOmin3:ixOmax3, iw) = &
                          tmp(ixOmin1:ixOmax1, ixOmin2, ixOmin3:ixOmax3) 
                   endif 
                  
                end do
              else
                do i = 0, ixOOmax2-ixOOmin2
                  if (i==ixOOmax2-ixOOmin2) then
                     w(ixOmin1:ixOmax1, ixOOmin2+i, ixOmin3:ixOmax3, iw) = &
                          4 * tmp(ixOmin1:ixOmax1, ixOmin2, ixOmin3:ixOmax3) - &
                          3* w(ixOmin1:ixOmax1, ix, ixOmin3:ixOmax3, iw)
                  elseif (i==ixOOmax2-ixOOmin2-1)then  
                    w(ixOmin1:ixOmax1, ixOOmin2+i, ixOmin3:ixOmax3, iw) = &
                          2 * tmp(ixOmin1:ixOmax1, ixOmin2, ixOmin3:ixOmax3) - &
                          w(ixOmin1:ixOmax1, ix, ixOmin3:ixOmax3, iw)
                   else 
                    w(ixOmin1:ixOmax1, ixOOmin2+i, ixOmin3:ixOmax3, iw) = &
                          tmp(ixOmin1:ixOmax1, ixOmin2, ixOmin3:ixOmax3) 
                   endif 
                  
                end do
              endif

            endif
          else
            do i = ixOOmin2,ixOOmax2
               w(ixOmin1:ixOmax1, i, ixOmin3:ixOmax3, iw) = &
                    tmp(ixOmin1:ixOmax1, ixOmin2, ixOmin3:ixOmax3)
            end do
          end if
       end do
       if(interp_phy_first_row) then
            if(boundary_data_primitive) then
              call phys_to_conserved(ixI^L,ixOmin1,ix,ixOmin3,ixOmax1,ix,ixOmax3,w,x)
            endif   
       endif 
    case (5, 6)
       if(interp_phy_first_row) then
            if (iB == 5) then
               ix = ixOmax3+1
            else
               ix = ixOmin3-1
            end if
            if(boundary_data_primitive) then
              call phys_to_primitive(ixI^L,ixOmin1,ixOmax1,ixOmin2,ixOmax2,ix,ix,w,x)
            endif   
       endif 
       if(bc_phy_first_row)then 
         if (iB == 5) then
           ixOOmax3=ixOOmax3+1 
         else
           ixOOmin3=ixOOmin3-1 
         endif
         if(boundary_data_primitive) then
           if (iB == 5) then
             ixOOmin3=ixOOmax3 
           else
             ixOOmax3=ixOOmin3
           endif
           call phys_to_primitive(ixI^L,ixOO^L,w,x)
           if (iB == 5) then
             ixOOmin3=ixOmin3 ! to be used later in to_conserved
           else
             ixOOmax3=ixOmax3
           endif
         endif
       endif
       do iw = 1, nwfluxbc
          n_bc = bc_data_ix(iw, iB)
          if (n_bc == -1) cycle

          tmp(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOmin3) = bc_data_get_3d(n_bc, &
                  x(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOmin3, 1), &
                  x(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOmin3, 2), qt)

          if(interp_phy_first_row) then
  
            if(interp_phy_first_row_same) then
              ! Approximate boundary value by linear interpolation to first ghost
              ! cell, rest of ghost cells contains the same value
              do i = 0, ixOOmax3-ixOOmin3
                 w(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOOmin3+i, iw) = &
                      2 * tmp(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOmin3) - &
                      w(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ix, iw)
              end do
            else
              if (iB==5) then
                do i = 0, ixOOmax3-ixOOmin3
                  if (i==0) then
                    w(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOOmin3+i, iw) = &
                        4 * tmp(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOmin3) - &
                        3 * w(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ix, iw)
                  elseif (i==1)then  
                    w(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOOmin3+i, iw) = &
                        2 * tmp(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOmin3) - &
                        w(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ix, iw)
                   else 
                    w(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOOmin3+i, iw) = &
                        tmp(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOmin3)
                   endif 
                  
                end do
              else
                do i = 0, ixOOmax3-ixOOmin3
                  if (i==ixOOmax3-ixOOmin3) then
                    w(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOOmin3+i, iw) = &
                        4 * tmp(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOmin3) - &
                        3 * w(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ix, iw)
                  elseif (i==ixOOmax3-ixOOmin3-1)then  
                    w(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOOmin3+i, iw) = &
                        2 * tmp(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOmin3) - &
                        w(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ix, iw)
                   else 
                    w(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOOmin3+i, iw) = &
                        tmp(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOmin3)
                   endif 
                  
                end do
              endif

            endif
          else
            do i = ixOOmin3,ixOOmax3
               w(ixOmin1:ixOmax1, ixOmin2:ixOmax2, i, iw) = &
                    tmp(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOmin3)
            end do
          end if
       end do
       if(interp_phy_first_row) then
            if(boundary_data_primitive) then
              call phys_to_conserved(ixI^L,ixOmin1,ixOmax1,ixOmin2,ixOmax2,ix,ix,w,x)
            endif   
       endif 
    case default
       call mpistop("bc_data_set: unknown iB")
    end select
    }

    if(boundary_data_primitive) then
      call phys_to_conserved(ixI^L,ixOO^L,w,x)
    endif

  end subroutine bc_data_set

  subroutine read_vtk_structured_points(fname, bc)
    character(len=*), intent(in)  :: fname
    type(bc_data_t), intent(out)  :: bc
    double precision, allocatable :: tmp_data(:, :, :, :)
    integer, parameter            :: max_variables = 10
    character(len=40)             :: tmp_names(max_variables)
    character(len=256)            :: line
    character(len=40)             :: word, typename
    integer, parameter            :: my_unit = 123
    integer                       :: n, n_points_total
    integer                       :: n_components

    open(my_unit, file=trim(fname), status="old", action="read")

    ! Header, e.g. # vtk DataFile Version 2.0
    read(my_unit, "(A)") line

    ! Dataset name
    read(my_unit, "(A)") line

    ! ASCII / BINARY
    read(my_unit, "(A)") line

    if (line /= "ASCII") then
       print *, "line: ", trim(line)
       error stop "read_vtk: not an ASCII file"
    end if

    ! DATASET STRUCTURED_POINTS
    read(my_unit, "(A)") line

    if (line /= "DATASET STRUCTURED_POINTS") then
       print *, "line: ", trim(line)
       error stop "read_vtk must have: DATASET STRUCTURED_POINTS"
    end if

    ! DIMENSIONS NX NY NZ
    read(my_unit, "(A)") line
    read(line, *) word, bc%n_points

    if (word /= "DIMENSIONS") then
       print *, "line: ", trim(line)
       error stop "read_vtk expects: DIMENSIONS"
    end if

    ! SPACING DX DY DZ
    read(my_unit, *) word, bc%dx

    if (word /= "SPACING") then
       print *, "line: ", trim(line)
       error stop "read_vtk expects: SPACING"
    end if

    ! ORIGIN 0 0 0
    read(my_unit, *) word, bc%origin
    if (word /= "ORIGIN") then
       print *, "line: ", trim(line)
       error stop "read_vtk expects: ORIGIN"
    end if

    ! POINT_DATA NPOINTS
    read(my_unit, *) word, n_points_total

    if (word /= "POINT_DATA") then
       print *, "line: ", trim(line)
       error stop "read_vtk expects: POINT_DATA n_points"
    end if

    if (n_points_total /= product(bc%n_points)) &
         error stop "read_vtk: n_points not equal to NX*NY*NZ"

    allocate(tmp_data(bc%n_points(1), bc%n_points(2), bc%n_points(3), &
         max_variables))

    ! Read all scalar variables
    do n = 1, max_variables

       ! SCALARS name type ncomponents
       read(my_unit, *, end=900) word, tmp_names(n), typename, n_components

       if (word /= "SCALARS") then
          print *, "line: ", trim(line)
          error stop "read_vtk expects: SCALARS name type ncomponents"
       end if

       if (n_components /= 1) error stop "read_vtk: ncomponents should be 1"

       ! LOOKUP_TABLE default
       read(my_unit, *) word, typename

       if (word /= "LOOKUP_TABLE") then
          print *, "line: ", trim(line)
          error stop "read_vtk expects: LOOKUP_TABLE name"
       end if

       ! Read list of data values
       read(my_unit, *) tmp_data(:, :, :, n)
    end do

    ! Done reading variables
900 continue

    close(my_unit)

    if (n == max_variables + 1) &
         error stop "read_vtk: increase max_variables"

    ! Loop index is one higher than number of variables
    bc%n_variables = n-1
    bc%values      = tmp_data(:, :, :, 1:n-1)
    bc%names       = tmp_names(1:n-1)
  end subroutine read_vtk_structured_points

end module mod_bc_data
