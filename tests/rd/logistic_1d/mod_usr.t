module mod_usr

  use mod_rd
  implicit none
    
  integer, parameter :: dp = kind(0.0d0)

contains

  subroutine usr_init()
     integer :: i

     call set_coordinate_system('Cartesian')
     call rd_activate()

     usr_init_one_grid => lg_init     
   
     if (iprob == 2) then
        usr_special_bc    => lg_bound
        do i = 1, 2*ndim
           rd_mg_bc(1, i)%boundary_cond => lg_bound_mg
        end do
     end if
  end subroutine usr_init

  subroutine lg_init(ixG^L,ix^L,w,x)
     integer, intent(in)             :: ixG^L, ix^L
     double precision, intent(in)    :: x(ixG^S,1:ndim)
     double precision, intent(inout) :: w(ixG^S,1:nw)
     double precision                :: urand(ix^S)
     double precision                :: ssu, ssv, ssw

     select case (iprob)
     case (1)
        ! From https://people.maths.ox.ac.uk/trefethen/pdectb/fisher2.pdf
        w(ix^S, u_) = 0.0d0
        where (x(ix^S, 1) < -10.0d0)
           w(ix^S, u_) = 1.0d0
        end where
        where (10.0d0 < x(ix^S, 1) .and. x(ix^S, 1) < 20.0d0)
           w(ix^S, u_) = 1.0d0 / 4.0d0
        end where
     case (2)
        ! Analytical solution of Fisher equation
        w(ix^S, u_) = lg_solution(x(ix^S, 1), 0.0_dp)
    
     case default
        call mpistop("Unknown iprob")
     end select
  end subroutine lg_init

  !> Boundary conditions for the analytical solution
  subroutine lg_bound(qt,ixG^L,ixO^L,iB,w,x)
     integer, intent(in)             :: ixG^L, ixO^L, iB
     double precision, intent(in)    :: qt, x(ixG^S,1:ndim)
     double precision, intent(inout) :: w(ixG^S,1:nw)

     w(ixO^S, u_) = lg_solution(x(ixO^S,1), qt)
  end subroutine lg_bound

  !> Boundary conditions for the analytical solution when using IMEX scheme
  subroutine lg_bound_mg(box, nc, iv, nb, bc_type, bc)
     use mod_multigrid_coupling
     type(mg_box_t), intent(in)    :: box
     integer, intent(in)           :: nc
     integer, intent(in)           :: iv      !< Index of variable
     integer, intent(in)           :: nb      !< Direction
     integer, intent(out)          :: bc_type !< Type of b.c.
     double precision, intent(out) :: bc(1)
     double precision              :: rr(2)
     integer                       :: i, ixs(ndim-1), nb_dim, n_bc

     ! Type of boundary condition
     bc_type = mg_bc_dirichlet

     ! Get the coordinates at the cell faces at the boundary
     call mg_get_face_coords(box, nb, nc, rr)

     bc = lg_solution(rr(1), global_time)
  end subroutine lg_bound_mg
  
  !> An analytical solution of the Fisher equation (for a certain wave speed)
  elemental function lg_solution(x,t) result(val)
     real(dp), intent(in) :: x, t
     real(dp)             :: val

     val = ( 1.0_dp / (1.0_dp + dexp((x+5.0_dp/dsqrt(6.0_dp)*t)/dsqrt(6.0_dp))) )**2  
  end function lg_solution

end module mod_usr
