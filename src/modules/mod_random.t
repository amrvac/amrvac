!> Module for pseudo random number generation. The internal pseudo random
!> generator is the xoroshiro128plus method.
module mod_random

#include "amrvac.h"

  implicit none
  private

  ! A 64 bit floating point type
  integer, parameter :: dp = kind(0.0d0)

  ! A 32 bit integer type
  integer, parameter :: i4 = selected_int_kind(9)

  ! A 64 bit integer type
  integer, parameter :: i8 = selected_int_kind(18)

  !> Random number generator type, which contains the state
  type rng_t
     !> The rng state (always use your own seed)
     integer(i8), private       :: s(2) = [123456789_i8, 987654321_i8]
   contains
     procedure, non_overridable :: set_seed    ! Seed the generator
     procedure, non_overridable :: jump        ! Jump function (see below)
     procedure, non_overridable :: int_4       ! 4-byte random integer
     procedure, non_overridable :: int_8       ! 8-byte random integer
     procedure, non_overridable :: unif_01     ! Uniform (0,1] real
     procedure, non_overridable :: unif_01_vec ! Uniform (0,1] real vector
     procedure, non_overridable :: normal      ! One normal(0,1) number
     procedure, non_overridable :: two_normals ! Two normal(0,1) samples
     procedure, non_overridable :: poisson     ! Sample from Poisson-dist.
     procedure, non_overridable :: circle      ! Sample on a circle
     procedure, non_overridable :: sphere      ! Sample on a sphere
     procedure, non_overridable :: next        ! Internal method
  end type rng_t

  !> Parallel random number generator type
  type prng_t
     type(rng_t), allocatable :: rngs(:)
   contains
     procedure, non_overridable :: init_parallel
  end type prng_t

  public :: rng_t
  public :: prng_t
  public :: Poisson_disk_sampling

contains


#if defined(__NVCOMPILER) ||  (defined(USE_INTRINSIC_SHIFT) && USE_INTRINSIC_SHIFT==0)

  !> added for nvidia compilers
  pure function shiftl(val, shift) result(res_value)
    integer(i8), intent(in) :: val
    integer, intent(in) :: shift
    integer(i8) :: res_value
    integer(i8) :: bit_mask1, bit_mask2

    ! cannot be initialized with b values in gnu, cannot have the big decimal numbers in nvidia 
    bit_mask1=huge(val)
    bit_mask2=-bit_mask1-1


    if(val<0) then
      res_value = ior(lshift(iand(val, bit_mask1), shift),bit_mask2)
    else
      res_value = lshift(val, shift)
    endif  

  end function shiftl

  pure function shiftr(val, shift) result(res_value)
    integer(i8), intent(in) :: val
    integer, intent(in) :: shift
    integer(i8) :: res_value
    integer(i8) :: bit_mask1, bit_mask2

    ! cannot be initialized with b values in gnu, cannot have the big decimal numbers in nvidia 
    bit_mask1=huge(val)
    bit_mask2=-bit_mask1-1

    if(val<0) then
      res_value = ior(rshift(iand(val, bit_mask1), shift),bit_mask2)
    else
      res_value = rshift(val, shift)
    endif  
    
  end function shiftr

#endif

  !> Initialize a collection of rng's for parallel use
  subroutine init_parallel(self, n_proc, rng)
    class(prng_t), intent(inout) :: self
    type(rng_t), intent(inout)   :: rng
    integer, intent(in)          :: n_proc
    integer                      :: n

    allocate(self%rngs(n_proc))
    self%rngs(1) = rng

    do n = 2, n_proc
       self%rngs(n) = self%rngs(n-1)
       call self%rngs(n)%jump()
    end do
  end subroutine init_parallel

  !> Set a seed for the rng
  subroutine set_seed(self, the_seed)
    class(rng_t), intent(inout) :: self
    integer(i8), intent(in)     :: the_seed(2)

    self%s = the_seed

    ! Simulate calls to next() to improve randomness of first number
    call self%jump()
  end subroutine set_seed

  ! This is the jump function for the generator. It is equivalent
  ! to 2^64 calls to next(); it can be used to generate 2^64
  ! non-overlapping subsequences for parallel computations.
  subroutine jump(self)
    class(rng_t), intent(inout) :: self
    integer                     :: i, b
    integer(i8)                 :: global_time(2), dummy

    ! The signed equivalent of the unsigned constants
    integer(i8), parameter      :: jmp_c(2) = &
         (/-4707382666127344949_i8, -2852180941702784734_i8/)

    global_time = 0
    do i = 1, 2
       do b = 0, 63
          if (iand(jmp_c(i), shiftl(1_i8, b)) /= 0) then
             global_time = ieor(global_time, self%s)
          end if
          dummy = self%next()
       end do
    end do

    self%s = global_time
  end subroutine jump

  !> Return 4-byte integer
  integer(i4) function int_4(self)
    class(rng_t), intent(inout) :: self
    int_4 = int(self%next(), i4)
  end function int_4

  !> Return 8-byte integer
  integer(i8) function int_8(self)
    class(rng_t), intent(inout) :: self
    int_8 = self%next()
  end function int_8

  !> Get a uniform [0,1) random real (double precision)
  real(dp) function unif_01(self)
    class(rng_t), intent(inout) :: self
    integer(i8)                 :: x
    real(dp)                    :: tmp

    x   = self%next()
    x   = ior(shiftl(1023_i8, 52), shiftr(x, 12))
    unif_01 = transfer(x, tmp) - 1.0_dp
  end function unif_01

  !> Fill array with uniform random numbers
  subroutine unif_01_vec(self, rr)
    class(rng_t), intent(inout) :: self
    real(dp), intent(out)       :: rr(:)
    integer                     :: i

    do i = 1, size(rr)
      rr(i) = self%unif_01()
    end do
  end subroutine unif_01_vec

  !> Return two normal random variates with mean 0 and variance 1.
  !> http://en.wikipedia.org/wiki/Marsaglia_polar_method
  function two_normals(self) result(rands)
    class(rng_t), intent(inout) :: self
    real(dp)                    :: rands(2), sum_sq

    do
       rands(1) = 2 * self%unif_01() - 1
       rands(2) = 2 * self%unif_01() - 1
       sum_sq = sum(rands**2)
       if (sum_sq < 1.0_dp .and. sum_sq > 0.0_dp) exit
    end do
    rands = rands * sqrt(-2 * log(sum_sq) / sum_sq)
  end function two_normals

  !> Single normal random number
  real(dp) function normal(self)
    class(rng_t), intent(inout) :: self
    real(dp)                    :: rands(2)

    rands  = self%two_normals()
    normal = rands(1)
  end function normal

  !> Return Poisson random variate with rate lambda. Works well for lambda < 30
  !> or so. For lambda >> 1 it can produce wrong results due to roundoff error.
  function poisson(self, lambda) result(rr)
    class(rng_t), intent(inout) :: self
    real(dp), intent(in)        :: lambda
    integer(i4)                 :: rr
    real(dp)                    :: expl, p

    expl = exp(-lambda)
    rr   = 0
    p    = self%unif_01()

    do while (p > expl)
       rr = rr + 1
       p = p * self%unif_01()
    end do
  end function poisson

  !> Sample point on a circle with given radius
  function circle(self, radius) result(xy)
    class(rng_t), intent(inout) :: self
    real(dp), intent(in)        :: radius
    real(dp)                    :: rands(2), xy(2)
    real(dp)                    :: sum_sq

    ! Method for uniform sampling on circle
    do
       rands(1) = 2 * self%unif_01() - 1
       rands(2) = 2 * self%unif_01() - 1
       sum_sq   = sum(rands**2)
       if (sum_sq <= 1) exit
    end do

    xy(1) = (rands(1)**2 - rands(2)**2) / sum_sq
    xy(2) = 2 * rands(1) * rands(2) / sum_sq
    xy    = xy * radius
  end function circle

  !> Sample point on a sphere with given radius
  function sphere(self, radius) result(xyz)
    class(rng_t), intent(inout) :: self
    real(dp), intent(in)        :: radius
    real(dp)                    :: rands(2), xyz(3)
    real(dp)                    :: sum_sq, tmp_sqrt

    ! Marsaglia method for uniform sampling on sphere
    do
       rands(1) = 2 * self%unif_01() - 1
       rands(2) = 2 * self%unif_01() - 1
       sum_sq   = sum(rands**2)
       if (sum_sq <= 1) exit
    end do

    tmp_sqrt = sqrt(1 - sum_sq)
    xyz(1:2) = 2 * rands(1:2) * tmp_sqrt
    xyz(3)   = 1 - 2 * sum_sq
    xyz      = xyz * radius
  end function sphere

  !> Interal routine: get the next value (returned as 64 bit signed integer)
  function next(self) result(res)
    class(rng_t), intent(inout) :: self
    integer(i8)                 :: res
    integer(i8)                 :: global_time(2)

    global_time         = self%s
    res       = global_time(1) + global_time(2)
    global_time(2)      = ieor(global_time(1), global_time(2))
    self%s(1) = ieor(ieor(rotl(global_time(1), 55), global_time(2)), shiftl(global_time(2), 14))
    self%s(2) = rotl(global_time(2), 36)
  end function next

  !> Helper function for next()
  pure function rotl(x, k) result(res)
    integer(i8), intent(in) :: x
    integer, intent(in)     :: k
    integer(i8)             :: res

    res = ior(shiftl(x, k), shiftr(x, 64 - k))
  end function rotl

  subroutine Poisson_disk_sampling(r_min,bmin1,bmin2,bmax1,bmax2,nmax,points_store,n_setpoints)
    use mod_constants
    ! a bigger number than expected total number of sampling points
    integer, intent(in) :: nmax
    double precision, intent(in)  :: r_min,bmin1,bmin2,bmax1,bmax2
    double precision, intent(out) :: points_store(nmax,2)
    integer, intent(out) :: n_setpoints

    double precision :: gsize,x1len,x2len,R2,threeR2,r_rand,a_rand,x_rand,y_rand
    double precision :: unif_random_number(2),candidate(2),distance,random01
    type(rng_t) :: rng
    ! the id of a point in a background grid, in which each cell can only contain at most one point
    integer, allocatable :: id_grid(:,:)
    integer :: nx1,nx2,ix1,ix2,n_active,id,j,k_candi
    integer :: id_active_points(nmax),point_index_in_grid(2,nmax)
    logical :: points_active(nmax)
    logical :: candidate_found

    k_candi=30
    id_active_points=0
    R2=r_min**2
    threeR2=3.d0*R2
    points_active=.false.
    ! cell size of the background grid
    gsize=r_min/sqrt(2.d0)
    ! size of region to sample
    x1len=bmax1-bmin1
    x2len=bmax2-bmin2
    nx1=ceiling(x1len/gsize)
    nx2=ceiling(x2len/gsize)
    allocate(id_grid(nx1,nx2))
    id_grid=0
    ! set the first point
    call rng%unif_01_vec(unif_random_number)
    points_store(1,1)=unif_random_number(1)*x1len
    points_store(1,2)=unif_random_number(2)*x2len
    ix1=ceiling(points_store(1,1)/gsize)
    ix2=ceiling(points_store(1,2)/gsize)
    id_grid(ix1,ix2)=1
    point_index_in_grid(1,1)=ix1
    point_index_in_grid(2,1)=ix2
    points_active(1)=.true.
    n_setpoints=1

    do while(any(points_active))
      ! generate active point list
      n_active=0
      do id=1,nmax
        if(points_active(id)) then
          n_active=n_active+1
          id_active_points(n_active)=id
        end if
      end do
      ! randomly select one point from the active point list
      id=id_active_points(ceiling(rng%unif_01()*dble(n_active)))
      do j=1,k_candi
        candidate_found=.true.
        r_rand=sqrt(rng%unif_01()*threeR2+R2)
        a_rand=2.d0*dpi*rng%unif_01()
        x_rand=points_store(id,1)+r_rand*cos(a_rand)
        y_rand=points_store(id,2)+r_rand*sin(a_rand)
        if(x_rand<0.d0 .or. x_rand> x1len .or. y_rand<0.d0 .or. y_rand > x2len) cycle
        do ix2=max(1,point_index_in_grid(2,id)-3),min(nx2,point_index_in_grid(2,id)+3)
          do ix1=max(1,point_index_in_grid(1,id)-3),min(nx1,point_index_in_grid(1,id)+3)
            if(id_grid(ix1,ix2)==id) cycle
            if(id_grid(ix1,ix2)/=0) then
              distance=sqrt((x_rand-points_store(id_grid(ix1,ix2),1))**2+(y_rand-points_store(id_grid(ix1,ix2),2))**2)
              if(distance<r_min) then
                candidate_found=.false.
                go to 12
              end if
            end if
          end do
        end do
        if(candidate_found) then
          n_setpoints=n_setpoints+1
          points_store(n_setpoints,1)=x_rand
          points_store(n_setpoints,2)=y_rand
          ix1=ceiling(points_store(n_setpoints,1)/gsize)
          ix2=ceiling(points_store(n_setpoints,2)/gsize)
          point_index_in_grid(1,n_setpoints)=ix1
          point_index_in_grid(2,n_setpoints)=ix2
          id_grid(ix1,ix2)=n_setpoints
          points_active(n_setpoints)=.true.
          exit
        end if
12      continue
      end do
      if(.not.candidate_found) then
        points_active(id)=.false.
      end if
    end do
    points_store(1:n_setpoints,1)=points_store(1:n_setpoints,1)+bmin1
    points_store(1:n_setpoints,2)=points_store(1:n_setpoints,2)+bmin2
    deallocate(id_grid)

  end subroutine Poisson_disk_sampling

end module mod_random
