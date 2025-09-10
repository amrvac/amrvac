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

contains


#if defined(__NVCOMPILER) ||  (defined(USE_INTRINSIC_SHIFT) && USE_INTRINSIC_SHIFT==0)

  !> added for nvidia compilers
   pure function shiftl(val, shift) result(res_value)
    integer(i8), intent(in) :: val
    integer, intent(in) :: shift
    integer(i8) :: res_value
!    integer(i8),parameter :: bit_mask1=9223372036854775807
!    !b'0111111111111111111111111111111111111111111111111111111111111111'
!    integer(i8),parameter :: bit_mask2=-9223372036854775808   
!    !b'1000000000000000000000000000000000000000000000000000000000000000'
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
    !    integer(i8),parameter :: bit_mask1=9223372036854775807
    !    !b'0111111111111111111111111111111111111111111111111111111111111111'
    !    integer(i8),parameter :: bit_mask2=-9223372036854775808   
    !    !b'1000000000000000000000000000000000000000000000000000000000000000'
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
    integer(i8), parameter      :: jmp_c(2) = (/-4707382666127344949_i8,&
        -2852180941702784734_i8/)

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
    self%s(1) = ieor(ieor(rotl(global_time(1), 55), global_time(2)),&
        shiftl(global_time(2), 14))
    self%s(2) = rotl(global_time(2), 36)
  end function next

  !> Helper function for next()
  pure function rotl(x, k) result(res)
    integer(i8), intent(in) :: x
    integer, intent(in)     :: k
    integer(i8)             :: res

    res = ior(shiftl(x, k), shiftr(x, 64 - k))
  end function rotl

end module mod_random
