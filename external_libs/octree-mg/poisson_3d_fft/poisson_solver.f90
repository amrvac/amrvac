
!! Copyright (C) 2002-2007 BigDFT group 
!! This file is distributed under the terms of the
!! GNU General Public License, see ~/COPYING file
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the list of contributors, see ~/AUTHORS 


!!****h* BigDFT/Poisson_Solver
!! NAME
!!    Poisson_Solver
!!
!! FUNCTION
!!    The module of the Poisson Solver.
!!    It must be used in the parent routine. 
!!
!! USAGE
!!    In the main routine in which the Poisson Solver is called
!!    1) The Poisson kernel must be declared as a pointer, then the 
!!       routine createKernel can be called. On exit, the kernel will be allocated and
!!       ready to use. See the documentation of the createKernel routine for more details
!!    2) The correct sizes for allocating the density/potential and the pot_ion arrays
!!       are given from the routine PS_dim4allocation (see routine documentation for details).
!!       Its argument MUST be in agreement with the arguments of the PSolver routine. 
!!       WARNING: No cross-check of the arguments is performed!
!!    3) The PSolver routine can then be called. On exit, the Hartree potential is computed
!!       and summed (following ixc value) to XC and external potential. 
!!       The input array is overwritten. Again, see routine documentation for details.
!!    4) QUICK INSTRUCTION FOR THE IMPATIENT:If you want to use the Poisson Solver in the 
!!       "ordinary" way, for a grid of dimensions nx,ny,nz and grid spacings hx,hy,hz, 
!!       just create the Kernel with
!!           call createKernel(geocode,nx,ny,nz,hx,hy,hz,14,0,1,kernel)
!!       where kernel is a pointer as described above; 
!!       geocode is 'F','S' or 'P' for Free, Surfaces of Periodic BC respectively.
!!       (Beware that for Surfaces BC the isolated direction is y!)
!!       After that you can calculate the potential with
!!           call PSolver(geocode,'G',0,1,nx,ny,nz,0,hx,hy,hz,&
!!                rhopot,kernel,fake_arr,energy,fake_exc,fake_vxc,0.d0,.false.,1)
!!       where:
!!         rhopot    is the density on input, and the electrostatic potential on output
!!                   (an array of dimension(nx,ny,nz))
!!         energy    is the result of 1/2 \int dx rho(x) pot(x)
!!         fake_arr  is an array of dimension(1), untouched
!!         fake_*xc  values of the XC energies, automatically zero in that case
!!
!!       Any other changment of the arguments require reading of the documentation.
!!       See documentations of the Public routines
!!
!!
!! WARNING
!!    This module REQUIRE the module of XC functional from ABINIT, defs_xc, which
!!    require defs_basis and defs_datatypes. 
!!    Such routines are provided inside the abinit directory of this bundle.
!!    They are based on XC functionals present in ABINIT 5.x
!!    If you want to use this Poisson Solver without the XC functionals, you can comment out
!!    the XC part in the PSolver routine
!!    Search for
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    February 2007
!!
!! SOURCE
!!
module Poisson_Solver

  private

  ! Set this to 0 to let processor 0 output some extra information
  integer, public :: iproc_verbose = -1

  !calculate the allocation dimensions
  public :: PS_dim4allocation
  !routine that creates the kernel
  public :: createKernel
  !calculate the poisson solver
  public :: PSolver
  !calculate the allocation dimensions
  public :: P_FFT_dimensions, S_FFT_dimensions, F_FFT_dimensions

contains

  include 'psolver_main.f90'
  include 'build_kernel.f90'
  include 'psolver_base.f90'
  include 'xcenergy.f90'
  include '3dgradient.f90'
  include 'fft3d.f90'
  include 'scaling_function.f90'

end module Poisson_Solver
