!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_basis
!! NAME
!! defs_basis
!!
!! FUNCTION
!! This module contains definitions for a number of named constants and
!! physical constants.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2006 ABINIT group (HM, XG,XW)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! Of the named constants,
!! by far the most important are those that define the 'kind' types of
!! virtually all the variables used in a (well-written) FORTRAN 90 code
!! the content of this file is derived from 'Numerical Recipes in Fortran 90'
!! W.H. Press et al., volume 2 of 'Fortran Numerical Recipes', Cambridge
!! University Press, Second Edition (1996), p. 937 and 1361
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module defs_basis

 implicit none

!Keyword 'integer' stands for default integer type
!and may be used whenever integer are presumed to be small

!nb of bytes related to an integer subtype n such as -10^(argument) < n < 10^(argument) (this is standard F90)
 integer, parameter :: i1b=selected_int_kind(2)
 integer, parameter :: i2b=selected_int_kind(4)
 integer, parameter :: i4b=selected_int_kind(9)

!nb of bytes related to default simple-precision real/complex subtypes
!(= 4 for many machine architectures, = 8 for Cray T3E for instance)
!integer, parameter :: sp=kind(1.0)          ! Single precision should not be used
!integer, parameter :: spc=kind((1.0,1.0))

!nb of bytes related to default double-precision real/complex subtypes
!(= 8 for many machine architectures)
 integer, parameter :: dp=kind(1.0d0)
 integer, parameter :: dpc=kind((1.0d0,1.0d0))  ! Complex should not be used presently
                                                ! except for use of libraries

!Example:
! integer, parameter :: urp=selected_real_kind((p=)12,(r=)50)
! real((kind=)urp) :: d
! d=5.04876_urp   ! for a real d with 12 significative digits
! and such as 10^-50 < |d| < 10^50

!To modify sp/spc and / or dp/dpc, insert instructions such as 'dp='
! but do not modify the other declarations in this module

!Default logical type
 integer, parameter :: lgt=kind(.true.)

!The default lengths
 integer, parameter :: fnlen=132    ! maximum length of file name variables
 integer, parameter :: strlen=32000 ! maximum length of input string

!Some constants:
 integer, parameter :: integer_not_used=0
 logical, parameter :: logical_not_used=.true.

!UNIX unit numbers : standard input, standard output, ab_out, and a number
!for temporary access to a file.
 integer, parameter :: std_in=5,ab_in=5  ! generally, the number 5 is directly used
 integer, parameter :: std_out=6         ! generally, the number 6 is directly used
 integer, parameter :: ab_out=7
 integer, parameter :: ab_xml_out = 50  ! this unit is used to print output into an XML file
 integer, parameter :: tmp_unit=9,tmp_unit2=10

!The 3x3 identity matrix
!WARNING : this seem not to work ?!
! integer, dimension(3,3), parameter :: &
!& identity3by3=reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))

!Real constants
 real(dp), parameter :: zero=0._dp
 real(dp), parameter :: one=1._dp
 real(dp), parameter :: two=2._dp
 real(dp), parameter :: three=3._dp
 real(dp), parameter :: four=4._dp
 real(dp), parameter :: five=5._dp
 real(dp), parameter :: six=6._dp
 real(dp), parameter :: seven=7._dp
 real(dp), parameter :: eight=8._dp
 real(dp), parameter :: nine=9._dp
 real(dp), parameter :: ten=10._dp

!Fractionary real constants
 real(dp), parameter :: half=0.50_dp
 real(dp), parameter :: third=one/three
 real(dp), parameter :: quarter=0.25_dp
 real(dp), parameter :: fifth=0.20_dp
 real(dp), parameter :: sixth=one/six
 real(dp), parameter :: seventh=one/seven
 real(dp), parameter :: eighth=0.125_dp
 real(dp), parameter :: ninth=one/nine
 real(dp), parameter :: two_thirds=two*third
 real(dp), parameter :: four_thirds=four*third
 real(dp), parameter :: five_thirds=five*third
 real(dp), parameter :: three_quarters=0.75_dp
 real(dp), parameter :: three_fifth=three/five

!Real constants related to the golden number
 real(dp), parameter :: gold=1.618033988749894848204586834365638117720309179_dp
 real(dp), parameter :: goldenratio=two-gold

!Real constants derived from pi
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 real(dp), parameter :: two_pi=two*pi
 real(dp), parameter :: four_pi=four*pi
 real(dp), parameter :: piinv=one/pi
!The following are not used
!real(dp), parameter :: rad_to_deg=180._dp/pi
!real(dp), parameter :: deg_to_rad=one/rad_to_deg
!real(dp), parameter :: half_pi=pi*half
!real(dp), parameter :: third_pi=pi*third
!real(dp), parameter :: quarter_pi=pi*quarter
!real(dp), parameter :: two_thirds_pi=two_thirds*pi


!Real precision
 real(dp), parameter :: greatest_real = huge(one)
 real(dp), parameter :: smallest_real = -greatest_real
 real(dp), parameter :: tol6= 0.000001_dp
 real(dp), parameter :: tol8= 0.00000001_dp
 real(dp), parameter :: tol10=0.0000000001_dp
 real(dp), parameter :: tol11=0.00000000001_dp
 real(dp), parameter :: tol12=0.000000000001_dp
 real(dp), parameter :: tol14=0.00000000000001_dp

!Real physical constants
!Revised fundamental constants from http://physics.nist.gov/cuu/Constants/index.html
!(from 2002 least squares adjustment)
 real(dp), parameter :: Bohr_Ang=0.5291772108_dp    ! 1 Bohr, in Angstrom
 real(dp), parameter :: Ha_cmm1=219474.6313705_dp  ! 1 Hartree, in cm^-1
 real(dp), parameter :: Ha_eV=27.2113845_dp ! 1 Hartree, in eV
 real(dp), parameter :: Ha_K=Ha_eV*11600.0_dp ! 1Hartree, in Kelvin
 real(dp), parameter :: Ha_THz=6579.683920720_dp ! 1 Hartree, in THz
 real(dp), parameter :: e_Cb=1.60217653d-19 ! minus the electron charge, in Coulomb
 real(dp), parameter :: kb_HaK=8.617343d-5/Ha_eV ! Boltzmann constant in Ha/K
 real(dp), parameter :: amu_emass=1.66053886d-27/9.1093826d-31 ! 1 atomic mass unit, in electronic mass
!This value is 1Ha/bohr^3 in 1d9 J/m^3
!real(dp), parameter :: HaBohr3_GPa=29421.033_dp ! 1 Ha/Bohr^3, in GPa
 real(dp), parameter :: HaBohr3_GPa=Ha_eV/Bohr_Ang**3*e_Cb*1.0d+21 ! 1 Ha/Bohr^3, in GPa
 real(dp), parameter :: Avogadro=6.0221415d23 ! per mole
!This value is 1 Ohm.cm in atomic units
 real(dp), parameter :: Ohmcm=two*pi*Ha_THz*ninth*ten
!real(dp), parameter :: eps0=8.854187817d-12 ! permittivity of free space in F/m
 real(dp), parameter :: eps0=one/(four_pi*0.0000001_dp*299792458.0_dp**2)
 real(dp), parameter :: AmuBohr2_Cm2=e_Cb*1.0d20/(Bohr_Ang*Bohr_Ang)

!Character constants
 character(len=1), parameter :: ch10 = char(10)

!Define fake communicator for sequential abinit
 integer, parameter :: abinit_comm_serial = -12345
 integer, parameter :: abinit_offset = 8
end module defs_basis
!!***
