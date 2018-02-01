module krome_user_commons
  implicit none

   ! *************************************************************
   !  This file has been generated with:
   !  KROME 14.08.dev on 2018-01-25 18:47:28
   !  Changeset 9a6e335
   !  see http://kromepackage.org
   ! 
   !  Written and developed by Tommaso Grassi and Stefano Bovino
   ! 
   !  Contributors:
   !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle, 
   !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini, 
   !  E.Tognelli
   !  KROME is provided "as it is", without any warranty.
   ! *************************************************************


  !user can add here the common variables needed
  !for rate calculation (e.g. optical depth, CR rate,
  !pressure, density, ...)
  real*8::tau,zrate,pah_size,gas_dust_ratio,krome_redshift
  real*8::krome_J21,krome_J21xr,user_tff,user_omega,user_nu
  real*8::krome_invdvdz !inverse of |velocity gradient| 1/abs(dv/dz)

  !$omp threadprivate(user_tff,krome_invdvdz)

contains

  !**********************
  !user can add here the functions he/she needs for
  !rate calculations (Kooij funcion provided as example)
  function kooij(kalpha,kbeta,kgamma,Tgas)
    real*8::kooij
    real*8,intent(in)::kalpha,kbeta,kgamma,Tgas
    kooij = kalpha*(Tgas/3d2)**kbeta*exp(-kgamma/Tgas)
  end function kooij


  !**********************
  !prototype for photo xsec kernel function.
  !this is called by xsec interpolator when resampling
  ! xsec produced by python into photobins at runtime
  function fxinterp(energy)
    implicit none
    real*8,intent(in)::energy
    real*8::fxinterp

    fxinterp = energy**2

  end function fxinterp

end module krome_user_commons


