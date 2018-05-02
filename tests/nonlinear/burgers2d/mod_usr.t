module mod_usr
  use mod_nonlinear

  implicit none

contains

  subroutine usr_init()
    use mod_variables

    usr_init_one_grid => initonegrid_usr

    call nonlinear_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    logical::          mask(ix^S)
    double precision:: rr2(ix^S)
    double precision:: xc^D,rad,nsig,rhosquare


   nsig=one
   rhosquare=1.0d0/exp(-nsig**2)
   w(ix^S,rho_) =rhosquare*exp(-nsig**2)
   rad=0.2d0
   xc1=2.5d0
   {^NOONED
   xc2=2.5d0
   }
   {^IFTHREED
   xc3=2.5d0
   }
   rr2(ix^S)={^D&(x(ix^S,^D)-xc^D)**2|+}
   mask(ix^S)=(rr2(ix^S)<=(nsig*rad)**2)
   where(mask(ix^S))
       w(ix^S,rho_)     = rhosquare*exp(-rr2(ix^S)/rad**2)
   endwhere
   xc1=1.5d0
   {^NOONED
   xc2=3.5d0
   }
   {^IFTHREED
   xc3=3.5d0
   }
   rr2(ix^S)={^D&(x(ix^S,^D)-xc^D)**2|+}
   mask(ix^S)=(rr2(ix^S)<=(nsig*rad)**2)
   where(mask(ix^S))
       w(ix^S,rho_)     = rhosquare*exp(-rr2(ix^S)/rad**2)
   endwhere
   xc1=3.5d0
   {^NOONED
   xc2=1.5d0
   }
   {^IFTHREED
   xc3=1.5d0
   }
   rr2(ix^S)={^D&(x(ix^S,^D)-xc^D)**2|+}
   mask(ix^S)=(rr2(ix^S)<=(nsig*rad)**2)
   where(mask(ix^S))
       w(ix^S,rho_)     = rhosquare*exp(-rr2(ix^S)/rad**2)
   endwhere

  end subroutine initonegrid_usr

end module mod_usr

