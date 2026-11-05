#:mute
#:include "../mod_gpu_directives.fpp"
#:endmute

module srmhd_con2prim
use mod_physics_vars

implicit none
real(dp), parameter, public :: h0=1.0_dp


!this module implements a conserved to primitive variable
!scheme from Kastaun et al. https://arxiv.org/pdf/2005.01821
!the code will be labled with equation and section numbers see the paper for reference

contains 

!=============================================================================
!
! This section implements brents root solver for use in the con2prim algorithm
!
!============================Brent Root Solver=================================
    pure function quad_inter(x1, f1, f2, f3) result(quad)
        ${GPU_ROUTINE_SEQ()}$
        real(dp), intent(in) :: x1, f1, f2, f3
        real(dp) :: quad

        quad=x1*f2*f3/(f1-f2)/(f1-f3)
    end function quad_inter

    pure subroutine swap(a,b)
        ${GPU_ROUTINE_SEQ()}$
        real(dp), intent(inout) :: a, b
        real(dp) :: temp
        temp = a
        a = b
        b = temp
    end subroutine swap 

    #: def brent_template(name, func, npar)
    pure function ${name}$(mu_plus, params) result(b)
        implicit none
        ${GPU_ROUTINE_SEQ()}$
        real(dp)  :: s
        real(dp), intent(in) :: params(${npar}$)
        real(dp), intent(in) :: mu_plus
        real(dp) :: a, b, c, d
        real(dp) :: fa, fb, fc, fs
        real(dp) :: tol
        logical :: bound_error
        logical :: mflag
        integer :: i, maxiter
        real(dp):: delta




        bound_error=.false.
        maxiter=100

        a=0
        b=mu_plus

        tol=1e-15_dp

        delta = 2.0_dp * epsilon(1.0_dp) * abs(b) + 0.5_dp * tol

        call ${func}$(a, params, fa)


        call ${func}$(b, params, fb)


        ! if (fa*fb>0) then
        !     b=b+1e-10_dp
        !     fb=${func}$(b, params)
        ! end if

        if (fa*fb>0) then
            bound_error=.true.
            if (abs(fb)>abs(fa)) then
                b=a
            end if
        end if

        if (abs(fa) < abs(fb)) then
            call swap(a, b)
            call swap(fa, fb)
        end if
        
        c=a
        fc=fa
        mflag=.true.

        do i=1, maxiter
            if (bound_error) then
                exit
            end if

            if (fa/=fc .and. fb/=fc) then
                s= quad_inter(a, fa, fb, fc) + quad_inter(b, fb, fa, fc) + quad_inter(c, fc, fb, fa)
            else
                s = b-fb*(b-a)/(fb-fa)
            end if

            if (.not. ( min((3.0_dp*a+b)/4.0_dp, b) < s .and. s < max((3.0_dp*a+b)/4.0_dp, b) ) ) then
                s=(a+b)/2
                mflag=.true.
            elseif (mflag) then

                if (abs(s-b)>=abs(c-b)/2) then
                    s=(b+a)/2
                    mflag=.true.
                else
                    mflag=.false.
                end if
            else

                if (abs(s-b) >= abs(c-d)/2 ) then
                    s=(b+a)/2
                    mflag=.true.
                else
                    mflag=.false.
                end if
            end if

            call ${func}$(s, params, fs)


            d=c
            c=b
            fc=fb

            


            if (fa*fs < 0) then
                b=s
                fb=fs
            else
                a=s
                fa=fs
            end if


            if (abs(fa) < abs(fb)) then
            call swap(a, b)
            call swap(fa, fb)
            end if


            delta = 2.0_dp * epsilon(1.0_dp) * abs(b) + 0.5_dp * tol

            if (fb==0 .or. abs(a-b) < delta*2) then
                exit
            end if

        end do


    end function ${name}$

    #:enddef

    $:brent_template("brent","master_function", 5)
    $:brent_template("aux_brent", "aux_f", 4)

!===========================================================================================
!
!           Implements Root solver TOMS748 https://doi.org/10.1145/210089.210111
!
!==================================TOMS748==================================================



    pure subroutine bracket(a, b, c, fa, fb, fc, d, fd)
        implicit none
        ${GPU_ROUTINE_SEQ()}$
        real(dp), intent(inout) :: a, b, fa, fb 
        real(dp), intent(in) :: c, fc
        real(dp), intent(out) :: d, fd

        if (sign(1.0_dp, fa) /= sign(1.0_dp, fc)) then
            d=b
            fd=fb
 
            b=c
            fb=fc      
        else
            d=a
            fd=fa

            a = c
            fa = fc
        end if
    end subroutine bracket

    pure function cubic_int(a, b, f1, f2, f3) result(q)
        implicit none
        ${GPU_ROUTINE_SEQ()}$
        real(dp), intent(in) :: a, b, f1, f2, f3
        real(dp) :: q
        q = (a-b)*f1/(f2-f3)
    end function cubic_int

    pure function secant(a, b, fa, fb) result(f)
        implicit none
        ${GPU_ROUTINE_SEQ()}$
        real(dp), intent(in) :: a, b, fa, fb
        real(dp) :: f
        f= (fb - fa)/(b - a)
    end function secant


    pure function ipzero(a, b, c, d, fa, fb, fc, fd) result(x_bar)
        implicit none
        ${GPU_ROUTINE_SEQ()}$
        real(dp), intent(in) :: a, b, c, d, fa, fb, fc, fd
        real(dp) :: q11, q21, q31, q22, q32, d21, d31, d32, q33
        real(dp) :: x_bar

        q11 = cubic_int(c, d, fc, fd, fc)

        q21 = cubic_int(b, c, fb, fc, fb)

        q31 = cubic_int(a, b, fa, fb, fa)

        d21 = cubic_int(b, c, fc, fc, fb)

        d31 = cubic_int(a, b, fb, fb, fa)

        q22 = cubic_int(d21, q11, fb, fd, fb)

        q32 = cubic_int(d31, q21, fa, fc, fa)

        d32 = cubic_int(d31, q21, fc, fc, fa)

        q33 = cubic_int(d32, q22, fa, fd, fa)

        x_bar = a + q31 + q32 + q33

    end function ipzero


    pure function newt_quad(a, b, d, fa, fb, fd, k) result(r)
        implicit none
        ${GPU_ROUTINE_SEQ()}$
        real(dp), intent(in) :: a, b, d, fa, fb, fd
        integer, intent(in)  :: k
        real(dp) :: r
        real(dp) :: big_a, big_b
        integer :: i
        real(dp) :: den

        big_b = secant(a, b, fa, fb)

        big_a = secant(a, d, big_b, secant(b, d, fb, fd))

        if (big_a == 0) then

            r = a - fa/big_b

        else
            if (big_a*fa>0) then
                r = a
            else
                r = b
            end if 

            do i=1, k 

                den=big_b + big_a*(2*r - a - b)
                if (den==0) then
                    r=(a+b)/2
                    exit
                end if

                r = r - (fa+ big_b*(r-a)+ big_a*(r-a)*(r-b))/(den)


            end do

        end if

    end function newt_quad

    pure function check_eq(fa, fb, fd, fe) result(equal)
        implicit none
        ${GPU_ROUTINE_SEQ()}$
        real(dp), intent(in) :: fa, fb, fd, fe
        logical :: equal

        equal=.false.

        if (fa==fb .or. fa==fd .or. fa==fe .or. fb==fd .or. fb==fe  .or. fd==fe) then
            equal=.true.
        end if    
    end function check_eq


    #: def TOMS748_template(name, func, npar)
    pure function ${name}$(mu_plus, params) result(out)
        implicit none
        ${GPU_ROUTINE_SEQ()}$
        real(dp), intent(in) :: params(${npar}$)
        real(dp), intent(in) :: mu_plus
        real(dp) :: a, b, c, d, e, u
        real(dp) :: fa, fb, fc, fd, fe, fu
        real(dp) :: tol
        logical :: bound_error
        logical :: success
        integer :: i, maxiter
        real(dp):: delta
        real(dp) :: old_a, old_b
        real(dp) :: eps
        real(dp) :: out

        eps=epsilon(1.0_dp)

        success=.false.
        bound_error=.false.
        maxiter=20

        a=0
        b=mu_plus

        e=0
        fe=0

        tol=1e-18_dp
        ! 4.1.1
        call ${func}$(a, params, fa)



        call ${func}$(b, params, fb)


        if (fa*fb>0) then
            bound_error=.true.
            if (abs(fb)>abs(fa)) then
                out=a
            else
                out=b
            end if
        
        elseif (fa==0) then
            out=a
            success=.true.

        elseif(fb==0) then
            out=b
            success=.true.
            
        else

            c = a - fa/secant(a, b, fa, fb)

            if (abs(c-a)<eps*32) then
                c = c + 100*eps
            elseif (abs(c-b)<eps*32) then
                c = c - 100*eps
            end if

            if ((c-a)*(c-b)>=0) then
                c = (a+b)/2
            end if

            call ${func}$(c, params, fc)


            ! 4.1.2
            if (fc == 0) then
                success = .true.
                out = c
            end if
            call bracket(a, b, c, fa, fb, fc, d, fd)

        end if



        do i=1, maxiter
            if (bound_error .or. success) then
                exit
            end if
            !4.2.3
            if (i==1 .or. check_eq(fa, fb, fd, fe)) then
                c = newt_quad(a, b, d, fa, fb, fd, 2)
            else
                c = ipzero(a, b, d, e, fa, fb, fd, fe)
                ! if step leaves bracket use newton step instead
                if ((c-a)*(c-b)>=0) then
                    c = newt_quad(a, b, d, fa, fb, fd, 2)
                end if
            end if 


            if ((c-a)*(c-b)>=0) then
                c = (a+b)/2
            end if

            call ${func}$(c, params, fc)


            old_a = a
            old_b = b

            !4.2.4
            e = d
            fe = fd

            if (fc == 0) then
                out=c 
                exit
            end if
            call bracket(a, b, c, fa, fb, fc, d, fd)

            if (abs(fa)<abs(fb)) then
                delta = 2*abs(a)*eps+tol
                out=a
            else
                delta = 2*abs(b)*eps+tol
                out=b
            end if
            if (abs(b-a)<=2*delta) exit

            !4.2.5
            if (check_eq(fa, fb, fd, fe))  then
                c = newt_quad(a, b, d, fa, fb, fd, 3)
            else
                c = ipzero(a, b, d, e, fa, fb, fd, fe)
                ! if step leaves bracket use newton step instead
                if ((c-a)*(c-b)>=0) then
                    c = newt_quad(a, b, d, fa, fb, fd, 3)
                end if
                
            end if 


            if ((c-a)*(c-b)>=0) then
                c = (a+b)/2
            end if

            call ${func}$(c, params, fc)


            !4.2.6

            if (fc == 0) then
                out=c 
                exit
            end if

            e = d
            fe = fd            
            call bracket(a, b, c, fa, fb, fc, d, fd)

            if (abs(fa)<abs(fb)) then
                delta = 2*abs(a)*eps+tol
                out=a
            else
                delta = 2*abs(b)*eps+tol
                out=b
            end if
            if (abs(b-a)<=2*delta) exit

            !4.1.5
            if (abs(fa)<abs(fb)) then
                u=b
                fu=fb

                if (abs(c-u)<eps*32) then
                    u = u - abs(u)*tol
                end if
            else 
                u=a
                fu=fa
                if (abs(c-u)<eps*32) then
                    u = u + abs(u)*tol
                end if
            end if



            !4.1.6
            c = u - 2 /secant(a, b, fa, fb)*fu

            !4.1.7
            if (abs(c-u)>0.5*(b-a)) then
                c = (a+b)/2
            end if


            if ((c-a)*(c-b)>=0) then
                c = (a+b)/2
            end if

            call ${func}$(c, params, fc)




            e = d
            fe = fd
            !4.1.8
            if (fc == 0) then
                out=c 
                exit
            end if
            call bracket(a, b, c, fa, fb, fc, d, fd)

            if (abs(fa)<abs(fb)) then
                delta = 2*abs(a)*eps+tol
                out=a
            else
                delta = 2*abs(b)*eps+tol
                out=b
            end if
            if (abs(b-a)<=2*delta) exit

            !4.1.9
            !if bracket hasnt closed enough default to bisection
            if (abs(b-a)>0.5*abs(old_b - old_a)) then

                e = d 
                fe = fd
                c=(a+b)/2
                call ${func}$(c, params, fc)


                if (fc == 0) then
                    out=c 
                    exit
                end if
                call bracket(a, b, c, fa, fb, fc, d, fd)

                if (abs(fa)<abs(fb)) then
                    delta = 2*abs(a)*eps+tol
                    out = a
                else
                    delta = 2*abs(b)*eps+tol
                    out = b
                end if
                if (abs(b-a)<=2*delta) exit
            end if
        end do
        
    end function ${name}$

    #:enddef

    
    $:TOMS748_template("toms", "master_function", 5)
    $:TOMS748_template("aux_toms", "aux_f", 4)
            



!=============================================================================
!
! This is the con2prim procedure
!
!==========================Kastaun con2prim===================================



    pure function ideal_eos(rho,eps) result(p)
        ${GPU_ROUTINE_SEQ()}$
        real(dp), intent(in) :: rho, eps
        real(dp) ::p

        p=(srmhd_gamma-1)*rho*eps
    end function ideal_eos

    ! see http://dx.doi.org/10.1016/j.jcp.2011.01.020 equation 5
    pure function synge_eos(rho,eps) result(p)
        ${GPU_ROUTINE_SEQ()}$
        real(dp), intent(in) :: rho, eps
        real(dp) ::p
        
        p=(srmhd_gamma-1)*rho/2*eps*(2 + eps)/(1 + eps)
    end function synge_eos



    pure subroutine master_function(mu, params, f) 
        implicit none
        ${GPU_ROUTINE_SEQ()}$
        real(dp), intent(in) :: mu
        real(dp) :: d, tau, s_sqr, b_sqr, s_dot_b
        real(dp), intent(in) :: params(5)
        real(dp)  ::  q, x
        real(dp) :: r_sqr, myb_sqr, r_dot_b
        real(dp) :: r_bar_sqr, q_bar
        real(dp) :: v0_sqr
        real(dp) :: v_hat_sqr, w_hat
        real(dp) :: rho_hat, eps_hat
        real(dp) :: p_hat
        real(dp) :: a_hat
        real(dp) :: nu_a, nu_b
        real(dp) :: nu_hat, mu_hat
        real(dp), intent(out) :: f

        
        d = params(1)
        tau = params(2)
        s_sqr = params(3)
        b_sqr = params(4) 
        s_dot_b = params(5)

        r_sqr=s_sqr/d**2 !(23)

        q=tau/d !(23)

        myb_sqr=b_sqr/d !(24)

        r_dot_b=s_dot_b !necessary for (25)

        v0_sqr=r_sqr/(h0**2 + r_sqr) ! (32)(33)

        x=1.0_dp/(1.0_dp + mu*myb_sqr) !(26)

        r_bar_sqr= r_sqr*x**2 + mu*x*(1+x)*r_dot_b**2/d**3 !(38)

        q_bar= q - 0.5_dp*myb_sqr - 0.5_dp*mu**2*x**2*(r_sqr*myb_sqr-r_dot_b**2/d**3) !(39)

        v_hat_sqr=min(v0_sqr, mu**2*r_bar_sqr) !(40)


        w_hat=1/sqrt(1 - v_hat_sqr) !(40)

        rho_hat=d/w_hat !(41)


        eps_hat=w_hat*(q_bar - mu*r_bar_sqr) + v_hat_sqr*w_hat**2/(1 + w_hat) !(42)

        !add a limiter for valid eps and rho
        if (eps_hat<0) then
            eps_hat=0
        end if


        if (.not. srmhd_eos) then 
            p_hat = ideal_eos(rho_hat, eps_hat) !(43)
        else
            p_hat = synge_eos(rho_hat, eps_hat) !(43)
        end if


        a_hat=p_hat/rho_hat/(1 + eps_hat) !(43)

        nu_a=(1 + a_hat)*(1 + eps_hat)/w_hat !(46)

        nu_b=(1 + a_hat)*(1 + q_bar - mu*r_bar_sqr) !(47)
         
        nu_hat=max(nu_a, nu_b) !(48)

        mu_hat = 1.0_dp/(nu_hat + mu*r_bar_sqr) !(45)

        f=mu - mu_hat!(44)

    end subroutine master_function



        


    pure subroutine aux_f(mu, params, f) 
        implicit none
        ${GPU_ROUTINE_SEQ()}$
        real(dp), intent(in) :: mu
        real(dp), intent(in) :: params(4)
        real(dp)             :: d, s_sqr, s_dot_b, b_sqr
        real(dp)             :: mymu
        real(dp)             :: r_dot_b, myb_sqr, myd
        real(dp)             :: r_bar_sqr, r_sqr, x
        real(dp), intent(out):: f
        real(dp) :: root

        d = params(1)
        s_sqr = params(2)
        b_sqr = params(3)
        s_dot_b = params(4)

        mymu=mu
        r_dot_b=s_dot_b
        myd=d
        myb_sqr=b_sqr/myd
        r_sqr=s_sqr/myd**2

        x=1/(1+mu*myb_sqr)

        r_bar_sqr= r_sqr*x**2 + mu*x* (1+x) *r_dot_b**2/myd**3 !(38)

        root = sqrt(h0**2+r_bar_sqr)
        
        f=mu*root-1

    end subroutine aux_f


    #:def con2prim_template(name, solver, aux_solver)
    pure subroutine ${name}$(mu, d, tau, s_sqr, b_sqr, s_dot_b)
        implicit none
        ${GPU_ROUTINE_SEQ()}$
        !real(dp), intent(out) :: lfac
        !real(dp), intent(out), optional:: rho, eps, p, h
        !logical, intent(out), optional :: success
        !integer, intent(out), optional :: iter, aux_iter
        !real(dp), intent(out) :: vi(3)
        real(dp), intent(in)  :: d, tau
        real(dp), intent(in)  :: s_sqr, b_sqr, s_dot_b
        real(dp)              :: myd, mytau
        !real(dp)              :: mysi(3), mybi(3)
        !real(dp)              :: s_sqr, b_sqr, s_dot_b
        real(dp), intent(inout) :: mu
        real(dp)              :: mu_plus
        real(dp)              :: q_bar, r_bar_sqr
        real(dp)              :: x
        real(dp)              :: params(5)
        real(dp)              :: aux_params(4)
        

        myd=d
        mytau=tau


        aux_params=[myd, s_sqr, b_sqr, s_dot_b]

        if (s_sqr/myd**2 < h0**2) then
            mu_plus=1/h0
        else
            mu_plus=${aux_solver}$(1/h0, aux_params)
        end if

        
        params=[myd, mytau, s_sqr, b_sqr, s_dot_b]

        mu = ${solver}$(mu_plus, params)


    end subroutine ${name}$
    #:enddef 



    $:con2prim_template("con2prim", "brent", "aux_brent")

    !$:con2prim_template("con2prim", "toms", "aux_toms")


   

    

end module srmhd_con2prim



