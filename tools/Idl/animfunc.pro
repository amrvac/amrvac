;===========================================================================
function animfunc,x,w,func,physics,eqpar,wnames
;
; Written by G. Toth for the Versatile Advection Code
;
; This is the function called by ".r animate" and ".r plotfunc".
;
; "x" array contains the "ndim" components of the coordinates for the grid.
; "w" array contains the "nw" conservative variables for the grid.
; "func" string describes the function to be returned.
; "physics" string defines the meaning of the variables in "w".
; "eqpar" array contains the equation parameters.
; "wnames" string array contains the names for the conservative variables.
;
; The "func" string is interpreted by the following rules:
;
;   If the first character is '-' it is stripped off and the function
;       is multiplied by -1 before returned.
;
;   A single number between 0 and nw-1 returns the variable indexed by
;       that number, e.g. '2' in 3D returns w(*,*,*,2). 
;
;   Variable names in the "wnames" array mean the appropriate variable.
;
;   Function names listed below are calculated and returned.
;
;   Expressions formed from the standard conservative variable names,
;     (i.e. rho m1 m2 m3 e b1 b2 b3), coordinate names (xx, yy, zz, r, phi, z),
;     and equation parameter names (gamma, adiab, eta) and any IDL function
;     and operator are evaluated and returned. 
;
; Examples for valid strings: '3', 'bphi', '-divb', 'r*m1/rho', ...
;
; One can use "animfunc" interactively too, e.g.
; 
; ekin=animfunc(x,w,'(m1^2+m2^2)/rho','hdadiab12')
; pneg=animfunc(x,w,'-pth','mhd22',eqpar)
;
; Note that "eqpar" is needed for the pressure but not for the kinetic energy, 
; while "wnames" is not needed in either case.
;===========================================================================

; In 1D x(n1), in 2D x(n1,n2,2), in 3D x(n1,n2,n3,3)
siz=size(x)
ndim=siz(0)-1
if ndim eq 0 then ndim=1
n1=siz(1)
if ndim gt 1 then n2=siz(2)
if ndim gt 2 then n3=siz(3)
; For 1 variable: w(n1),   w(n1,n2),   w(n1,n2,n3)
; for more      : w(n1,nw),w(n1,n2,nw),w(n1,n2,n3,nw)
siz=size(w)
if siz(0) eq ndim then nw=1 else nw=siz(ndim+1)

if n_elements(wnames) eq 0 then wnames=strarr(nw)

; Check for a negative sign in func
if strmid(func,0,1) eq '-' then begin 
   f=strmid(func,1,strlen(func)-1)
   sign=-1
endif else begin
   f=func
   sign=1
endelse

; Check if f is among the variable names listed in wnames, or if it is a number
for iw=0,nw-1 do $
   if f eq strtrim(string(iw),2) or f eq wnames(iw) then $
   case ndim of
      1:result=w(*,iw)
      2:result=w(*,*,iw)
      3:result=w(*,*,*,iw)
   endcase

if n_elements(result) gt 0 then return,sign*result

; Extract variables from w based on physics
if not keyword_set(physics) then begin
   print,'Error in animfunc: physics=',physics,' is missing'
   retall
endif
physics=strcompress(physics, /remove_all)
case physics of
   'rho11':rho=w
   'rho22':rho=w
   'rho33':rho=w
   'hdadiab11':begin
      rho=w(*,0) & m1=w(*,1)
   end
   'hdadiab12':begin
      rho=w(*,0) & m1=w(*,1) & m2=w(*,2)
   end
   'hdadiab13':begin
      rho=w(*,0) & m1=w(*,1) & m2=w(*,2) & m3=w(*,3)
   end
   'hdadiab22':begin
      rho=w(*,*,0) & m1=w(*,*,1) & m2=w(*,*,2)
   end
   'hdadiab23':begin
      rho=w(*,*,0) & m1=w(*,*,1) & m2=w(*,*,2) & m3=w(*,*,3)
   end
   'hdadiab33':begin
      rho=w(*,*,*,0) & m1=w(*,*,*,1) & m2=w(*,*,*,2) & m3=w(*,*,*,3)
   end
   'hd11':begin
      rho=w(*,0) & m1=w(*,1) & e=w(*,2)
   end
   'hd12':begin
      rho=w(*,0) & m1=w(*,1) & m2=w(*,2) & e=w(*,3)
   end
   'hd13':begin
      rho=w(*,0) & m1=w(*,1) & m2=w(*,2) & m3=w(*,3) & e=w(*,4)
   end
   'hd22':begin
      rho=w(*,*,0) & m1=w(*,*,1) & m2=w(*,*,2) & e=w(*,*,3)
   end
   'hd23':begin
      rho=w(*,*,0) & m1=w(*,*,1) & m2=w(*,*,2) & m3=w(*,*,3) & e=w(*,*,4)
   end
   'hd33':begin
      rho=w(*,*,*,0) & m1=w(*,*,*,1) & m2=w(*,*,*,2) & m3=w(*,*,*,3)
      e=w(*,*,*,4)
   end
   'mhdiso11':begin
      rho=w(*,0) & m1=w(*,1) & b1=w(*,2)
   end
   'mhdiso12':begin
      rho=w(*,0) & m1=w(*,1) & m2=w(*,2) & b1=w(*,3) & b2=w(*,4)
   end
   'mhdiso13':begin
      rho=w(*,0) & m1=w(*,1) & m2=w(*,2) & m3=w(*,3)
      b1=w(*,4) & b2=w(*,5) & b3=w(*,6)
   end
   'mhdiso22':begin
      rho=w(*,*,0) & m1=w(*,*,1) & m2=w(*,*,2) & b1=w(*,*,3) & b2=w(*,*,4)
   end
   'mhdiso23':begin
      rho=w(*,*,0) & m1=w(*,*,1) & m2=w(*,*,2) & m3=w(*,*,3)
      b1=w(*,*,4) & b2=w(*,*,5) & b3=w(*,*,6)
   end
   'mhdiso33':begin
      rho=w(*,*,*,0) & m1=w(*,*,*,1) & m2=w(*,*,*,2) & m3=w(*,*,*,3)
      b1=w(*,*,*,4) & b2=w(*,*,*,5) & b3=w(*,*,*,6)
   end
   'mhd11':begin
      rho=w(*,0) & m1=w(*,1) & e=w(*,2) & b1=w(*,3)
   end
   'mhd12':begin
      rho=w(*,0) & m1=w(*,1) & m2=w(*,2) & e=w(*,3) & b1=w(*,4) & b2=w(*,5)
   end
   'mhd13':begin
      rho=w(*,0) & m1=w(*,1) & m2=w(*,2) & m3=w(*,3) & e=w(*,4)
      b1=w(*,5) & b2=w(*,6) & b3=w(*,7)
   end
   'mhd22':begin
      rho=w(*,*,0) & m1=w(*,*,1) & m2=w(*,*,2) & e=w(*,*,3)
      b1=w(*,*,4) & b2=w(*,*,5)
   end
   'mhd23':begin
      rho=w(*,*,0) & m1=w(*,*,1) & m2=w(*,*,2) & m3=w(*,*,3) & e=w(*,*,4)
      b1=w(*,*,5) & b2=w(*,*,6) & b3=w(*,*,7)
   end
   'mhd33':begin
      rho=w(*,*,*,0) & m1=w(*,*,*,1) & m2=w(*,*,*,2) & m3=w(*,*,*,3)
      e=w(*,*,*,4) & b1=w(*,*,*,5) & b2=w(*,*,*,6) & b3=w(*,*,*,7)
   end
   else: begin
     print,'Error in animfunc: physics=',physics,' is not known'
     retall
   end
endcase

; Extract coordinate variables from x:
; xx,yy,zz for Cartesian and r,phi/z for polar
case ndim of
   1:begin
      xx=x & r=x
   end
   2:begin
      xx=x(*,*,0) & yy=x(*,*,1) & r=xx & z=yy & phi=yy
   end
   3:begin
      xx=x(*,*,*,0) & yy=x(*,*,*,1) & zz=x(*,*,*,2) & r=xx & phi=yy
   end
endcase

; Extract phys (physics without the numbers) and ndir from physics

l=strlen(physics)
phys=strmid(physics,0,l-2)
ndir=long(strmid(physics,l-1,1))
if ndim ne long(strmid(physics,l-2,1)) then begin
   print,'Error in animfunc: physics=',physics,$
         ' is inconsistent with ndim=',ndim
   retall
endif

; Put elements of eqpar into simple variables based on phys
if n_elements(eqpar) ne 0 then case phys of
   'rho': v1=eqpar(0)
   'hdadiab': begin
      gamma=eqpar(0) & adiab=eqpar(1)
   end
   'hd': gamma=eqpar(0)
   'mhdiso': begin
      gamma=eqpar(0) & adiab=eqpar(1) & csound2=adiab & eta=eqpar(2)
   end
   'mhd': begin
      gamma=eqpar(0) & eta=eqpar(1)
   end
endcase

; define momentum and magnetic field squared if appropriate
if n_elements(m1) gt 0 then case ndir of 
    1: mm=m1^2
    2: mm=m1^2+m2^2
    3: mm=m1^2+m2^2+m3^2
endcase
if n_elements(b1) gt 0 then case ndir of 
    1: bb=b1^2
    2: bb=b1^2+b2^2
    3: bb=b1^2+b2^2+b3^2
endcase

;==== Put your function definition(s) below using the variables and eq. params
;     extracted from x, w and eqpar, and select cases by the function name "f"
;     and the "phys, ndim, ndir" variables extracted from "physics".

case 1 of
  ; Velocities
  f eq 'vx' or f eq 'v1' or f eq 'v': result=m1/rho/1.e5
  f eq 'vy' or f eq 'v2'            : result=m2/rho
  f eq 'vz' or f eq 'v3'            : result=m3/rho
  f eq 'calfven1'                   : result=b1/sqrt(rho)
  f eq 'calfven2'                   : result=b2/sqrt(rho)
  f eq 'calfven3'                   : result=b3/sqrt(rho)
  f eq 'calfven'                    : result=sqrt(bb/rho)
  f eq 'maeseta'                    : result=0.1*0.125*xx*sqrt(bb/rho)*exp(-2.0*yy^2/0.125^2/xx^2)
  f eq 'Malfven1'                   : result=m1/b1/sqrt(rho)
  f eq 'Malfven2'                   : result=m2/b2/sqrt(rho)
  f eq 'Malfven3'                   : result=m3/b3/sqrt(rho)
  f eq 'Malfven'                    : result=sqrt(mm/bb/rho)

  ;Thermal pressure, temperature, sound/fast speed or (Alfvenic)Mach number
  f eq 'p' or f eq 'pth' or f eq 'pbeta' or f eq 'T' or f eq 's' or $
  f eq 'logpth' or f eq 'n' or f eq 'logn' or $
  f eq 'logs' or f eq 'logT' or $
  f eq 'eps' or f eq 'csound' or strmid(f,0,4) eq 'mach' or $
  strmid(f,0,5) eq 'cfast' or strmid(f,0,5) eq 'cslow' or $
  strmid(f,0,5) eq 'Mfast' or strmid(f,0,5) eq 'Mslow' $
  :begin
     if phys eq 'rho' then begin
        print,'Error in animfunc: func=',f,$
              ' is not defined for physics=',physics
        retall
     endif

     ; Calculate thermal pressure
     case phys of
        'hdadiab': p=adiab*rho^gamma
        'mhdiso':  p=adiab*rho^gamma
        'hd':      p=(gamma-1)*(e-mm/rho/2)
        'mhd':     p=(gamma-1)*(e-(mm/rho+bb)/2)
     endcase

     ; Calculate gamma*p+bb if that is needed
     if strmid(f,1,4) eq 'fast' or strmid(f,1,4) eq 'slow' $
        or f eq 'machA' then cc=gamma*p+bb

     ;Calculate f from p=p_thermal and aa
     case f of
         'n'     :result=rho/1.67262d-24/1.4d0
         'logn'  :result=alog10(rho/1.67262d-24/1.4d0)
         'pth'   :result=p
         'logpth'   :result=alog10(p)
         'p'     :if phys eq 'mhd' or phys eq 'mhdiso' then result=p+bb/2 $
                  else result=p
         'pbeta' :result=2*p/bb
         'T'     :result=p/rho/1.356e8
         'logT'  :result=alog10(p/rho/1.356e8)
         's'     :result=p/rho^gamma
         'logs'  :result=alog10(p/rho^gamma)
         'eps'   :result=p/(rho*(gamma-1.0))
         'csound':result=sqrt(gamma*p/rho)
         'mach'  :result=sqrt(mm/(gamma*p*rho))
         'mach1' :result=m1/sqrt(gamma*p*rho)
         'mach2' :result=m2/sqrt(gamma*p*rho)
         'mach3' :result=m3/sqrt(gamma*p*rho)
         'machr' :result=(xx*m1+yy*m2)/sqrt((xx^2+yy^2)*gamma*p*rho)
	 'cfast' :result=sqrt(cc/rho)
	 'cfast1':result=sqrt((cc+sqrt(cc^2-4*gamma*p*b1^2))/2/rho)
	 'cfast2':result=sqrt((cc+sqrt(cc^2-4*gamma*p*b2^2))/2/rho)
	 'cfast3':result=sqrt((cc+sqrt(cc^2-4*gamma*p*b3^2))/2/rho)
	 'cslow1':result=sqrt((cc-sqrt(cc^2-4*gamma*p*b1^2))/2/rho)
	 'cslow2':result=sqrt((cc-sqrt(cc^2-4*gamma*p*b2^2))/2/rho)
	 'cslow3':result=sqrt((cc-sqrt(cc^2-4*gamma*p*b3^2))/2/rho)
	 'Mfast' :result=sqrt(mm/cc/rho)
         'machA' :result=sqrt(mm/cc/rho)
	 'Mfast1':result=m1/sqrt((cc+sqrt(cc^2-4*gamma*p*b1^2))/2*rho)
	 'Mfast2':result=m2/sqrt((cc+sqrt(cc^2-4*gamma*p*b2^2))/2*rho)
	 'Mfast3':result=m3/sqrt((cc+sqrt(cc^2-4*gamma*p*b3^2))/2*rho)
	 'Mslow1':result=m1/sqrt((cc-sqrt(cc^2-4*gamma*p*b1^2))/2*rho)
	 'Mslow2':result=m2/sqrt((cc-sqrt(cc^2-4*gamma*p*b2^2))/2*rho)
	 'Mslow3':result=m3/sqrt((cc-sqrt(cc^2-4*gamma*p*b3^2))/2*rho)
      endcase
  end
  f eq 'jz': case ndim of
     1: result=deriv(x,b2)
     2: result=curl(b1,b2,xx,yy)
     3: begin
          print,'Error in animfunc: jz is not implemented for 3D yet'
          retall
        end
  endcase

  ; The functions below now work in 2D only. 
  ; The capitalized versions are valid for Cartesian coordinates only
  f eq 'curlv':  result=curl(m1/rho,m2/rho,xx,yy)
  f eq 'CURLV':  result=diff2(1,m2/rho,xx)-diff2(2,m1/rho,yy)
  f eq 'divv':   result=div(m1/rho,m2/rho,xx,yy)
  f eq 'DIVV':   result=diff2(1,m1/rho,xx)+diff2(2,m2/rho,yy)
  f eq 'j':      result=curl(b1,b2,xx,yy)
  f eq 'J':      result=diff2(1,b2,xx)-diff2(2,b1,yy)
  f eq 'j_rz':   result=curl_rz(b1,b2,r,z)
  f eq 'j_RZ':   result=-curl(b1,b2,r,z)
  f eq 'J_rz':   result=diff2(2,b1,z)-diff2(1,b2,r)
  f eq 'j_rp':   result=(diff2(1,r*b2,r)-diff2(2,b1,phi))/r
  f eq 'divb':   result=div(b1,b2,xx,yy)
  f eq 'DIVB':   result=diff2(1,b1,xx)   +diff2(2,b2,yy)
  f eq 'divb4':  result=diff4(1,b1,xx)   +diff4(2,b2,yy)
  f eq 'divb_rz':result=div_rz(b1,b2,r,z)
  f eq 'DIVB_rz':result=diff2(1,b1*r,r)/r+diff2(2,b2,z)
  f eq 'divb_rp':result=diff2(1,b1*r,r)/r+diff2(2,b2,phi)/r
  f eq 'divb_CT':begin
      ; Cell corner centered div B definition, 0 for Constrained Transport
      result=dblarr(n1,n2)
      result(0:n1-2,0:n2-2)=(b1(1:n1-1,0:n2-2)+b1(1:n1-1,1:n2-1)  $
                            -b1(0:n1-2,0:n2-2)-b1(0:n1-2,1:n2-1)) $
                           /(xx(1:n1-1,0:n2-2)-xx(0:n1-2,0:n2-2)) $
                           +(b2(0:n1-2,1:n2-1)+b2(1:n1-1,1:n2-1)  $
                            -b2(0:n1-2,0:n2-2)-b2(1:n1-1,0:n2-2)) $
                           /(yy(0:n1-2,1:n2-1)-yy(0:n1-2,0:n2-2))
      result(n1-1,*)=result(n1-2,*)
      result(*,n2-1)=result(*,n2-2)
      result=result/2
  end
  f eq 'divb_CD':begin
      ;Central diff. div B definition on generalized grid, zero for CD schemes
      result=dblarr(n1,n2)
      ;b_xi = dy/deta*B_x-dx/deta*B_y
      b_xi = (yy(*,2:n2-1)-yy(*,0:n2-3))*b1(*,1:n2-2) - $
             (xx(*,2:n2-1)-xx(*,0:n2-3))*b2(*,1:n2-2)
      ;b_eta=-dy/dxi*B_x+dx/dxi*B_y
      b_eta=-(yy(2:n1-1,*)-yy(0:n1-3,*))*b1(1:n1-2,*) + $
             (xx(2:n1-1,*)-xx(0:n1-3,*))*b2(1:n1-2,*)

      ; divb_CD=(db_xi/dxi+db_eta/deta)/(dx/dxi*dy/deta-dy/dxi*dx/deta)
      result(1:n1-2,1:n2-2)=(b_xi(2:n1-1,*) -b_xi(0:n1-3,*)   $
                            +b_eta(*,2:n2-1)-b_eta(*,0:n2-3)) $
         /((xx(2:n1-1,1:n2-2)-xx(0:n1-3,1:n2-2)) $
          *(yy(1:n1-2,2:n2-1)-yy(1:n1-2,0:n2-3)) $
          -(yy(2:n1-1,1:n2-2)-yy(0:n1-3,1:n2-2)) $
          *(xx(1:n1-2,2:n2-1)-xx(1:n1-2,0:n2-3)))
 
      ; boundaries
      result(0,*)=result(1,*)
      result(*,0)=result(*,1)
      result(n1-1,*)=result(n1-2,*)
      result(*,n2-1)=result(*,n2-2)
  end
  ; Magnetic vector potential A in slab or A*r in cylindrical symmetry,
  ; integrating Bx,By or r*Br,r*Bz respectively. The contour lines of
  ; A (or A*r in cylindrical symmetry) are parallel to the magntic field.
  ; The density of contourlines is proportional to B (or B*r).
  f eq 'A' or f eq 'A_r' or f eq 'AA' or f eq 'AA_r' $
  or f eq 'B' or f eq 'B_r':begin
     bx=b1
     by=b2
     if f eq 'A_r' or f eq 'AA_r' or f eq 'B_r' then begin
        ; For axial symmetry Cartesian_curl(r*A)=(B_r*r,B_z*r)
        bx=r*bx
        by=r*by
     endif
     result=dblarr(n1,n2)
     if f eq 'A' or f eq 'A_r' or f eq 'B_r' then begin
             ; Integrate along the first row
             for i1=1,n1-1 do result(i1,0)=result(i1-1,0) $
               +(by(i1,0)+by(i1-1,0))*(xx(i1,0)-xx(i1-1,0))*0.5 $
               -(bx(i1,0)+bx(i1-1,0))*(yy(i1,0)-yy(i1-1,0))*0.5
             ; Integrate all columns vertically
             for i2=1,n2-1 do result(*,i2)=result(*,i2-1) $
               +(by(*,i2)+by(*,i2-1))*(xx(*,i2)-xx(*,i2-1))*0.5 $
               -(bx(*,i2)+bx(*,i2-1))*(yy(*,i2)-yy(*,i2-1))*0.5
     endif else begin
             ; Integrate first column vertically
             for i2=1,n2-1 do result(0,i2)=result(0,i2-1) $
                -(bx(0,i2)+bx(0,i2-1))*(yy(0,i2)-yy(0,i2-1))*0.5 $
                +(by(0,i2)+by(0,i2-1))*(xx(0,i2)-xx(0,i2-1))*0.5
             ; Integrate all rows
             for i1=1,n1-1 do result(i1,*)=result(i1-1,*)$
                -(bx(i1,*)+bx(i1-1,*))*(yy(i1,*)-yy(i1-1,*))*0.5 $
                +(by(i1,*)+by(i1-1,*))*(xx(i1,*)-xx(i1-1,*))*0.5
     endelse
     ; Smooth results for B
     if f eq 'B' or f eq 'B_r' then result=smooth(result,4)
  end

  ; Some functions used in specific problems
  f eq 'schlier':  begin
       result=sqrt((diff2(1,rho,xx))^2+(diff2(2,rho,yy))^2)
  end
  ; Some functions used in specific problems
  f eq 'schlier2':  begin
       gradrho=sqrt((diff2(1,rho,xx))^2+(diff2(2,rho,yy))^2)
       kk=5.0
       kk0=0.01
       kk1=1.0
       grhomax=300.0
       kk=5.0
       kk0=0.01
       kk1=1.0
       grhomax=900.0
       result=exp(-kk*(gradrho-kk0*grhomax)/(kk1*grhomax-kk0*grhomax))
  end
  ; Some functions used in specific problems
  f eq 'schlier2p':  begin
       p=(gamma-1)*(e-(mm/rho+bb)/2)
       gradrho=sqrt((diff2(1,p,xx))^2+(diff2(2,p,yy))^2)
       kk=5.0
       kk0=0.01
       kk1=1.0
       grhomax=900.0
       result=exp(-kk*(gradrho-kk0*grhomax)/(kk1*grhomax-kk0*grhomax))
  end

  ; Powell's source terms in 2 and 2.5D MHD

  f eq 'v1divb': result=div(b1,b2,xx,yy)*m1/rho
  f eq 'v2divb': result=div(b1,b2,xx,yy)*m2/rho
  f eq 'v3divb': result=div(b1,b2,xx,yy)*m3/rho
  f eq 'b1divb': result=div(b1,b2,xx,yy)*b1
  f eq 'b2divb': result=div(b1,b2,xx,yy)*b2
  f eq 'b3divb': result=div(b1,b2,xx,yy)*b3
  f eq 'bv2divb': result=div(b1,b2,xx,yy)*(m1*b1+m2*b2)/rho
  f eq 'bv3divb': result=div(b1,b2,xx,yy)*(m1*b1+m2*b2+m3*b3)/rho

  ; These can be used for checking the velocity defined Courant condition
  f eq 'vx/dx':  begin
     dx=fltarr(n1,n2)
     dx(1:n1-2,*)=(xx(2:n1-1,*)-xx(0:n1-3,*))/2
     dx(0,*)=dx(1,*) & dx(n1-1,*)=dx(n1-2,*)
     result=m1/rho/dx
  end
  f eq 'vy/dy':  begin
     dy=fltarr(n1,n2)
     dy(*,1:n2-2)=(yy(*,2:n2-1)-yy(*,0:n2-3))/2
     dy(*,0)=dy(*,1) & dy(*,n2-1)=dy(*,n2-2)
     result=m2/rho/dy
  end

  f eq 'logrho':result=alog10(rho)
  f eq 'bern2': result=gamma*e/rho+(1-gamma)/2*(m1^2+m2^2)/rho^2-0.5/(r-1)
  f eq 'bern3': result=gamma*e/rho+(1-gamma)/2*(m1^2+m2^2+m3^2)/rho^2-0.5/(r-1)

  f eq 'Ez':begin
     result=(m1*b2-m2*b1)/rho-eqpar(1)*(diff2(1,b2,xx)-diff2(2,b1,yy))
  end
  f eq 'Erat':begin
     result=(((m1*b2-m2*b1)/rho)/(-eqpar(1)*(diff2(1,b2,xx)-diff2(2,b1,yy))))
  end
  f eq 'EvB':begin
     result=(m1*b2-m2*b1)/rho
  end
  f eq 'EetaJ':begin
     result=-eqpar(1)*(diff2(1,b2,xx)-diff2(2,b1,yy))
  end

  f eq 'dpdz':begin
     p=(gamma-1)*(e-0.5*(m1^2+m2^2+m3^2)/rho)
     result= diff2(2,p,z)
  end
  f eq 'gz':begin
     Phi=0.5/(sqrt(r^2+z^2)-1)
     result= diff2(2,Phi,z)
  end
  f eq 'dmzdt':begin
     p=(gamma-1)*(e-0.5*(m1^2+m2^2+m3^2)/rho)
     Phi=0.5/(sqrt(r^2+z^2)-1)
     result= -diff2(2,p,z)+rho*diff2(2,Phi,z)
  end
  f eq 'DMZDT':begin
     br=b1
     bz=b2
     pth=(gamma-1)*(e-(br^2+bz^2)/2)
     j=diff2(1,bz,r)-diff2(2,br,z)
     result= -diff2(2,pth,z)+eqpar(3)*w(*,*,0)-j*br
  end
  f eq 'dmrdt':begin
     br=b1
     bz=b2
     pth=(gamma-1)*(e-(br^2+bz^2)/2)
     j=diff2(1,bz,r)-diff2(2,br,z)
     result= -diff2(1,pth,r)+j*bz
  end
  f eq 'TT':begin
     rho=w(*,*,0)
     pr=w(*,*,3)
     result= pr/rho
  end

  ;If "f" has not matched any function yet, try evaluating it as an expression
  else:begin
     if not execute('result='+f) then begin
        print,'Error in animfunc: cannot evaluate function=',func
        retall
     endif
  end
endcase

if n_elements(result) gt 0 then begin
   return,sign*result
endif else begin
   print,'Error in animfunc: function=',func,' was not calculated ?!'
   retall
endelse

end
