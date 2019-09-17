;=============================================================================
PRO VEL,U,V,XX,YY,NVECS=nvecs,MAXVAL=maxval,LENGTH=length,HEAD=head,$
		  NSTEPS=nsteps,X0=x0,DYNAMIC=dynamic,NOERASE=noerase,$
		  XXOLD=xxold,YYOLD=yyold,TRIANGLES=triangles,SEED=seed
;+
; NAME:
;	VEL
;
; PURPOSE:
;	Draw a velocity (flow) field with arrows following the field 
;	proportional in length to the field strength.  Arrows are composed 
;	of a number of small segments that follow the streamlines.
;	When called repeatedly the arrows can move along the streamlines.
;	Also works with irregular grids using triangulation.
;
; CATEGORY:
;	Graphics, two-dimensional.
;
; CALLING SEQUENCE:
;	VEL, U, V [, XX, YY]
;
; INPUTS:
;	U:	The X component at each point of the vector field.  This 
;		parameter must be a 2D array.
;
;	V:	The Y component at each point of the vector field.  This 
;		parameter must have the same dimensions as U.
;
; OPTIONAL INPUT PARAMETERS:
;       XX:     X-coordinate array (2D), default is uniform grid with dx=1.
;
;       YY:     Y-coordinate array (2D), default is uniform grid with dy=1.
;
; KEYWORD PARAMETERS:
;	NVECS:	The number of vectors (arrows) to draw.  Default is 200.
;
;	MAXVAL:	The value of the highest velocity that can possibly occur.
;		Default is calculated as max(sqrt(U^2+V^2)).
;
;	LENGTH:	The length of the arrow corresponding to MAXVAL as a fraction 
;		of the diagonal length of the X-Y box. The default is 0.04.
;
;	HEAD:	The size of the arrow head relative to LENGTH. Default is 0.3.
;
;	NSTEPS:	The number of line_segments forming each arrow. Default is 5.
;
;	X0:	Positions of the arrows. Default is random.
;
;       DYNAMIC:The number of segments by which X0 is shifted. Default is 0.
;
;	SEED: Random number generator seed for proper randomization when VEL
;	      is called repaetedly. It has to be a NAMED VARIABLE.
;
;	The following 3 keyword parameters should be supplied together
;	when the XX,YY parameters describe a nonuniform grid. Do not initialize
;	them! Their purpose is to save on the triangulation time for the next
;	call of VEL, if the grid does not change.
;
;	XXOLD:  XX in previous call, on output XX in this call.
;
;	YYOLD:  YY in previous call, on output YY in this call.
;
;       TRIANGLES: Array to store the triangles from the triangulization.
;
;	NOERASE: NOERASE=1 or /NOERASE does not allow erase for VEL.
;
; OUTPUTS:
;	A velocity field graph is drawn on the current graphics device.
;
;	If X0 is given it returns the position for the next set of arrows.
;	If XXOLD,YYOLD,TRIANGLES are given, they return the old grid
;	coordinates and the triangles from their triangulation.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	A plot is drawn on the current graphics device.
;
; RESTRICTIONS:
;	None.
;
; PROCEDURE:
;	NVECS random points within (X,Y) are selected. For each "shot" the 
;	(U,V) field (as bilinearly interpolated) at each point is followed 
;	with NSTEPS line segments. An arrow head is drawn at the end.
;	The length of the arrow corresponding to a velocity MAXVAL is LENGTH.
;	If X0 is provided it is used for the location of the arrows instead
;	of random positions and X0 is advansed by DYNAMIC segments.
;
; MODIFICATION HISTORY:
;	12/2/92	- modified to handle !p.multi (jiy-RSI)
;	Neal Hurlburt, April, 1988.
;       7/12/94 HJM - Fixed error in weighting factors in function
;                     vel_mybi() which produced incorrect velocity vectors.
;       30/5/96 - R. Keppens changed to allow for dynamic velocity plot.
;		  Added X and Y parameters for correct axis labeling.
;                 X0 stores starting positions, DYNAMIC the index for shifts,
;                 MAXVAL the common normalization for multiple calls.
;       03/6/96 - G. Toth changed a few more things. 
;                 No TITLE keyword parameter, just use !p.title.
;                 The use of !radeg is corrected (used to be the inverse).
;                 Simplified ARRHEAD and ARROW, and merged both into VEL.
;		  Aspect ratio is taken into account for drawing arrow heads.
;		  Use NOCLIP=0 for clipping arrows outside box.
;                 Replaced vel_mybi() with IDL's interpolate() function.
;		  LENGTH is made to be relative to the diagonal of the box.
;		  HEAD parameter added for arrowhead size. 0 for streamlines.
;		  NSTEPS now means the number of segments (instead of points).
;	25/11/96- G. Toth added /NOERASE keyword
;	02/12/96- G. Toth and R. Keppens fixed the plotting error for
;		  differrent number of grid points in the X and Y directions.
;       03/27/97- G. Toth and D. Innes generalized VEL for non-uniform grids.
;	04/07/00- G. Toth added second order integration useful from accurate
;		  stream line drawing. It is used when NSTEPS>5
;-

;Return to caller if an error occurs
;
;on_error,2

sizes=size(U)
nx=sizes(1)
ny=sizes(2)
; Default values for optional parameters
;
if n_elements(XX)       eq 0 then xx      = [0,nx-1]
if n_elements(YY)       eq 0 then yy      = [0,ny-1]
if n_elements(NVECS)    eq 0 then nvecs   = 200
if not keyword_set(MAXVAL)   then maxval  = max(sqrt(u^2+v^2))
if n_elements(LENGTH)   eq 0 then length  = 0.04
if n_elements(HEAD)     eq 0 then head    = 0.3
if n_elements(NSTEPS)   eq 0 then nsteps  = 5
if n_elements(DYNAMIC)  eq 0 then dynamic = 0

; Derived parameters and derived defaults
;
xmin=min(XX)   &   ymin=min(YY)   &   xmax=max(XX)   &   ymax=max(YY)
xdel=xmax-xmin &   ydel=ymax-ymin

if dynamic lt 0      then dynamic=0
if dynamic gt nsteps then dynamic=nsteps

; Check is X0 is provided
random = n_elements(X0) ne 2*nvecs

if random and keyword_set(X0) then print,$
   'VEL: Initial position array X0 has incorrect size. Using random positions.'

if random then begin
   ; Random positions within (xmin-xmax,ymin-ymax) if x0 is not defined
   ;
   x0=fltarr(nvecs,2)
   x0(*,0)=xmin+xdel*randomu(seed,nvecs)
   x0(*,1)=ymin+ydel*randomu(seed,nvecs)
endif else begin
   ; Put the x0 arrow positions inside the box by applying a mod function
   ;
   x0(*,0)=x0(*,0)-xdel*floor((x0(*,0)-xmin)/xdel)
   x0(*,1)=x0(*,1)-ydel*floor((x0(*,1)-ymin)/ydel)
endelse

; Assume first that U and V are defined on a regular grid
irregular=0
ureg=u
vreg=v
; Check if there is a grid passed in XX, YY
if n_elements(xx) gt 2 and n_elements(yy) gt 2 then begin
   ; Check if the grid is uniform by calculating maximum grid spacing
   if min(xx(1:nx-1,*)-xx(0:nx-2,*)) lt xdel/(nx-0.99999) or $
      min(yy(*,1:ny-1)-yy(*,0:ny-2)) lt ydel/(ny-0.99999) then begin
      ; Check if the triangulation has already been done
      newx=1
      if keyword_set(triangles) and $
         n_elements(xxold) eq nx*ny and n_elements(yyold) eq nx*ny then $
            if max(abs(xx-xxold))+max(abs(yy-yyold)) eq 0 then newx=0
      if newx then begin
          print,'Triangulating in VEL ...'
          triangulate,xx,yy,triangles
          xxold=xx
          yyold=yy
      endif
      ; Interpolate velocity arrays onto an nx*ny sized REGULAR grid
      dx=xdel/(nx-1.00001) & dy=ydel/(ny-1.00001)
      ureg=trigrid(xx,yy,u,triangles,[dx,dy])
      vreg=trigrid(xx,yy,v,triangles,[dx,dy])
      irregular=1
      minval=maxval*1e-4
   endif
endif

; Calculate the coordinates forming the arrows: X(nvecs,nsteps+4,2)
; --> nvecs arrows to be drawn as 
; --> nsteps straight line segments that trace out the vector field
;     with 3 additional line segments forming the arrow head. 
;     The nsteps+3 segments are defined by altogether nsteps+4 points.
; --> third dimension is 0 and 1 for X and Y

; Starting positions
;
x=fltarr(nvecs,nsteps+4,2)
x(*,0,*)=x0

; The time step based on the length of the longest arrow and the maxval speed
;
dt=length*sqrt(xdel^2+ydel^2)/nsteps/maxval

; Integrate the velocity field for nsteps
;
for i=1,nsteps do begin
   xt=(nx-1)*(x(*,i-1,0)-xmin)/xdel
   yt=(ny-1)*(x(*,i-1,1)-ymin)/ydel
   ut=interpolate(ureg,xt,yt)
   vt=interpolate(vreg,xt,yt)
   if irregular and i eq 1 then for ivec=0,nvecs-1 do begin
      if random and abs(ut(ivec)) lt minval and abs(vt(ivec)) lt minval then $
      begin
         xivec=randomu(seed)         & yivec=randomu(seed)
         x(ivec,0,0)=xmin+xdel*xivec & x(ivec,0,1)=ymin+ydel*yivec
         xivec=xivec*(nx-1)          & yivec=yivec*(ny-1)
         ut(ivec)=interpolate(ureg,xivec,yivec)
         vt(ivec)=interpolate(vreg,xivec,yivec)
      endif
   endfor

   ; Second order for stream line integration
   if nsteps gt 5 then begin
      xt1=xt+0.5*dt*ut*(nx-1)/xdel
      yt1=yt+0.5*dt*vt*(ny-1)/ydel
      ut=interpolate(ureg,xt1,yt1)
      vt=interpolate(vreg,xt1,yt1)
   endif

   x(*,i,0)=x(*,i-1,0)+ut*dt
   x(*,i,1)=x(*,i-1,1)+vt*dt
endfor

; Put the position of the arrowhead into the elements nsteps+1..nsteps+3.
; The arrow head consists of two "wings" starting from the tip of the arrow.
; The wings are at a 30 degree angle (taking into account the aspect ratio).
; The projected length of the wings onto the arrow is head*nsteps*last_segment.
; The last segment is between x(*,nsteps,*) and x(*,nsteps-1,*).

aspect = xdel/ydel
tant   = tan(30/!radeg)
scal   = head*nsteps

i = nsteps   &   h = nsteps-1

dx = scal*(x(*,i,0)-x(*,h,0))	   ; Vector of the projected wings
dy = scal*(x(*,i,1)-x(*,h,1))
xm = x(*,i,0)-dx                   ; The coordinates of the midpoint
ym = x(*,i,1)-dy                   ; between the two 'wings'
dx = dx*tant/aspect
dy = dy*tant*aspect                ; Contract by tan(theta) and use aspect.

x(*,i+1,0) = xm-dy                 ; (dx,dy) is rotated to (-dy,dx) and added
x(*,i+2,0) = x(*,i,0)		   ; to and subtracted from (xm,ym).
x(*,i+3,0) = xm+dy

x(*,i+1,1) = ym+dx
x(*,i+2,1) = x(*,i,1)
x(*,i+3,1) = ym-dx

; Put new arrow positions into X0 for next call based on dynamic [0..nsteps]

x0(*,*) = x(*,dynamic,*)

; Draw box
plot,[xmin,xmax],[ymin,ymax],/nodata,XSTYLE=1,YSTYLE=1,NOERASE=noerase

; Draw arrows
for i=0,nvecs-1 do plots,x(i,*,0),x(i,*,1),NOCLIP=0

return
end
;=============================================================================

