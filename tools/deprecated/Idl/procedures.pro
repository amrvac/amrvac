; Written by G. Toth for the Versatile Advection Code
;
; Procedures for 
; 
; reading ascii and binary data produced by VAC and VACINI:
;    openfile,gettype,gethead,getpict,getpict_asc,getpict_bin
; reading numbers and strings from input:
;    asknum, askstr, str2arr, readplotpar, readlimits
; transforming initial data:
;    regulargrid, polargrid, unpolargrid, spheregrid, getaxes
; calculating functions of the data
;    getfunc, getlimits
; plotting
;    plotfunc, plotgrid
; calculating cell corners and cell volumes for general 2D grids
;    gengrid
; comparing two w or x arrays for relative differences
;    compare
; checking index ranges for functions quadruplet and triplet
;    checkdim
; procedure "quit" as an alias for "exit"
;    quit
;
; Functions for
;    
; calculating derivatives in 2D for Cartesian grids to 2nd,3rd,4th order
;    diff2,diff3,diff4
; calculating derivatives in 2D for general grids
;    grad,div,curl,grad_rz,div_rz,curl_rz, filledge,intedge,intedge_rz
; taking a part of an array or coarsen an array
;    triplet, quadruplet, coarse

;==========================================
pro openfile,unit,filename,filetype
;==========================================
   on_error,2

   close,unit
   case filetype of
       'ascii':	openr,unit,filename
       'binary':openr,unit,filename,/f77_unf
       'real4' :openr,unit,filename,/f77_unf
       else: 	print,'Openfile: unknown filetype:',filetype
   endcase
end

;==========================================
pro gettype,filenames,filetypes,npictinfiles
;==========================================
   on_error,2

   filetypes=filenames
   npictinfiles=intarr(n_elements(filenames))
   for ifile=0,n_elements(filenames)-1 do begin
      ; Obtain filetype based on the length info in the first 4 bytes
      openr,10,filenames(ifile)
      len=long(1)
      readu,10,len
      if len ne 79 then ftype='ascii' else begin
         ; The length of the 2nd line decides between real4 and binary
         ; since it contains the time, which is real*8 or real*4
         head=bytarr(79+4)
         readu,10,head,len
         case len of
            20: ftype='real4'
            24: ftype='binary'
            else: begin
               print,'Error in GetType: strange unformatted file:',$
                     filenames(ifile)
               retall
            end
         endcase
      endelse
      close,10

      ; Obtain file size
      openfile,1,filenames(ifile),ftype
      status=fstat(1)
      fsize=status.size

      ; Obtain size of a single snapshot
      gethead,1,ftype,phys, $
         headline,it,time,gencoord,ndim,neqpar,nw,nx
      nxs=1
      for idim=1,ndim do nxs=nxs*nx(idim-1)
      case ftype of
       'ascii': pictsize=1+79 + 1+7+13+9 + 1+ndim*4 + 1+neqpar*13 + 1+79 + $
                  (18*(ndim+nw)+1)*nxs
       'binary':pictsize=8+79 + 8+4*4+8  + 8+ndim*4 + 8+neqpar*8  + 8+79 + $
		  8*(1+nw)+8*(ndim+nw)*nxs
       'real4':pictsize=8+79 + 5*4+8  + 8+ndim*4 + 8+neqpar*4  + 8+79 + $
		  8*(1+nw)+4*(ndim+nw)*nxs
      endcase
      close,1

      ; Calculate number of snapshots in the file
      npictinfiles(ifile)= fsize/pictsize
      filetypes(ifile)=ftype
   endfor
end

;==========================================
pro gethead,unit,filetype,physics,headline,it,time,gencoord, $
            ndim,neqpar,nw,nx,eqpar,variables
;==========================================
   on_error,2

;Type definitions
   headline='                                                                               '
   it=long(1)
   ndim=long(1)
   neqpar=long(1)
   nw=long(1)
   varname='                                                                               '
;Read header
   case filetype of
      'ascii': begin
                  time=double(1)
		  readf,unit,headline
		  readf,unit,it,time,ndim,neqpar,nw
		  gencoord=(ndim lt 0)
		  ndim=abs(ndim)
		  nx=lonarr(ndim)
		  readf,unit,nx
		  eqpar=dblarr(neqpar)
		  readf,unit,eqpar
		  readf,unit,varname
	       end
      'binary':begin
                  time=double(1)
		  readu,unit,headline
		  readu,unit,it,time,ndim,neqpar,nw
		  gencoord=(ndim lt 0)
		  ndim=abs(ndim)
		  nx=lonarr(ndim)
		  readu,unit,nx
		  eqpar=dblarr(neqpar)
		  readu,unit,eqpar
		  readu,unit,varname
               end
      'real4': begin
                  time=float(1)
		  readu,unit,headline
		  readu,unit,it,time,ndim,neqpar,nw
		  gencoord=(ndim lt 0)
		  ndim=abs(ndim)
		  nx=lonarr(ndim)
		  readu,unit,nx
		  eqpar=fltarr(neqpar)
		  readu,unit,eqpar
		  readu,unit,varname
	       end
      else: print,'Gethead: unknown filetype',filetype
   endcase
   variables=str_sep(strtrim(strcompress(varname),2),' ')
   tmp=str_sep(strtrim(headline,2),'_')
   if n_elements(tmp) eq 2 then begin
      headline=tmp(0)
      physics=tmp(1)
   endif
end

;==========================================
pro getpict,unit,filetype,npict,$
    x,w,headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,error
;==========================================
   on_error,2

   ipict=0
   while ipict lt npict and not eof(unit) do begin
      ipict=ipict+1
      case filetype of
         'ascii': getpict_asc,unit,npict,$
             x,w,headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables
         'binary': getpict_bin,unit,npict,$
             x,w,headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables
         'real4': getpict_real,unit,npict,$
             x,w,headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables
         else:     begin
                      print,'Getpict: unknown filetype:',filetype
		      ipict=0
                      close,unit
                   end
      endcase
   endwhile
   error=(ipict LT npict)
end

;==========================================
pro getpict_asc,unit,npict,$
   x,w,headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables
;==========================================
  on_error,2                      
  ;----------------------------------------
  ;Read header information
  ;----------------------------------------
  gethead,unit,'ascii',physics,$
     headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables
  ndim=abs(ndim)
  ;----------------------------------------
  ; Read coordinates and values row by row
  ;----------------------------------------
  wrow=dblarr(nw)
  xrow=dblarr(ndim)
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
    x=dblarr(nx(0),ndim)
    w=dblarr(nx(0),nw)
    for x0=0,nx(0)-1 do begin
      readf,unit,xrow,wrow
      x(x0,0:ndim-1)=xrow(0:ndim-1)
      w(x0,0:nw-1)  =wrow(0:nw-1)
    endfor
  end
  ;-------------- 2D ----------------------
  2: begin
    x=dblarr(nx(0),nx(1),ndim)
    w=dblarr(nx(0),nx(1),nw)
    for x1=0,nx(1)-1 do begin
      for x0=0,nx(0)-1 do begin
        readf,unit,xrow,wrow
        x(x0,x1,0:ndim-1)=xrow(0:ndim-1)
        w(x0,x1,0:nw-1)  =wrow(0:nw-1)
      endfor
    endfor
  end
  ;-------------- 3D ----------------------
  3: begin
    x=dblarr(nx(0),nx(1),nx(2),ndim)
    w=dblarr(nx(0),nx(1),nx(2),nw)
    for x2=0,nx(2)-1 do begin
      for x1=0,nx(1)-1 do begin
        for x0=0,nx(0)-1 do begin
          readf,unit,xrow,wrow
          x(x0,x1,x2,0:ndim-1)=xrow(0:ndim-1)
          w(x0,x1,x2,0:nw-1)=wrow(0:nw-1)
	endfor
      endfor
    endfor
  end
  endcase
end

;==========================================
pro getpict_bin,unit,npict,$
   x,w,headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables
;==========================================
  on_error,2                      
  ;----------------------------------------
  ;Read header information
  ;----------------------------------------
  gethead,unit,'binary',physics,$
     headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables
  ;----------------------------------------
  ; Read coordinates and values 
  ;----------------------------------------
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
    n1=nx(0)
    x=dblarr(n1,ndim)
    w=dblarr(n1,nw)
    wi=dblarr(n1)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,iw)=wi
    endfor
  end
  ;-------------- 2D ----------------------
  2: begin
    n1=nx(0)
    n2=nx(1)
    x=dblarr(n1,n2,ndim)
    w=dblarr(n1,n2,nw)
    wi=dblarr(n1,n2)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,*,iw)=wi
    endfor
  end
  ;-------------- 3D ----------------------
  3: begin
    n1=nx(0)
    n2=nx(1)
    n3=nx(2)
    x=dblarr(n1,n2,n3,ndim)
    w=dblarr(n1,n2,n3,nw)
    wi=dblarr(n1,n2,n3)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,*,*,iw)=wi
    endfor
  end
  endcase
end

;==========================================
pro getpict_real,unit,npict,$
   x,w,headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables
;==========================================
  on_error,2                      
  ;----------------------------------------
  ;Read header information
  ;----------------------------------------
  gethead,unit,'real4',physics,$
     headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables
  ;----------------------------------------
  ; Read coordinates and values 
  ;----------------------------------------
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
    n1=nx(0)
    x=fltarr(n1,ndim)
    w=fltarr(n1,nw)
    wi=fltarr(n1)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,iw)=wi
    endfor
  end
  ;-------------- 2D ----------------------
  2: begin
    n1=nx(0)
    n2=nx(1)
    x=fltarr(n1,n2,ndim)
    w=fltarr(n1,n2,nw)
    readu,unit,x
    wi=fltarr(n1,n2)
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,*,iw)=wi
    endfor
  end
  ;-------------- 3D ----------------------
  3: begin
    n1=nx(0)
    n2=nx(1)
    n3=nx(2)
    x=fltarr(n1,n2,n3,ndim)
    w=fltarr(n1,n2,n3,nw)
    wi=fltarr(n1,n2,n3)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,*,*,iw)=wi
    endfor
  end
  endcase
end

;==========================================
pro asknum,prompt,var,doask
;==========================================
   on_error,2

   if var eq 0 then read,PROMPT=prompt+'? ',var $
   else begin
      if doask then begin
         tmp=''
         read,PROMPT=prompt+'='+strtrim(string(var),2)+' ? ',tmp
         if tmp ne '' then reads,tmp,var
      endif else print,prompt,'=',var
   endelse
end

;==========================================
pro askstr,prompt,var,doask
;==========================================
   on_error,2                      

   if var eq '' then read,PROMPT=prompt+'? ',var $
   else begin
      if doask then begin
         tmp=''
         read,PROMPT=prompt+'='+var+' ? ',tmp
         if tmp ne '' then var=tmp
      endif else print,prompt,'=',var
   endelse
end

;==========================================
pro str2arr,s,a,n,sep
;==========================================
; Split string s at the sep characters (default is space) into array a
; If n is 0, it will be the size of the array on output
; If n is defined, fill up the rest of the array with the last defined element
; If n is defined but smaller than the number of elements in s, print a warning

on_error,2

if keyword_set(sep) then $
   a0=str_sep(s,sep)     $
else                     $
   a0=str_sep(strtrim(strcompress(s),2),' ')

n0=n_elements(a0)

if not keyword_set(n) then begin
   a=a0
   n=n0
endif else if n ge n0 then begin
   a=strarr(n)
   a(0:n0-1)=a0
   if n0 lt n then a(n0:n-1)=a0(n0-1)
endif else begin
   a=strarr(n)
   a=a0(0:n-1)
   print,'Warning: more than',n,' values defined by string: ',s,$
       FORMAT='(a,i3,a,a)'
endelse

end

;===========================================================================
pro readplotpar,ndim,cut,cut0,plotdim,nfunc,func,funcs,funcs1,funcs2,$
   plotmode,plotmodes,plottitle,plottitles,autorange,autoranges,doask
;===========================================================================
   on_error,2

   ; Determine dimension of plots based on cut or ndim, 
   ; calculate reduced cut0 array by eliminating degenerate dimensions
   siz=size(cut)
   if siz(0) eq 0 then begin
      plotdim=ndim
      cut0=0
   endif else begin
      plotdim=0
      for i=1,siz(0) do if siz(i) gt 1 then plotdim=plotdim+1
      cut0=reform(cut)
      if siz(0) eq 3 and siz(2) eq 1 then cut0=reform(cut,siz(1),siz(3))
   endelse

   askstr,'func(s) (e.g. rho p m1;m2 r*m1 -T) ',func,doask
   if plotdim eq 1 then begin
      print,'1D plotmode: plot'
      plotmode='plot' 
   endif else begin
      if plotmode eq 'plot' then plotmode=''
     print,'2D plotmodes: contour/contlabel/contfill/shade/shadeirr/surface/tv/'
      print,'              stream/vel/velovect/ovelovect'
      askstr,'plotmode(s)                ',plotmode,doask
   endelse
   askstr,'plottitle(s) (e.g. B [G];J)',plottitle,doask
   askstr,'autorange(s) (y/n)         ',autorange,doask

   nfunc=0
   str2arr,func,funcs,nfunc
   str2arr,plotmode,plotmodes,nfunc
   str2arr,plottitle,plottitles,nfunc,';'
   str2arr,autorange,autoranges,nfunc

   funcs1=strarr(nfunc)
   funcs2=strarr(nfunc)
   for ifunc=0,nfunc-1 do begin
      func12=str_sep(funcs(ifunc),';')
      funcs1(ifunc)=func12(0)
      if n_elements(func12) eq 2 then funcs2(ifunc)=func12(1)
   endfor

end
;===========================================================================
pro readtransform,ndim,nx,gencoord,transform,nxreg,wregpad,$
    physics,nvector,vectors,grid,doask
;===========================================================================
   on_error,2

   if (gencoord or transform eq 'unpolar') and ndim eq 2 then begin
      if transform eq '' then begin
         transform='none'
         askstr,"transform (regular/polar/unpolar/none)",transform,1
      endif else $
         askstr,"transform (regular/polar/unpolar/none)",transform,doask
      case 1 of
        transform eq 'regular':begin
           print,'Generalized coordinates, dimensions for regular grid'
           nxreg0=nxreg(0)
           nxreg1=nxreg(1)
           asknum,'nxreg(0)',nxreg0,doask
           asknum,'nxreg(1)',nxreg1,doask
           grid=lindgen(nxreg0,nxreg1)
           nxreg=[nxreg0,nxreg1]
           wregpad=0
	end
	transform eq 'polar' or transform eq 'unpolar':begin
	   getvectors,physics,nvector,vectors
	   grid=lindgen(nx(0),nx(1))
	end
	transform eq 'none':grid=lindgen(nx(0),nx(1))
	else: print,'Unknown value for transform:',transform
      endcase
   endif else if gencoord and ndim eq 3 then begin
      if transform eq '' then begin 
         transform="none" & askstr,"transform (sphere/none)",transform,1
      endif else $
         askstr,"transform (sphere/none)",transform,doask
      if transform eq 'sphere' then getvectors,physics,nvector,vectors
      grid=lindgen(nx(0),nx(1),nx(2))
   endif else case ndim of
      1: grid=lindgen(nx(0))
      2: grid=lindgen(nx(0),nx(1))
      3: grid=lindgen(nx(0),nx(1),nx(2))
   endcase

   ;===== GRID HELPS TO CREATE A CUT, E.G.: cut=grid(*,4)

   help,grid
end

;===========================================================================
pro getvectors,physics,nvector,vectors
;===========================================================================
   physic=strtrim(physics)
   phys=strmid(physic,0,strlen(physic)-2)
   ndir=0
   reads,strmid(physic,strlen(physic)-1,1),ndir
   case phys of
   'rho':nvector=0
   'hd' :begin
	 nvector=1
	 vectors=1
	 end
   'hdadiab':begin
	 nvector=1
	 vectors=1
	 end
   'mhdiso':begin
	 nvector=2
	 vectors=[1,ndir+1]
	 end
   'mhd':begin
	 nvector=2
	 vectors=[1,ndir+2]
	 end
   else:begin
      if nvector eq 0 then begin
	 print,'Unrecognised physics: ',physics
	 print,'Vector variables to transform for WREG'
	 asknum,'nvector',nvector,doask
	 if nvector gt 0 then begin
	    vectors=intarr(nvector)
	    read,PROMPT='Indices of first components in w? ',vectors
	 endif
      endif
      end
   endcase
end

;===========================================================================
pro readlimits,nfunc,funcs,autoranges,noautorange,fmax,fmin,doask
;===========================================================================
   on_error,2

   if n_elements(fmax) ne nfunc then fmax=dblarr(nfunc)
   if n_elements(fmin) ne nfunc then fmin=dblarr(nfunc)
   ; check if there is any function for which autorange is 'y'
   noautorange=1
   for ifunc=0,nfunc-1 do noautorange=noautorange and autoranges(ifunc) eq 'n'

   if(noautorange)then begin
      for ifunc=0,nfunc-1 do begin
         f_min=fmin(ifunc)
         f_max=fmax(ifunc)
         asknum,'Min value for '+funcs(ifunc),f_min,doask
         asknum,'Max value for '+funcs(ifunc),f_max,doask
         fmin(ifunc)=f_min
         fmax(ifunc)=f_max
      endfor
   endif

end
;===========================================================================
pro getlimits,first,nfunc,funcs,funcs1,funcs2,autoranges,fmax,fmin,doask,$
                x,w,xreg,wreg,usereg,physics,eqpar,wnames,cut
;===========================================================================
   on_error,2                      

   for ifunc=0,nfunc-1 do begin
      if autoranges(ifunc) eq 'n' then begin
         if first then begin
            f_min=fmin(ifunc)
            f_max=fmax(ifunc)
            asknum,'Min value for '+funcs(ifunc),f_min,doask
            asknum,'Max value for '+funcs(ifunc),f_max,doask
            fmin(ifunc)=f_min
            fmax(ifunc)=f_max
         endif
      endif else begin
         if usereg then getfunc,f,f1,f2,funcs1(ifunc),funcs2(ifunc),   $
                            xreg,wreg,physics,eqpar,wnames,cut $
         else           getfunc,f,f1,f2,funcs1(ifunc),funcs2(ifunc),   $
                            x,   w,   physics,eqpar,wnames,cut

         f_max=max(f)
         f_min=min(f)
         if first then begin
            fmax(ifunc)=f_max
            fmin(ifunc)=f_min
         endif else begin
            if f_max gt fmax(ifunc) then fmax(ifunc)=f_max
            if f_min lt fmin(ifunc) then fmin(ifunc)=f_min
         endelse
      endelse
   endfor
end

;===========================================================================
pro regulargrid,x_old,nxreg_old,x,xreg,nxreg,dxreg,w,wreg,nw,wregpad,triangles
;
;    Regularize grid and interpolate "w" via triangulate and trigrid.
;    The original "w" data is interpolated into "wreg", for points outside
;    the convex hull of the irregular grid the "wregpad(nw)" array is used.
;
;    If "x_old" and "x" or "nxreg_old" and "nxreg" are different 
;    a triangulization is done first and a regular coordinate array 
;    "xreg" is created. The size of the "xreg" array is in "nxreg(2)", 
;    "dxreg(2)" is the spacing. The triangles are saved in "triangles".
; 
;    "q" can be interpolated from the irregular grid to the regular one by:
;
;    qreg(*,*)=trigrid(x(*,*,0),x(*,*,1),q,triangles,dxreg)
;
;===========================================================================
   on_error,2                      

   ;Floating underflow is not a real error, the message is suppressed
   err=check_math(1,1)

   xx=x(*,*,0)
   yy=x(*,*,1)
   ; Check if nxreg==nxreg_old and x==x_old
   newx=1
   if n_elements(nxreg_old) eq n_elements(nxreg) then $
      if max(abs(nxreg_old-nxreg)) eq 0 then $
         if n_elements(x_old) eq n_elements(x) then $
            if max(abs(x_old-x)) eq 0 then newx=0

   if newx then begin
      x_old=x
      nxreg_old=nxreg
      dxreg=[(max(xx)-min(xx))/(nxreg(0)-1.00001), $
             (max(yy)-min(yy))/(nxreg(1)-1.00001)]

      xreg=dblarr(nxreg(0),nxreg(1),2)
      for i=0,nxreg(1)-1 do xreg(*,i,0)=dxreg(0)*indgen(nxreg(0))+min(xx)
      for i=0,nxreg(0)-1 do xreg(i,*,1)=dxreg(1)*indgen(nxreg(1))+min(yy)

      triangulate,xx,yy,triangles
      wreg=dblarr(nxreg(0),nxreg(1),nw)
  endif
  if not keyword_set(wregpad) then begin
     wregpad=dblarr(nw)
     for iw=0,nw-1 do begin
        wmax=max(w(*,*,iw))
        wmin=min(w(*,*,iw))
        if wmax*wmin lt 0 then wregpad(iw)=0 $
        else                   wregpad(iw)=wmin-0.1*(wmax-wmin)
     endfor
  endif
  for iw=0,nw-1 do $
     wreg(*,*,iw)=trigrid(xx,yy,w(*,*,iw),triangles,dxreg,missing=wregpad(iw))

   err=check_math(0,0)
   ;Floating underflow is not a real error, the message is suppressed
   if err ne 32 and err ne 0 then print,'Math error in regulargrid:',err

end

;===========================================================================
pro polargrid,nvector,vectors,x,w,xreg,wreg
;
;    Transform vector variables from x and y to radial and phi components
;
;===========================================================================
  on_error,2                      

  xreg=x
  xreg(*,*,0)=sqrt(x(*,*,0)^2+x(*,*,1)^2)
  xreg(*,*,1)=atan(x(*,*,1),x(*,*,0))
  phi=xreg(*,*,1)
  wreg=w
  for i=1,nvector do begin
     ivx=vectors(i-1)
     ivy=ivx+1
     wreg(*,*,ivx)=  w(*,*,ivx)*cos(phi)+w(*,*,ivy)*sin(phi)
     wreg(*,*,ivy)= -w(*,*,ivx)*sin(phi)+w(*,*,ivy)*cos(phi)
  endfor

  ;Remove 2*pi jumps from phi
  pi2=8*atan(1) & sz=size(phi) & nx2=sz(2)
  for ix2=1,nx2-1 do while phi(1,ix2-1) gt phi(1,ix2) do $
     phi(*,ix2)=phi(*,ix2)+pi2

  xreg(*,*,1)=phi
end

;===========================================================================
pro spheregrid,nvector,vectors,x,w,xreg,wreg
;
;    Transform vector variables from x,y,z to radial,phi,z components
;
;===========================================================================
  on_error,2                      

  xreg=x
  xreg(*,*,*,0)=sqrt(x(*,*,*,0)^2+x(*,*,*,1)^2+x(*,*,*,2)^2)
  xreg(*,*,*,2)=-atan(x(*,*,*,2),x(*,*,*,0))
  xreg(*,*,*,1)=atan(x(*,*,*,1),sqrt(x(*,*,*,0)^2+x(*,*,*,2)^2))
  phi=xreg(*,*,*,2)
  theta=xreg(*,*,*,1)
  wreg=w
  sinphi=sin(phi)
  cosphi=cos(phi)
  sintheta=sin(theta)
  costheta=cos(theta)
  for i=1,nvector do begin
     ivx=vectors(i-1)
     ivy=ivx+1
     ivz=ivy+1
     wreg(*,*,*,ivx)=(w(*,*,*,ivx)*cosphi-w(*,*,*,ivz)*sinphi)*costheta $
                     +w(*,*,*,ivy)*sintheta
     wreg(*,*,*,ivz)=-w(*,*,*,ivx)*sinphi-w(*,*,*,ivz)*cosphi
     wreg(*,*,*,ivy)=(-w(*,*,*,ivx)*cosphi+w(*,*,*,ivz)*sinphi)*sintheta $
                     +w(*,*,*,ivy)*costheta
  endfor

  ;Remove 2*pi jumps from phi
  pi=4*atan(1) & pi2=2*pi & sz=size(phi) & nx2=sz(2) & nx3=sz(3)
  for ix3=1,nx3-1 do while phi(1,1,ix3-1) gt phi(1,1,ix3) do $
     phi(*,*,ix3)=phi(*,*,ix3)+pi2

  ;Remove turn over from theta
  for ix2=1,nx2-1 do $
  if theta(1,ix2-1,1) ge theta(1,ix2,1) then begin
     if theta(1,ix2,1) lt 0 then $
          theta(*,ix2-1,*)=-pi-theta(*,ix2-1,*) $
     else $
          theta(*,ix2,*)=pi-theta(*,ix2,*)
  endif

  xreg(*,*,*,2)=phi
  xreg(*,*,*,1)=theta
end

;===========================================================================
pro unpolargrid,nvector,vectors,x,w,xreg,wreg
;
;    Transform vector variables from x and y to radial and phi components
;
;===========================================================================
  on_error,2                      

  xreg=x
  phi=x(*,*,1)
  xreg(*,*,0)=x(*,*,0)*cos(phi)
  xreg(*,*,1)=x(*,*,0)*sin(phi)

  wreg=w
  for i=1,nvector do begin
     ivx=vectors(i-1)
     ivy=ivx+1
     wreg(*,*,ivx)=  w(*,*,ivx)*cos(phi)-w(*,*,ivy)*sin(phi)
     wreg(*,*,ivy)=  w(*,*,ivx)*sin(phi)+w(*,*,ivy)*cos(phi)
  endfor
end

;===========================================================================
pro getaxes,ndim,x,xx,yy,zz,cut,cut0,plotdim,variables
;===========================================================================
on_error,2
case ndim of
  1: xx=x
  2: begin
        xx=x(*,*,0)
        yy=x(*,*,1)
     end
  3: begin
       xx=x(*,*,*,0)
       yy=x(*,*,*,1)
       zz=x(*,*,*,2)
     end
endcase

if keyword_set(cut0) then begin
   xx=xx(cut0)
   if ndim gt 1 then yy=yy(cut0)
   if ndim gt 2 then zz=zz(cut0)
endif

!x.title=variables(0)
if plotdim gt 1 then !y.title=variables(1)
if plotdim gt 2 then !z.title=variables(2)

; Cut with fixed X value?
siz=size(cut)
; in 2D
if siz(0) eq 2 and siz(1) eq 1 then begin
   xx=yy
   !x.title=variables(1)
endif
; in 3D
if siz(0) eq 3 then begin
   case 1 of
   plotdim eq 1: begin
      xx=zz
      !x.title=variables(2)
   end
   siz(1) eq 1: begin
         xx=yy
         yy=zz
         !x.title=variables(1)
         !y.title=variables(2)
   end
   siz(2) eq 1: begin
      yy=zz
      !y.title=variables(2)
   end
   else: print,'internal error in getaxes'
   endcase
endif

end

;===========================================================================
pro getfunc,f,f1,f2,func1,func2,x,w,physics,eqpar,wnames,cut
;===========================================================================
on_error,2

f1=animfunc(x,w,func1,physics,eqpar,wnames)

if keyword_set(cut) then f1=f1(cut)

if func2 eq '' then f=f1 else begin

   f2=animfunc(x,w,func2,physics,eqpar,wnames)

   if keyword_set(cut) then f2=f2(cut)

   ; Calculate f=sqrt(f1^2+f2^2)
   f=sqrt(f1^2+f2^2)
endelse

end

;===========================================================================
pro plotfunc,x,w,xreg,wreg,usereg,ndim,physics,eqpar,$
  variables,wnames,axistype,plotmodes,plottitles,$
  ax,az,contourlevel,linestyle,$
  velvector,velspeed,velseed,velpos,velx,vely,veltri,$
  cut,cut0,plotdim,$
  nfunc,multix,multiy,plotix,plotiy,funcs,funcs1,funcs2,fmin,fmax,f
;===========================================================================
   on_error,2

   plotsizex=!d.x_size/multix
   plotsizey=!d.y_size/multiy
   tvsizex=plotsizex-2*!d.x_ch_size
   tvsizey=plotsizey-4*!d.y_ch_size

   if axistype eq 'coord' then begin
      if usereg then getaxes,ndim,xreg,xx,yy,zz,cut,cut0,plotdim,variables $
      else           getaxes,ndim,x   ,xx,yy,zz,cut,cut0,plotdim,variables
   endif

   for ifunc=0,nfunc-1 do begin
      !p.title=plottitles(ifunc)
      if !p.title eq 'default' then !p.title=funcs(ifunc)

      plotmod=plotmodes(ifunc)

      ; Calculate the next p.multi(0) explicitly
      if !p.multi(0) gt 0 then multi0=!p.multi(0)-1 $
      else multi0=!p.multi(1)*!p.multi(2)-1

      ; Calculate subplot position indices
      if !p.multi(4) then begin
         plotix=multix-1-multi0/multiy
         plotiy=multi0-multiy*(multi0/multiy)
      endif else begin
         plotix=multix*(multi0/multix+1)-multi0-1
         plotiy=multi0/multix
      endelse

      if usereg then getfunc,f,f1,f2,funcs1(ifunc),funcs2(ifunc),   $
                             xreg,wreg,physics,eqpar,wnames,cut0 $
      else           getfunc,f,f1,f2,funcs1(ifunc),funcs2(ifunc),   $
                             x,  w,   physics,eqpar,wnames,cut0
      f_min=fmin(ifunc)
      f_max=fmax(ifunc)
      if f_max eq f_min then begin
         f_max=f_max+1
         f_min=f_min-1
      endif

      if nfunc gt multix*multiy then linestyle=ifunc $
                                else linestyle=!p.linestyle

      if strmid(plotmod,0,4) eq 'cont' then $
         levels=findgen(contourlevel)/contourlevel*(f_max-f_min)+f_min

      if plotmod eq 'tv' then begin
         ; Calculate plotting position and print out title for TV plotmode
         plotx=plotsizex*plotix
         ploty=plotsizey*plotiy
         xyouts,plotx+plotsizex/2,ploty+plotsizey-1.5*!d.y_ch_size,$
            !p.title,/DEV,ALIGNMENT=0.5
         plotx=plotx+5
	 ploty=ploty+20
         if !d.name eq 'PS' then f=congrid(f,200,200) $
         else                    f=congrid(f,tvsizex,tvsizey)
         f=bytscl(f,MIN=f_min,MAX=f_max,TOP=!D.TABLE_SIZE)
      endif

      case axistype of
      'cells': case plotmod of
	 'contour': contour,f>f_min,LEVELS=levels,XSTYLE=1,YSTYLE=1,/NOERASE
	 'contlabel': contour,f>f_min,LEVELS=levels,/FOLLOW,$
						  XSTYLE=1,YSTYLE=1,/NOERASE
	 'contfill':contour,f>f_min,LEVELS=levels,/FILL,$
						  XSTYLE=1,YSTYLE=1,/NOERASE
	 'plot'     :plot,f,YRANGE=[f_min,f_max],XSTYLE=18,ystyle=18, $
						 LINE=linestyle,/NOERASE
	 'shade'    :shade_surf,f>f_min,ZRANGE=[f_min,f_max],$
			XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
         'surface'  :surface,f>f_min,ZRANGE=[f_min,f_max],$
			XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
         'tv'       :tv,f,plotx,ploty,XSIZE=tvsizex,YSIZE=tvsizey
         'vel'      :vel,f1,f2,NVECS=velvector,MAXVAL=f_max,$
			DYNAMIC=velspeed,SEED=velseed,X0=velpos,/NOERASE
	 'stream'   :begin
			eps=1.e-30
			v1=f1/sqrt(f1^2+f2^2+eps) & v2=f2/sqrt(f1^2+f2^2+eps)

			vel,v1,v2,NVECS=velvector,MAXVAL=1.,$
			NSTEP=6,LENGTH=0.06,HEAD=0.1,$
			DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE

			vel,v1,v2,NVECS=velvector,MAXVAL=1.,$
	                NSTEP=100,LENGTH=1.,HEAD=0.,$
			DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE

			v1=-v1 & v2=-v2
			vel,v1,v2,NVECS=velvector,MAXVAL=1.,$
	                NSTEP=100,LENGTH=1.,HEAD=0.,$
			DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
		    end
         'velovect' :velovect,f1,f2,/NOERASE
         'ovelovect':velovect,f1,f2,/NOERASE,$
	    XRANGE=[0,n_elements(f1(*,0))-1],YRANGE=[0,n_elements(f1(0,*))-1]
         endcase
      'coord': case plotmod of
	 'contour'  :contour,f>f_min,xx,yy,$
			LEVELS=levels,XSTYLE=1,YSTYLE=1,/NOERASE
	 'contlabel':contour,f>f_min,xx,yy,$
			LEVELS=levels,/FOLLOW,XSTYLE=1,YSTYLE=1,/NOERASE
	 'contfill' :contour,f>f_min,xx,yy, $
			LEVELS=levels,/FILL,XSTYLE=1,YSTYLE=1,/NOERASE
	 'plot'    :plot,xx,f,YRANGE=[f_min,f_max],XSTYLE=18,YSTYLE=18,$
						   LINE=linestyle,/NOERASE
	 'shade'    :shade_surf,f>f_min,xx,yy,ZRANGE=[f_min,f_max],$
			XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
	 'shadeirr' :begin
			shade_surf_irr,f>f_min,xx,yy,AX=ax,AZ=az
			shade_surf,f>f_min,xx,yy,AX=ax,AZ=az,/NODATA,/NOERASE
		     end
         'surface'  :surface,f>f_min,xx,yy,ZRANGE=[f_min,f_max],$
			XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
         'tv'       :tv,f,plotx,ploty,XSIZE=tvsizex,YSIZE=tvsizey
         'vel'      :vel,f1,f2,xx,yy,XXOLD=velx,YYOLD=vely,TRIANGLES=veltri,$
                        NVECS=velvector,MAXVAL=f_max,$
                        DYNAMIC=velspeed,SEED=velseed,X0=velpos,/NOERASE
	 'stream'   :begin
			eps=1.e-30
			v1=f1/sqrt(f1^2+f2^2+eps) & v2=f2/sqrt(f1^2+f2^2+eps)

			vel,v1,v2,xx,yy,NVECS=velvector,MAXVAL=1.,$
			XXOLD=velx,YYOLD=vely,TRIANGLES=veltri,$
			NSTEP=6,LENGTH=0.06,HEAD=0.1,$
			DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE

			vel,v1,v2,xx,yy,NVECS=velvector,MAXVAL=1.,$
			XXOLD=velx,YYOLD=vely,TRIANGLES=veltri,$
	                NSTEP=100,LENGTH=1.,HEAD=0.,$
			DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE

			v1=-v1 & v2=-v2
			vel,v1,v2,xx,yy,NVECS=velvector,MAXVAL=1.,$
			XXOLD=velx,YYOLD=vely,TRIANGLES=veltri,$
	                NSTEP=100,LENGTH=1.,HEAD=0.,$
			DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
		    end
         'velovect' :velovect,f1,f2,xx(*,0),yy(0,*),/NOERASE
         'ovelovect':velovect,f1,f2,xx(*,0),yy(0,*),/NOERASE,$
			XRANGE=[min(xx),max(xx)],YRANGE=[min(yy),max(yy)]
         endcase
      else:print,'Unknown axistype:',axistype
      endcase
      !p.multi(0)=multi0
   endfor
end
;===========================================================================
pro putbottom,multix,multiy,ix,iy,ninfo,nx,it,time

on_error,2

if ninfo lt 1 then return
info=''
if ninfo gt 2 then info='nx='+string(nx,format='(3(i4,","))')+' '
if ninfo gt 1 then info=info+'it='+string(it,format='(i6)')+', '
info=info+'time='+string(time,format='(g12.5)')
xyouts,5+(ix*!d.x_size)/multix,8+(iy*!d.y_size)/multiy,/DEV,info,FONT=1

end
;===========================================================================
pro putheader,multix,multiy,ix,iy,ninfo,headline,nx

on_error,2
   
if ninfo lt 1 then return
info=strtrim(headline,2)
if ninfo gt 1 then info=info+' (nx='+string(nx,format='(3(i4))')+')'
xyouts,5+(ix*!d.x_size)/multix,-12+((iy+1)*!d.y_size)/multiy,/DEV,info,FONT=1

end
;===========================================================================
function diff2,direction,a,x
;
; Take derivative of "a" with respect to "x" in the direction "direction" 
; using 2nd order centered differencing
;
;===========================================================================
on_error,2

siz=size(a)
if siz(0) ne 2 then begin
   print,'Function diff2 is intended for 2D arrays only'
   retall
endif

n1=siz(1)
n2=siz(2)

if direction eq 1 then begin
   ind1=indgen(n1)
   jnd1=ind1+1
   jnd1(n1-1)=n1
   hnd1=ind1-1
   hnd1(0)=0
   dadx=(a(jnd1,*)-a(hnd1,*))/(x(jnd1,*)-x(hnd1,*))
endif
if direction eq 2 then begin
   ind2=indgen(n2)
   jnd2=ind2+1
   jnd2(n2-1)=n2
   hnd2=ind2-1
   hnd2(0)=0
   dadx=(a(*,jnd2)-a(*,hnd2))/(x(*,jnd2)-x(*,hnd2))
endif

return,dadx

end

;===========================================================================
function diff4,direction,a,x
;
; Take derivative of "a" with respect to "x" in the direction "direction" 
; using 4th order centered differencing
;
;===========================================================================
on_error,2

siz=size(a)
if siz(0) ne 2 then begin
   print,'Function diff4 is intended for 2D arrays only'
   retall
endif

n1=siz(1)
n2=siz(2)

dadx=a

if direction eq 1 then begin
   if n1 lt 5 then begin
      print,'Cannot take 4th order X gradient of grid with less than 5 columns'
      retall
   endif
   dadx(2:n1-3,*)=(a(4:n1-1,*)-8*a(3:n1-2,*)+8*a(1:n1-4,*)-a(0:n1-5,*)) $
                 /(x(3:n1-2,*)-x(1:n1-4,*))/6;
   dadx(0,*)   =dadx(2,*)
   dadx(1,*)   =dadx(2,*)
   dadx(n1-2,*)=dadx(n1-3,*)
   dadx(n1-1,*)=dadx(n1-3,*)
endif
if direction eq 2 then begin
   if n2 lt 5 then begin
      print,'Cannot take 4th order Y gradient of grid with less than 5 rows'
      retall
   endif
   dadx(*,2:n2-3)=(a(*,4:n2-1)-8*a(*,3:n2-2)+8*a(*,1:n2-4)-a(*,0:n2-5)) $
                 /(x(*,3:n2-2)-x(*,1:n2-4))/6;
   dadx(*,0)   =dadx(*,2)
   dadx(*,1)   =dadx(*,2)
   dadx(*,n2-2)=dadx(*,n2-3)
   dadx(*,n2-1)=dadx(*,n2-3)
endif

return,dadx

end

;===========================================================================
function diff3,direction,a,x
;
; Take derivative of "a" with respect to "x" in the direction "direction" 
; using IDL's 1D deriv() function
;
;===========================================================================
on_error,2

siz=size(a)
if siz(0) ne 2 then begin
   print,'Function diff3 is intended for 2D arrays only'
   retall
endif

dadx=a

if direction eq 1 then for i2=0,siz(2)-1 do dadx(*,i2)=deriv(x(*,i2),a(*,i2))
if direction eq 2 then for i1=0,siz(1)-1 do dadx(i1,*)=deriv(x(i1,*),a(i1,*))

return,dadx

end

;===========================================================================
function filledge,a

; On the edges use copy of closest cells

siz=size(a)
n1=siz(1)
n2=siz(2)

result=a
result(0,*)   =result(1,*)
result(*,0)   =result(*,1)
result(n1-1,*)=result(n1-2,*)
result(*,n2-1)=result(*,n2-2)

return,result
end

;===========================================================================
pro gengrid,name,x,y,xc,yc,vol2,u,v
;
; From cell center coordinates x,y calculate cell corner coordinates xc,yc,
; cell volumes. Check for array sizes of the optional u,v arguments.
; The name of the calling function is shown for error messages.
;===========================================================================

siz=size(x)
if siz(0) ne 2 then begin
   print,'Function ',name,' is for 2D arrays only'
   retall
endif

n1=siz(1)
n2=siz(2)

error=''
siz=size(y)
if siz(0) ne 2 or siz(1) ne n1 or siz(2) ne n2 then error='2nd coord'
if keyword_set(u) then begin
   siz=size(u)
   if siz(0) ne 2 or siz(1) ne n1 or siz(2) ne n2 then error='1st func'
endif
if keyword_set(v) then begin
   siz=size(v)
   if siz(0) ne 2 or siz(1) ne n1 or siz(2) ne n2 then error='2nd func'
endif
if error ne '' then begin
  print,'In function ',name,' the first argument does not match the ',error,'.'
  retall
endif

; Coordinates for cell corners
xc=(x(0:n1-2,0:n2-2)+x(0:n1-2,1:n2-1)+x(1:n1-1,0:n2-2)+x(1:n1-1,1:n2-1))/4
yc=(y(0:n1-2,0:n2-2)+y(0:n1-2,1:n2-1)+y(1:n1-1,0:n2-2)+y(1:n1-1,1:n2-1))/4

; Calculate 2*volume=(diagonal_1 X diagonal_2)
vol2=dblarr(n1,n2)+1
vol2(1:n1-2,1:n2-2)= $
 ((xc(1:n1-2,1:n2-2)-xc(0:n1-3,0:n2-3))*(yc(0:n1-3,1:n2-2)-yc(1:n1-2,0:n2-3)) $
 -(yc(1:n1-2,1:n2-2)-yc(0:n1-3,0:n2-3))*(xc(0:n1-3,1:n2-2)-xc(1:n1-2,0:n2-3)))

end

;===========================================================================
function intedge,f,xc
;
; Integrate the neighbouring values of "f" for the four edges described by "xc"
; The size of "f", "xc", and the result are n1*n2, (n1-1)*(n2-1), and n1*n2
; respectively, but only the inner (n1-2)*(n2-2) points are calculated, the
; edge values are 0-s.
;===========================================================================

siz=size(f)
n1=siz(1)
n2=siz(2)

intf=dblarr(n1,n2)
intf(1:n1-2,1:n2-2)=-(xc(1:n1-2,1:n2-2)-xc(0:n1-3,1:n2-2))*f(1:n1-2,2:n2-1) $
                    -(xc(1:n1-2,0:n2-3)-xc(1:n1-2,1:n2-2))*f(2:n1-1,1:n2-2) $
                    -(xc(0:n1-3,0:n2-3)-xc(1:n1-2,0:n2-3))*f(1:n1-2,0:n2-3) $
                    -(xc(0:n1-3,1:n2-2)-xc(0:n1-3,0:n2-3))*f(0:n1-3,1:n2-2)

return,intf

end

;===========================================================================
function intedge_rz,f,rc,zc
;
; Integrate r_edge*f_neighbour*dz for the four cell edges.
; assuming axial symmetry in the ignored direction.
; Only the inner (n1-2)*(n2-2) points are calculated, the edge values are 0-s.
;
;===========================================================================

siz=size(f)
n1=siz(1)
n2=siz(2)

intf=dblarr(n1,n2)
intf(1:n1-2,1:n2-2)= $
    -f(1:n1-2,2:n2-1)*(rc(1:n1-2,1:n2-2)+rc(0:n1-3,1:n2-2)) $
                     *(zc(1:n1-2,1:n2-2)-zc(0:n1-3,1:n2-2)) $
    -f(2:n1-1,1:n2-2)*(rc(1:n1-2,0:n2-3)+rc(1:n1-2,1:n2-2)) $
                     *(zc(1:n1-2,0:n2-3)-zc(1:n1-2,1:n2-2)) $
    -f(1:n1-2,0:n2-3)*(rc(0:n1-3,0:n2-3)+rc(1:n1-2,0:n2-3)) $
                     *(zc(0:n1-3,0:n2-3)-zc(1:n1-2,0:n2-3)) $
    -f(0:n1-3,1:n2-2)*(rc(0:n1-3,1:n2-2)+rc(0:n1-3,0:n2-3)) $
                     *(zc(0:n1-3,1:n2-2)-zc(0:n1-3,0:n2-3))

return,intf

end

;===========================================================================
function grad,idir,f,x,y
;
; Take gradient of "f" in direction "idir" on the "x,y" structured 2D grid.
; Gradient is the contour integral of edge_normal_idir*f_edge_averaged 
; divided by cell_volume for each cells. The cell corners are at the
; averaged coordinates of the four neighboring cell centers.
; However there is no need for edge averaging since the contribution of
; the value in the cell center cancels.
; Gradient can be calculated for inner points only, edge values are 
; copies of inner neighbors.
;===========================================================================

if n_elements(ndir) eq 0 or n_elements(f) eq 0 $
   or n_elements(x) eq 0 or n_elements(y) eq 0 then begin
   print,'Missing arguments in function grad'
   retall
endif

gengrid,'grad',x,y,xc,yc,vol2,f

if idir eq 1 then return,filledge( intedge(f,yc)/vol2) $
else              return,filledge(-intedge(f,xc)/vol2)

end

;===========================================================================
function grad_rz,idir,f,r,z
;
; Take gradient of "f" in direction "idir" on the "r,z" structured 2D grid
; assuming axial symmetry in the ignored direction.
; Gradient is the contour integral of edge_normal_idir*f*R_edge_averaged 
; divided by R*cell_volume - f/r for each cells. The cell corners are at the
; averaged coordinates of the four neighboring cell centers.
; However there is no need for edge averaging since the contribution of
; the value in the cell center cancels for idir=2, or equals +f/2R for idir=1.
; Gradient can be calculated for inner points only, edge values are 
; copies of inner neighbors.
;===========================================================================

if n_elements(idir) eq 0 or n_elements(f) eq 0 $
   or n_elements(r) eq 0 or n_elements(z) eq 0 then begin
   print,'Missing arguments in function grad_rz'
   retall
endif

gengrid,'grad_rz',r,z,rc,zc,vol2,f

if idir eq 1 then return,filledge( (intedge_rz(f,rc,zc)/vol2 - f)/2/r ) $
else              return,filledge( -intedge(f,rc^2)/vol2/2/r)

end

;===========================================================================
function div,u,v,x,y
;
; Take divergence of "u,v" vector with respect to "x,y" on a structured 2D grid
; Divergence is the contour integral of edge_normal.(u,v)_edge_averaged 
; divided by cell_volume for each cells. The cell corners are at the
; averaged coordinates of the four neighboring cell centers.
; However there is no need for edge averaging since the contribution of
; the value in the cell center cancels.
; Divergence can be calculated for inner points only, edge values are 
; copies of inner neighbors.
;===========================================================================

if n_elements(u) eq 0 or n_elements(v) eq 0 $
   or n_elements(x) eq 0 or n_elements(y) eq 0 then begin
   print,'Missing arguments in function div'
   retall
endif

gengrid,'div',x,y,xc,yc,vol2,u,v

return,filledge((intedge(u,yc)-intedge(v,xc))/vol2)

end

;===========================================================================
function div_rz,u,v,r,z
;
; Take divergence of "u,v" vector with respect to "r,z" on a structured 2D grid
; assuming axial symmetry in the ignored direction.
; Divergence is the contour integral of edge_normal.(u,v)*R_edge_averaged 
; divided by R*cell_volume for each cells. The cell corners are at the
; averaged coordinates of the four neighboring cell centers.
; However there is no need for edge averaging since the contribution of
; the value in the cell center is simply u/(2R).
; Divergence can be calculated for inner points only, edge values are 
; copies of inner neighbors.
;===========================================================================

if n_elements(u) eq 0 or n_elements(v) eq 0 $
   or n_elements(r) eq 0 or n_elements(z) eq 0 then begin
   print,'Missing arguments in function div_rz'
   retall
endif

gengrid,'div_rz',r,z,rc,zc,vol2,u,v

return,filledge(((intedge_rz(u,rc,zc)-intedge(v,rc^2))/vol2 + u)/2/r)

end

;===========================================================================
function curl,u,v,x,y
;
; Take curl of "u,v" vector with respect to "x,y" on a structured 2D grid.
; Curl is the contour integral of edge_vector.(u,v)_edge_averaged 
; divided by cell_volume for each cells. See also comments for div function.
;
;===========================================================================

if n_elements(u) eq 0 or n_elements(v) eq 0 $
   or n_elements(x) eq 0 or n_elements(y) eq 0 then begin
   print,'Missing arguments in function curl'
   retall
endif

gengrid,'curl',x,y,xc,yc,vol2,u,v

return,filledge((intedge(u,xc)+intedge(v,yc))/vol2)

end

;===========================================================================
function curl_rz,u,v,r,z
;
; Take curl of "u,v" vector with respect to "r,z" on a structured 2D grid
; with axial symmetry in the ignored direction.
; Curl is the contour integral of edge_vector.(u,v)*R_edge_averaged 
; divided by R*cell_volume for each cells - v/R. 
; See also comments for the div_rz function on edge average and edge cells.
;
;===========================================================================

if n_elements(u) eq 0 or n_elements(v) eq 0 $
   or n_elements(r) eq 0 or n_elements(z) eq 0 then begin
   print,'Missing arguments in function curl_rz'
   retall
endif

gengrid,'curl',r,z,rc,zc,vol2,u,v

return,filledge(-((intedge_rz(v,rc,zc)+intedge(u,rc^2))/vol2 - v)/2/r)

end

;==========================================
function coarse,a,boxsize
;
; Produce a coarser array from "a" by averaging out cells in a box.
; The box size can be defined by a scalar (n long interval, n*n squarle, 
; or ,n*n*n cube) or as an
; array of the same dimension as "a" (n1*n2 rectangle or n1*n2*n3 brick)

;on_error,2

if(n_elements(a) eq 0 or n_elements(boxsize) eq 0)then begin
   print,'Calling sequence is: array_co=coarse(array, boxsize)'
   retall
endif

siz=size(a)
ndim=siz(0)

if(ndim eq 0 or ndim gt 3)then begin
   print,'coarse requires a 1,2 or 3D array for the 1st argument'
   retall
endif
nx=siz(1:ndim)

siz=size(box)
if(siz(0) eq 0)then begin
   n=intarr(ndim)+boxsize
endif else if siz(0) eq ndim then begin
   n=boxsize
endif else begin
   print,'boxsize should either be a scalar, or an array '
   print,'of the same dimension as the number of dimensions of the array'
   retall
endelse

case ndim of
   1: begin
      result=dblarr(nx(0)/n(0))
      for ix=0,(nx(0)-1)/n(0) do $
	for i=0,n(0)-1 do $
           result(ix)=result(ix)+a(ix*n(0)+i)
      result=result/n(0)
   end
   2: begin
      result=dblarr(nx(0)/n(0),nx(1)/n(1))
      for ix=0,(nx(0)-1)/n(0) do $
      for iy=0,(nx(1)-1)/n(1) do $
	for i=0,n(0)-1 do $
	for j=0,n(1)-1 do $
           result(ix,iy)=result(ix,iy)+a(ix*n(0)+i,iy*n(1)+j)
      result=result/n(0)/n(1)
   end
   3: begin
      result=dblarr(nx(0)/n(0),nx(1)/n(1),nx(2)/n2)
      for ix=0,(nx(0)-1)/n(0) do $
      for iy=0,(nx(1)-1)/n(1) do $
      for iz=0,(nx(2)-1)/n(2) do $
	for i=0,n(0)-1 do $
	for j=0,n(1)-1 do $
	for k=0,n(2)-1 do $
           result(ix,iy,iz)=result(ix,iy,iz)+a(ix*n(0)+i,iy*n(1)+j,iz*n(2)+k)
      result=result/n(0)/n(1)/n(2)
   end
endcase
return,result
end

;===========================================================================
function quadruplet,nx,x0,x1,dx,ny,y0,y1,dy,nz,z0,z1,dz,nw,w0,w1,dw
;
; Produce an index array corresponding to the Fortran 90 triplet notation
;
; Usage: cut=quadruplet(100,0,30,2,100,30,40,1)
;
;        velvector=25*25  &  velpos=dblarr(velvector,2)
;        velpos(*,*)=x(quadruplet(100,0,99,4,100,30,69,2,2,0,1,1))
;===========================================================================

if keyword_set(dx) then begin
   checkdim,1,nx,x0,x1,dx
   all=lindgen(x1+1)
   sub=all(x0:x1)
   ind=sub(where(sub mod dx eq x0 mod dx))
end
if keyword_set(dy) then begin
   checkdim,2,ny,y0,y1,dy
   ixs=ind
   all=lindgen(y1+1)
   sub=all(y0:y1)
   iys=sub(where(sub mod dy eq y0 mod dy))
   ind=(ixs # (0*iys+1)) + ((0*ixs+nx) # iys)
end
if keyword_set(dz) then begin
   checkdim,3,nz,z0,z1,dz
   ixys=ind
   nxy=long(nx)*long(ny)
   all=lindgen(z1+1)
   sub=all(z0:z1)
   izs=sub(where(sub mod dz eq z0 mod dz))
   ind=lonarr(n_elements(ixs),n_elements(iys),n_elements(izs))
   for iz=0,n_elements(izs)-1 do ind(*,*,iz)=ixys + izs(iz)*nxy
end
if keyword_set(dw) then begin
   checkdim,4,nw,w0,w1,dw
   ixyzs=ind
   nxyz=long(nx)*long(ny)*long(nz)
   all=lindgen(w1+1)
   sub=all(w0:w1)
   iws=sub(where(sub mod dw eq w0 mod dw))
   ind=lonarr(n_elements(ixs),n_elements(iys),n_elements(izs),n_elements(iws))
   for iw=0,n_elements(iws)-1 do ind(*,*,*,iw)=ixyzs + iws(iw)*nxyz
end

return,ind
end

;===========================================================================
function triplet,x0,x1,dx,y0,y1,dy,z0,z1,dz,w0,w1,dw
;
; Produce an index array corresponding to the Fortran 90 triplet notation
;
; Usage: cut=triplet(0,99,2,0,99,2)
;
;        velvector=25*25  &  velpos=dblarr(velvector,2)
;        velpos(*,*)=x(triplet(0,99,4,0,99,4,0,1,1))
;
; Note: the resulting indices are valid for an array of size 
;
;      (x1+1)*(y1+1)*(z1+1)*(w1+1)
;
;===========================================================================

if keyword_set(dw) then $
   return,quadruplet(x1+1,x0,x1,dx,y1+1,y0,y1,dy,z1+1,z0,z1,dz,w1+1,w0,w1,dw)

if keyword_set(dz) then $
   return,quadruplet(x1+1,x0,x1,dx,y1+1,y0,y1,dy,z1+1,z0,z1,dz)

if keyword_set(dy) then $
   return,quadruplet(x1+1,x0,x1,dx,y1+1,y0,y1,dy)

if keyword_set(dx) then $
   return,quadruplet(x1+1,x0,x1,dx)

print,'Error in TRIPLET: All strides are 0!'
retall

end

;===========================================================================
pro checkdim,idim,nx,x0,x1,dx

; Check quadruplet for conditions nx>x1>=x0>=0 and dx>0
;===========================================================================

   if nx le 0 then begin
      print,'Size must be positive for dimension',idim
      retall
   endif
   if x1 ge nx then begin
      print,'Maximum index must be less than size for dimension',idim
      retall
   endif
   if x0 lt 0 then begin
      print,'Minimum index must be greater than 0 for dimension',idim
      retall
   endif
   if x0 gt x1 then begin
      print,'Minimum index must be less than maximum index for dimension',idim
      retall
   endif
   if dx le 0 then begin
      print,'Stride must be a positive integer for dimension',idim
      retall
   endif

return
end

;==========================================
pro plotgrid,x
;==========================================

on_error,2

s=size(x)
if s(0) ne 3 then begin
   print,'Error: plotgrid works for 2D grids only, like x(nx,ny,2)'
   retall
endif
if s(3) ne 2 then begin
   print,'Error: the third dimension should go from 0 to 1, like x(nx,ny,2)'
   retall
endif

plot,x(*,*,0),x(*,*,1),XSTYLE=1,YSTYLE=1,/NOERASE,/NODATA
for ix=0,s(1)-1 do begin
   oplot,x(ix,*,0),x(ix,*,1)
endfor
for iy=0,s(2)-1 do begin
   oplot,x(*,iy,0),x(*,iy,1)
endfor

return
end

;==========================================
pro compare,w1,w2

; Compare all variables in w1 and w2 by calculating 
; relative difference in the 1st norm.
;==========================================

on_error,2

sizew=size(w1)
ndim=sizew(0)-1
if ndim eq 0 then begin
   ndim=1
   nw=1
endif else $
   nw=sizew(ndim+1)

print,'iw max(|w1-w2|)/max(|w1|+|w2|) sum(|w1-w2|)/sum(|w1|+|w2|)'

for iw=0,nw-1 do begin
   case ndim of
   1: begin
      wsum=max(abs(w1(*,iw))+abs(w2(*,iw)))
      wdif=max(abs(w1(*,iw)-w2(*,iw)))
      wsum1=total(abs(w1(*,iw))+abs(w2(*,iw)))
      wdif1=total(abs(w1(*,iw)-w2(*,iw)))
      end
   2: begin
      wsum=max(abs(w1(*,*,iw))+abs(w2(*,*,iw)))
      wdif=max(abs(w1(*,*,iw)-w2(*,*,iw)))
      wsum1=total(abs(w1(*,*,iw))+abs(w2(*,*,iw)))
      wdif1=total(abs(w1(*,*,iw)-w2(*,*,iw)))
      end
   3: begin
      wsum=max(abs(w1(*,*,*,iw))+abs(w2(*,*,*,iw)))
      wdif=max(abs(w1(*,*,*,iw)-w2(*,*,*,iw)))
      wsum1=total(abs(w1(*,*,*,iw))+abs(w2(*,*,*,iw)))
      wdif1=total(abs(w1(*,*,*,iw)-w2(*,*,*,iw)))
      end
   endcase

   if wsum eq 0. then print,iw,' wsum=0' $
   else               print,iw,wdif/wsum,wdif1/wsum1
endfor
end

;==========================================
pro quit
   exit
end
;==========================================

;===========================================================================
function interp_orth,f,x,y,z,nxreg,tri_xy

; Interpolate in 3D for a grid, which is irregular in the X-Y planes but
; all X-Y planes are identical with each other and orthogonal to the Z 
; direction. The grid spacing in direction Z can be non-uniform.
;
; f      is the 3D function to interpolate
; x, y   are 2D arrays for the X and Y coordinates
; z      is a 1D array for the Z coordinate
; nxreg  is a 3 element array containing the size for the regular grid
; tri_xy is the 2D triangulation for any of the X-Y planes

nx=nxreg(0) & ny=nxreg(1) & nz=nxreg(2)

freg=fltarr(nx,ny,nz)

zmax=max(z) & zmin=min(z)

dx=(max(x)-min(x))/(nx-1.000001)
dy=(max(y)-min(y))/(ny-1.000001)
dz=(zmax  -zmin  )/(nz-1.000001)

iz=0                                  ; Initialize index for lower plane
jz=-1                                 ; Initialize index for upper plane
for izreg=0,nz-1 do begin
   zreg = zmin + izreg*dz             ; Z coordinate in regular grid

   while zreg gt z(iz+1) do iz=iz+1   ; Find elements of z enclosing zreg

   if iz gt jz-1 then begin           ; Check if iz has changed at all
      zi=z(iz)                        ; New Lower plane: zi<=zreg

      if jz eq iz then begin
         fi=fj                        ; Old upper plane is the new lower plane
      endif else begin
                                      ; Interpolate f in new lower plane
         fi=trigrid(x,y,f(*,iz,*),tri_xy, [dx,dy], missing=0) 
      endelse

      jz=iz+1 & zj=z(jz)              ; New Upper plane: zj>=zreg

                                      ; Interpolate f in new upper plane
      fj=trigrid(x,y,f(*,jz,*),tri_xy, [dx,dy], missing=0) 
   endif

   freg(*,*,izreg)=(fi*(zj-zreg)+fj*(zreg-zi))/(zj-zi)  ; Interpolate in Z
endfor

return,freg
end
;===========================================================================

pro PROJVOLUME, vol0, opaque0, AX=ax, AY=ay, AZ=az, $
    WINDOW=window, XSIZE=xsize, YSIZE=ysize, $
    XRES=xres, YRES=yres, ZRES=zres

; Default rotation angles
if not keyword_set(ax) then ax=0
if not keyword_set(ay) then ay=0
if not keyword_set(az) then az=0

; Default window parameters
if not keyword_set(window) then window=0
if not keyword_set(xsize) then xsize=512
if not keyword_set(ysize) then ysize=512

; Default resolution
if not keyword_set(xres) then xres=64
if not keyword_set(yres) then yres=64
if not keyword_set(zres) then zres=64

; Convert intensity to byte
vol=bytscl(vol0)
; Convert opacity to 0 to 25 range
opaque=bytscl(opaque0,top=25b)

; Set up viewing angle
s=size(vol)
xmin=0 & ymin = 0 & zmin = 0 & xmax=s(1) & ymax=s(2) & zmax=s(3)

scale3,xrange=[0,xmax-1],yrange=[0,ymax-1],zrange=[0,zmax-1],ax=ax,az=az

;!X.S = [-xmin, 1.0] / (xmax - xmin)
;!Y.S = [-ymin, 1.0] / (ymax - ymin)
;!Z.S = [-zmin, 1.0] / (zmax - zmin)
;T3D, /RESET 
;T3D, TRANSLATE=[-0.5, -0.5, -0.5] 
;T3D, SCALE=[0.7, 0.7, 0.7] 
;T3D, ROTATE=[ax,ay,az] 
;T3D, TRANSLATE=[0.5, 0.5, 0.5]

WINDOW, window, XSIZE=xsize, YSIZE=ysize

img = PROJECT_VOL(vol, xres, yres, zres, OPAQUE=opaque, TRANS=(!P.T))
help,img

TVSCL, img

end
;===========================================================================

pro  SHOWVOLUME, vol, thresh, COLOR = color, LOW = low

s = SIZE(vol)

IF s[0] NE 3 THEN begin
   print, 'Showvolume works for 3D arrays only'
   retall
endif

SCALE3, XRANGE=[0, S[1]], YRANGE=[0, S[2]], ZRANGE=[0, S[3]]

IF N_ELEMENTS(low) EQ 0 THEN low = 0

erase

if(keyword_set(COLOR))then begin
    col=color
    SHADE_VOLUME, vol, thresh, v, p, SHADES = col, LOW = low
    TVSCL, POLYSHADE(v,p,SHADES=col,/T3D)
endif else begin
    SHADE_VOLUME, vol, thresh, v, p, LOW = low
    TVSCL, POLYSHADE(v,p,/T3D)
endelse

end

;===========================================================================
