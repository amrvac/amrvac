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


;
; Changes by R. Keppens for AMR filetypes in
;         openfile/gettype/gethead/getpict/readtransform/readplotpar/plotfunc
; New procedures specific for AMR: 
;         getpict_asc_amr/getgrid_amr/getdomain_amr/getpict_bin_amr
;         cut2dto1d_amr
;
;==========================================
pro openfile,unit,filename,filetype
;==========================================
   on_error,2

   close,unit
   case filetype of
       ; ***AMR: so far only ascii or binary input
       'ascii_amr':openr,unit,filename
       'bin_amr'  :openr,unit,filename,/f77_unf
       'ascii'    :openr,unit,filename
       'binary'   :openr,unit,filename,/f77_unf
       'real4'    :openr,unit,filename,/f77_unf
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
      if len ne 131 then ftype='ascii' else begin
         ; The length of the 2nd line decides between real4 and binary
         ; since it contains the time, which is real*8 or real*4
         head=bytarr(131+4)
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
      ; ***AMR: gethead will reset ftype
      gethead,1,ftype,phys, $
         headline,it,time,gencoord,ndim,neqpar,nw,nx
      nxs=1
      for idim=1,ndim do nxs=nxs*nx(idim-1)
      case ftype of
       ; ***AMR: pictsize can not be determined a priori, and can
       ;         vary from snapshot to snapshot
       'ascii_amr': pictsize=-fsize
       'bin_amr': pictsize=-fsize
       'ascii': pictsize=1+79 + 1+7+13+9 + 1+ndim*4 + 1+neqpar*13 + 1+79 + $
                  (18*(ndim+nw)+1)*nxs
       'binary':pictsize=8+79 + 8+4*4+8  + 8+ndim*4 + 8+neqpar*8  + 8+79 + $
		  8*(1+nw)+8*(ndim+nw)*nxs
       'real4':pictsize=8+79 + 5*4+8  + 8+ndim*4 + 8+neqpar*4  + 8+79 + $
		  8*(1+nw)+4*(ndim+nw)*nxs
      endcase
      close,1

      ; Calculate number of snapshots in the file
      ; ***AMR: phony value of -1 for npictinfiles
      ;print,'ftype now=',ftype
      npictinfiles(ifile)= fsize/pictsize
      filetypes(ifile)=ftype
   endfor
end

;==========================================
pro gethead,unit,filetype,physics,headline,it,time,gencoord, $
            ndim,neqpar,nw,nx,eqpar,variables, $
            nxamr,corneramr,infoamr,ngrids
;==========================================
   on_error,2

;Type definitions
   headline='                                                                               '
   it=long(1)
   ndim=long(1)
   neqpar=long(1)
   nw=long(1)
   varname='                                                                                '
; ***AMR: when called from gettype: filetype is still
;         ascii/binary/real4: amr case identified further on in
;         the header by a negative grid size nx, which for amr files
;         is minus the number of grids (repeated ndim times)
;         when called from getpict or animate: filetype has been
;         reset to include _amr: undo for case statement
   filetypeamr=filetype
   if filetype eq 'ascii_amr' then filetype='ascii'
   if filetype eq 'bin_amr'   then filetype='binary'
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
                  ; ***AMR: negative nx
                  if nx(0) le 0  then begin
                    filetypeamr='ascii_amr'
                    ngrids   =abs(nx(0))
                    ;print,'ngrids=',ngrids
                    nxamr    =lonarr(ndim,ngrids)
                    infoamr  =lonarr(ngrids)
                    corneramr=dblarr(2*ndim,ngrids)
                  endif else begin
                    nxamr    =lonarr(ndim,1)
                    corneramr=dblarr(2*ndim,1)
                  endelse
                  ; ***AMR: end ****
		  eqpar=dblarr(neqpar)
		  readf,unit,eqpar
		  readf,unit,varname
                  ; ***AMR: read sizes of all grids, grid numbers, level
                  ;         and corners
                  if filetypeamr eq 'ascii_amr' then begin
                    nxgrid=lonarr(ndim)
                    cornergrid=dblarr(2*ndim)
                    for igrid=0,ngrids-1 do begin
                     ;print,'igrid=',igrid
                     readf,unit,nxgrid
                     nxamr(*,igrid)=nxgrid(*)
                     ;print,'sizes=',nxamr(*,igrid)
                     readf,unit,level
                     infoamr(igrid)=level
                     ;print,'level=',infoamr(igrid)
                     readf,unit,cornergrid
                     corneramr(*,igrid)=cornergrid(*)
                     ;print,'corners=',corneramr(*,igrid)
                    endfor
                  endif
                  ; ***AMR: end ****
	       end
      'binary':begin
                  time=double(1)
		  readu,unit,headline
		  readu,unit,it,time,ndim,neqpar,nw
		  gencoord=(ndim lt 0)
		  ndim=abs(ndim)
		  nx=lonarr(ndim)
		  readu,unit,nx
                  ; ***AMR: negative nx
                  if nx(0) le 0  then begin
                    filetypeamr='bin_amr'
                    ngrids   =abs(nx(0))
                    ;print,'ngrids=',ngrids
                    nxamr    =lonarr(ndim,ngrids)
                    infoamr  =lonarr(ngrids)
                    corneramr=dblarr(2*ndim,ngrids)
                  endif else begin
                    nxamr    =lonarr(ndim,1)
                    corneramr=dblarr(2*ndim,1)
                  endelse
                  ; ***AMR: end ****
		  eqpar=dblarr(neqpar)
		  readu,unit,eqpar
		  readu,unit,varname
                  ; ***AMR: read sizes of all grids, grid numbers, level
                  ;         and corners
                  if filetypeamr eq 'bin_amr' then begin
                    nxgrid=lonarr(ndim)
                    cornergrid=dblarr(2*ndim)
                    for igrid=0,ngrids-1 do begin
                     ;print,'igrid=',igrid
                     readu,unit,nxgrid
                     nxamr(*,igrid)=nxgrid(*)
                     ;print,'sizes=',nxamr(*,igrid)
                     level=lonarr(1)
                     readu,unit,level
                     ;help,level
                     ;print,'level=',level
                     infoamr(igrid)=level
                     ;print,'level=',infoamr(igrid)
                     readu,unit,cornergrid
                     corneramr(*,igrid)=cornergrid(*)
                     ;print,'corners=',corneramr(*,igrid)
                    endfor
                  endif
                  ; ***AMR: end ****
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
                  nxamr    =lonarr(ndim,1)
                  corneramr=dblarr(2*ndim,1)
	       end
      else: print,'Gethead: unknown filetype',filetype
   endcase
; ***AMR: reset filetype
   filetype=filetypeamr
   variables=str_sep(strtrim(strcompress(varname),2),' ')
   tmp=str_sep(strtrim(headline,2),'_')
   if n_elements(tmp) eq 2 then begin
      headline=tmp(0)
      physics=tmp(1)
   endif
end

;==========================================
pro getpict,unit,filetype,npict,$
    x,w,headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,error, $
    nxamr,corneramr,infoamr,xamr,wamr,ngrids,infoprint
;==========================================
   on_error,2

   ipict=0
   while ipict lt npict and not eof(unit) do begin
      ipict=ipict+1
      case filetype of
         ; ***AMR: instead of filling x/w, we fill xamr/wamr which
         ;         contain all grids
         'ascii_amr': getpict_asc_amr,unit,npict,$
             x,w,headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,$
             nxamr,corneramr,infoamr,xamr,wamr,ngrids,infoprint
         'bin_amr'  : getpict_bin_amr,unit,npict,$
             x,w,headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,$
             nxamr,corneramr,infoamr,xamr,wamr,ngrids,infoprint
         'ascii': getpict_asc,unit,npict,$
             x,w,headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,$
             nxamr,cornernamr
         'binary': getpict_bin,unit,npict,$
             x,w,headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,$
             nxamr,corneramr
         'real4': getpict_real,unit,npict,$
             x,w,headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,$
             nxamr,corneramr
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
pro getpict_asc_amr,unit,npict,$
   x,w,headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,$
   nxamr,corneramr,infoamr,xamr,wamr,ngrids,infoprint
;==========================================
  on_error,2                      
  ;----------------------------------------
  ;Read header information
  ;----------------------------------------
  ftypedummy='ascii'
  gethead,unit,ftypedummy,physics,$
     headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,$
     nxamr,corneramr,infoamr,ngrids
  if infoprint eq 1 then begin
     print
     print,' ngrids=',ngrids
  endif
  ndim=abs(ndim)
  ;----------------------------------------
  ; Read coordinates and values row by row
  ;----------------------------------------
  wrow=dblarr(nw)
  xrow=dblarr(ndim)

  ;print,'getpict_asc_amr: ngrids=',ngrids
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
    nxmax0=max(nxamr(0,*))
    ;print,'nxmax0=',nxmax0
    xamr=dblarr(nxmax0,ndim,ngrids)
    wamr=dblarr(nxmax0,nw,ngrids)
    for igrid=0,ngrids-1 do begin
      for x0=0,nxamr(0,igrid)-1 do begin
        readf,unit,xrow,wrow
        xamr(x0,0:ndim-1,igrid)=xrow(0:ndim-1)
        wamr(x0,0:nw-1,igrid)  =wrow(0:nw-1)
      endfor
    endfor
  end
  ;-------------- 2D ----------------------
  2: begin
    nxmax0=max(nxamr(0,*))
    nxmax1=max(nxamr(1,*))
    ;print,'nxmax0 and nxmax1=',nxmax0,nxmax1
    xamr=dblarr(nxmax0,nxmax1,ndim,ngrids)
    wamr=dblarr(nxmax0,nxmax1,nw,ngrids)
    for igrid=0,ngrids-1 do begin
      for x1=0,nxamr(1,igrid)-1 do begin
        for x0=0,nxamr(0,igrid)-1 do begin
          readf,unit,xrow,wrow
          xamr(x0,x1,0:ndim-1,igrid)=xrow(0:ndim-1)
          wamr(x0,x1,0:nw-1,igrid)  =wrow(0:nw-1)
        endfor
      endfor
    endfor
  end
  ;-------------- 3D ----------------------
  3: begin
    nxmax0=max(nxamr(0,*))
    nxmax1=max(nxamr(1,*))
    nxmax2=max(nxamr(2,*))
    ;print,'nxmax0 nxmax1 nxmax2=',nxmax0,nxmax1,nxmax2
    xamr=dblarr(nxmax0,nxmax1,nxmax2,ndim,ngrids)
    wamr=dblarr(nxmax0,nxmax1,nxmax2,nw,ngrids)
    for igrid=0,ngrids-1 do begin
      for x2=0,nxamr(2,igrid)-1 do begin
        for x1=0,nxamr(1,igrid)-1 do begin
          for x0=0,nxamr(0,igrid)-1 do begin
            readf,unit,xrow,wrow
            xamr(x0,x1,x2,0:ndim-1,igrid)=xrow(0:ndim-1)
            wamr(x0,x1,x2,0:nw-1,igrid)=wrow(0:nw-1)
	  endfor
	endfor
      endfor
    endfor
  end
  endcase
end

;==========================================
pro getpict_bin_amr,unit,npict,$
   x,w,headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,$
   nxamr,corneramr,infoamr,xamr,wamr,ngrids,infoprint
;==========================================
  on_error,2                      
  ;----------------------------------------
  ;Read header information
  ;----------------------------------------
  ftypedummy='binary'
  gethead,unit,ftypedummy,physics,$
     headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,$
     nxamr,corneramr,infoamr,ngrids
  if infoprint eq 1 then begin
     print
     print,' ngrids=',ngrids
  endif
  ndim=abs(ndim)
  ;----------------------------------------
  ; Read coordinates and values 
  ;----------------------------------------
  wrow=dblarr(nw)
  xrow=dblarr(ndim)

  ;print,'getpict_bin_amr: ngrids=',ngrids
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
    nxmax0=max(nxamr(0,*))
    ;print,'nxmax0=',nxmax0
    xamr=dblarr(nxmax0,ndim,ngrids)
    wamr=dblarr(nxmax0,nw,ngrids)
    for igrid=0,ngrids-1 do begin
      n1g=nxamr(0,igrid)
      x1g=dblarr(n1g,ndim)
      w1g=dblarr(n1g,nw)
      wi=dblarr(n1g)
      readu,unit,x1g
      for iw=0,nw-1 do begin
         readu,unit,wi
         w1g(*,iw)=wi
      endfor
      xamr(0:n1g-1,0:ndim-1,igrid)=x1g(0:n1g-1,0:ndim-1)
      wamr(0:n1g-1,0:nw-1,igrid)  =w1g(0:n1g-1,0:nw-1)
    endfor
  end
  ;-------------- 2D ----------------------
  2: begin
    nxmax0=max(nxamr(0,*))
    nxmax1=max(nxamr(1,*))
    ;print,'nxmax0 and nxmax1=',nxmax0,nxmax1
    xamr=dblarr(nxmax0,nxmax1,ndim,ngrids)
    wamr=dblarr(nxmax0,nxmax1,nw,ngrids)
    for igrid=0,ngrids-1 do begin
      ;print,'doing igrid:',igrid
      n1g=nxamr(0,igrid)
      n2g=nxamr(1,igrid)
      ;print,'sizes:',n1g,n2g
      ;print,'ndim:',ndim
      x1g=dblarr(n1g,n2g,ndim)
      ;help,x1g
      w1g=dblarr(n1g,n2g,nw)
      wi=dblarr(n1g,n2g)
      ;print,'HERE1: start reading x'
      readu,unit,x1g
      ;print,'read x'
      for iw=0,nw-1 do begin
         readu,unit,wi
         ;print,'read w field number',iw
         w1g(*,*,iw)=wi
      endfor
      xamr(0:n1g-1,0:n2g-1,0:ndim-1,igrid)=x1g(0:n1g-1,0:n2g-1,0:ndim-1)
      wamr(0:n1g-1,0:n2g-1,0:nw-1,igrid)  =w1g(0:n1g-1,0:n2g-1,0:nw-1)
    endfor
  end
  ;-------------- 3D ----------------------
  3: begin
    nxmax0=max(nxamr(0,*))
    nxmax1=max(nxamr(1,*))
    nxmax2=max(nxamr(2,*))
    ;print,'nxmax0 nxmax1 nxmax2=',nxmax0,nxmax1,nxmax2
    xamr=dblarr(nxmax0,nxmax1,nxmax2,ndim,ngrids)
    wamr=dblarr(nxmax0,nxmax1,nxmax2,nw,ngrids)
    for igrid=0,ngrids-1 do begin
       n1g=nxamr(0,igrid)
       n2g=nxamr(1,igrid)
       n3g=nxamr(2,igrid)
       x1g=dblarr(n1g,n2g,n3g,ndim)
       w1g=dblarr(n1g,n2g,n3g,nw)
       wi=dblarr(n1g,n2g,n3g)
       readu,unit,x1g
       for iw=0,nw-1 do begin
          readu,unit,wi
          w1g(*,*,*,iw)=wi
       endfor
       xamr(0:n1g-1,0:n2g-1,0:n3g-1,0:ndim-1,igrid)= $
                   x1g(0:n1g-1,0:n2g-1,0:n3g-1,0:ndim-1)
       wamr(0:n1g-1,0:n2g-1,0:n3g-1,0:nw-1,igrid)=   $
                   w1g(0:n1g-1,0:n2g-1,0:n3g-1,0:nw-1)
    endfor
  end
  endcase
end

;==========================================
pro cut2dto1d_amr,ndimin,nxamrin,corneramrin,infoamrin,xamrin,wamrin,ngridsin, $
     cutx,vcutx,cuty,vcuty,nw,ndim,nxamr,corneramr,infoamr,xamr,wamr,ngrids
;==========================================
  on_error,2                      
  ;----------------------------------------
  if ndimin ne 2 then begin
     print,'Error: cut2dto1d for 2D data only!'
     retall
  endif
  if (cutx+cuty) ne 1 then begin
     print,'Error: cut2dto1d must have 1 cutting dimension set!'
     print,'  cutx=',cutx,' cuty=',cuty
  endif
  domain=dblarr(2*ndimin)
  domain(0)=min(corneramrin(0,*)) 
  domain(1)=min(corneramrin(1,*)) 
  domain(2)=max(corneramrin(2,*)) 
  domain(3)=max(corneramrin(3,*)) 
  if cutx eq 1 then begin
     ; cut along y-direction
     print,'Will do cut along y at constant x-value=',vcutx
     ; check whether cut value is interior to domain x-range
     if vcutx lt domain(0) or vcutx gt domain(2) then begin
        print,'Cut value exterior to domain range:',domain(0),domain(2)
        retall
     endif
     intersect=lonarr(ngridsin)
     ngrids=0
     for igrid=0,ngridsin-1 do begin
        if vcutx ge corneramrin(0,igrid) and vcutx lt corneramrin(2,igrid) $
           then begin
           ngrids=ngrids+1
           intersect(igrid)=1
        endif else begin
           intersect(igrid)=0
        endelse
     endfor
     idimset=1
     idimcut=0
  endif else begin
     ; cut along x-direction
     print,'Will do cut along x at constant y-value=',vcuty
     ; check whether cut value is interior to domain y-range
     if vcuty lt domain(1) or vcuty gt domain(3) then begin
        print,'Cut value exterior to domain range:',domain(1),domain(3)
        retall
     endif
     intersect=lonarr(ngridsin)
     ngrids=0
     for igrid=0,ngridsin-1 do begin
        if vcuty ge corneramrin(1,igrid) and vcuty lt corneramrin(3,igrid) $
           then begin
           ngrids=ngrids+1
           intersect(igrid)=1
        endif else begin
           intersect(igrid)=0
        endelse
     endfor
     idimset=0
     idimcut=1
  endelse

  nxamr=lonarr(1,ngrids)
  infoamr=lonarr(ngrids)
  corneramr=dblarr(2,ngrids)

  print,'Found 1D intersection with ',ngrids,' grids from total of ',ngridsin
  nxmax=max(nxamrin(idimset,*))
  print,'nxmax=',nxmax
  xamr=dblarr(nxmax,1,ngrids)
  wamr=dblarr(nxmax,nw,ngrids)
  igrid1d=0
  for igrid=0,ngridsin-1 do begin
    if intersect(igrid) eq 1 then begin
       print,'intersect with grid:',igrid
       nxamr(0,igrid1d)=nxamrin(idimset,igrid)
       nn1d=nxamr(0,igrid1d)
       infoamr(igrid1d)=infoamrin(igrid)
       corneramr(0,igrid1d)=corneramrin(idimset,igrid) 
       corneramr(1,igrid1d)=corneramrin(idimset+2,igrid) 
       print,'grid ',igrid1d,' range:',corneramr(0,igrid1d),corneramr(1,igrid1d)
       if cutx eq 1 then begin
          dxdim=xamrin(1,0,0,igrid)-xamrin(0,0,0,igrid)
          print,'dxdim=',dxdim,' vcutx=',vcutx
          index=(vcutx-corneramrin(idimcut,igrid))/dxdim
          print,'intersect at index:',index
          xamr(0:nn1d-1,0,igrid1d)=xamrin(index,0:nn1d-1,1,igrid)
          wamr(0:nn1d-1,0:nw-1,igrid1d)=wamrin(index,0:nn1d-1,0:nw-1,igrid)
       endif
       if cuty eq 1 then begin
          dxdim=xamrin(0,1,1,igrid)-xamrin(0,0,1,igrid)
          print,'dxdim=',dxdim,' vcuty=',vcuty
          index=(vcuty-corneramrin(idimcut,igrid))/dxdim
          print,'intersect at index:',index
          xamr(0:nn1d-1,0,igrid1d)=xamrin(0:nn1d-1,index,0,igrid)
          wamr(0:nn1d-1,0:nw-1,igrid1d)=wamrin(0:nn1d-1,index,0:nw-1,igrid)
       endif
       igrid1d=igrid1d+1
    endif
  endfor
  print,'Resetting ndim=1'
  ndim=1
end

;==========================================
pro cut3dto2d_amr,ndimin,nxamrin,corneramrin,infoamrin,xamrin,wamrin,ngridsin, $
     cutx,vcutx,cuty,vcuty,cutz,vcutz, $
     nw,ndim,nxamr,corneramr,infoamr,xamr,wamr,ngrids
;==========================================
  on_error,2                      
  ;----------------------------------------
  if ndimin ne 3 then begin
     print,'Error: cut3dto2d for 3D data only!'
     retall
  endif
  if (cutx+cuty+cutz) ne 1 then begin
     print,'Error: cut3dto2d must have 1 cutting dimension set!'
     print,'  cutx=',cutx,' cuty=',cuty,' cutz=',cutz
  endif
  domain=dblarr(2*ndimin)
  domain(0)=min(corneramrin(0,*)) 
  domain(1)=min(corneramrin(1,*)) 
  domain(2)=min(corneramrin(2,*)) 
  domain(3)=max(corneramrin(3,*)) 
  domain(4)=max(corneramrin(4,*)) 
  domain(5)=max(corneramrin(5,*)) 
  if cutx eq 1 then begin
     ; cut along yz-plane
     print,'Will do cut in y-z plane at constant x-value=',vcutx
     ; check whether cut value is interior to domain x-range
     if vcutx lt domain(0) or vcutx gt domain(3) then begin
        print,'Cut value exterior to domain range:',domain(0),domain(3)
        retall
     endif
     intersect=lonarr(ngridsin)
     ngrids=0
     for igrid=0,ngridsin-1 do begin
        if vcutx ge corneramrin(0,igrid) and vcutx lt corneramrin(3,igrid) $
           then begin
           ngrids=ngrids+1
           intersect(igrid)=1
        endif else begin
           intersect(igrid)=0
        endelse
     endfor
     idimcut=0
     idimset1=1
     idimset2=2
  endif 
  if cuty eq 1 then begin
     ; cut along xz-plane
     print,'Will do cut in x-z plane at constant y-value=',vcuty
     ; check whether cut value is interior to domain y-range
     if vcuty lt domain(1) or vcuty gt domain(4) then begin
        print,'Cut value exterior to domain range:',domain(1),domain(4)
        retall
     endif
     intersect=lonarr(ngridsin)
     ngrids=0
     for igrid=0,ngridsin-1 do begin
        if vcuty ge corneramrin(1,igrid) and vcuty lt corneramrin(4,igrid) $
           then begin
           ngrids=ngrids+1
           intersect(igrid)=1
        endif else begin
           intersect(igrid)=0
        endelse
     endfor
     idimcut=1
     idimset1=0
     idimset2=2
  endif
  if cutz eq 1 then begin
     ; cut along xy-plane
     print,'Will do cut in x-y plane at constant z-value=',vcutz
     ; check whether cut value is interior to domain z-range
     if vcutz lt domain(2) or vcutz gt domain(5) then begin
        print,'Cut value exterior to domain range:',domain(2),domain(5)
        retall
     endif
     intersect=lonarr(ngridsin)
     ngrids=0
     for igrid=0,ngridsin-1 do begin
        if vcutz ge corneramrin(2,igrid) and vcutz lt corneramrin(5,igrid) $
           then begin
           ngrids=ngrids+1
           intersect(igrid)=1
        endif else begin
           intersect(igrid)=0
        endelse
     endfor
     idimcut=2
     idimset1=0
     idimset2=1
  endif

  nxamr=lonarr(2,ngrids)
  infoamr=lonarr(ngrids)
  corneramr=dblarr(4,ngrids)

  print,'Found 2D intersection with ',ngrids,' grids from total of ',ngridsin
  nxmax1=max(nxamrin(idimset1,*))
  nxmax2=max(nxamrin(idimset2,*))
  print,'nxmax1=',nxmax1,' nxmax2=',nxmax2
  xamr=dblarr(nxmax1,nxmax2,2,ngrids)
  wamr=dblarr(nxmax1,nxmax2,nw,ngrids)
  igrid2d=0
  for igrid=0,ngridsin-1 do begin
    if intersect(igrid) eq 1 then begin
       nxamr(0,igrid2d)=nxamrin(idimset1,igrid)
       nxamr(1,igrid2d)=nxamrin(idimset2,igrid)
       nn2d1=nxamr(0,igrid2d)
       nn2d2=nxamr(1,igrid2d)
       infoamr(igrid2d)=infoamrin(igrid)
       corneramr(0,igrid2d)=corneramrin(idimset1,igrid) 
       corneramr(1,igrid2d)=corneramrin(idimset2,igrid) 
       corneramr(2,igrid2d)=corneramrin(idimset1+3,igrid) 
       corneramr(3,igrid2d)=corneramrin(idimset2+3,igrid) 
       print,'grid ',igrid2d,' range(1):', $
              corneramr(0,igrid2d),corneramr(2,igrid2d)
       print,'                 range(2):', $
              corneramr(1,igrid2d),corneramr(3,igrid2d)
       if cutx eq 1 then begin
          dxdim=xamrin(1,0,0,0,igrid)-xamrin(0,0,0,0,igrid)
          print,'dxdim=',dxdim,' vcutx=',vcutx
          index=(vcutx-corneramrin(idimcut,igrid))/dxdim
          print,'intersect at index:',index
          xamr(0:nn2d1-1,0:nn2d2-1,0:1,igrid2d)= $
              xamrin(index,0:nn2d1-1,0:nn2d2-1,1:2,igrid)
          wamr(0:nn2d1-1,0:nn2d2-1,0:nw-1,igrid2d)= $
              wamrin(index,0:nn2d1-1,0:nn2d2-1,0:nw-1,igrid)
       endif
       if cuty eq 1 then begin
          dxdim=xamrin(0,1,0,1,igrid)-xamrin(0,0,0,1,igrid)
          print,'dxdim=',dxdim,' vcuty=',vcuty
          index=(vcuty-corneramrin(idimcut,igrid))/dxdim
          print,'intersect at index:',index
          xamr(0:nn2d1-1,0:nn2d2-1,0,igrid2d)= $
              xamrin(0:nn2d1-1,index,0:nn2d2-1,0,igrid)
          xamr(0:nn2d1-1,0:nn2d2-1,1,igrid2d)= $
              xamrin(0:nn2d1-1,index,0:nn2d2-1,2,igrid)
          wamr(0:nn2d1-1,0:nn2d2-1,0:nw-1,igrid2d)= $
              wamrin(0:nn2d1-1,index,0:nn2d2-1,0:nw-1,igrid)
       endif
       if cutz eq 1 then begin
          dxdim=xamrin(0,0,1,2,igrid)-xamrin(0,0,0,2,igrid)
          print,'dxdim=',dxdim,' vcutz=',vcutz
          index=(vcutz-corneramrin(idimcut,igrid))/dxdim
          print,'intersect at index:',index
          xamr(0:nn2d1-1,0:nn2d2-1,0:1,igrid2d)= $
              xamrin(0:nn2d1-1,0:nn2d2-1,index,0:1,igrid)
          wamr(0:nn2d1-1,0:nn2d2-1,0:nw-1,igrid2d)= $
              wamrin(0:nn2d1-1,0:nn2d2-1,index,0:nw-1,igrid)
       endif
       igrid2d=igrid2d+1
    endif
  endfor
  print,'Resetting ndim=2'
  ndim=2
end

;==========================================
pro getdomain_amr,ndim,domainamr,corneramr,noautodomain,domainset
;==========================================
  on_error,2                      
  ;----------------------------------------
  if domainset eq 0 then begin
     domainset=1
  domainamr=dblarr(2*ndim)
  if noautodomain eq 1 then begin
     print,'give domain [2*ndim floats: min(s) -- max(s)]'    
     read,domainamr
     print,'-->Domain set to:',domainamr
  endif else begin
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
     domainamr(0)=min(corneramr(0,*)) 
     domainamr(1)=max(corneramr(1,*)) 
  end
  ;-------------- 2D ----------------------
  2: begin
     domainamr(0)=min(corneramr(0,*)) 
     domainamr(1)=min(corneramr(1,*)) 
     domainamr(2)=max(corneramr(2,*)) 
     domainamr(3)=max(corneramr(3,*)) 
  end
  ;-------------- 3D ----------------------
  3: begin
     domainamr(0)=min(corneramr(0,*)) 
     domainamr(1)=min(corneramr(1,*)) 
     domainamr(2)=min(corneramr(2,*)) 
     domainamr(3)=max(corneramr(3,*)) 
     domainamr(4)=max(corneramr(4,*)) 
     domainamr(5)=max(corneramr(5,*)) 
  end
  endcase
  endelse
  endif
  print,'domain size =',domainamr
end

;==========================================
pro getgrid_amr,ndim,nw,nx,x,w,nxamr,xamr,wamr,igrid
;==========================================
  on_error,2                      
  ;----------------------------------------
  ;print,'getting grid=',igrid
  nx=lonarr(ndim)
  nx(*)=nxamr(*,igrid)
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
     x=dblarr(nx(0),ndim)
     w=dblarr(nx(0),nw)
     x(*,*)=xamr(0:nx(0)-1,*,igrid) 
     w(*,*)=wamr(0:nx(0)-1,*,igrid)
  end
  ;-------------- 2D ----------------------
  2: begin
     x=dblarr(nx(0),nx(1),ndim)
     w=dblarr(nx(0),nx(1),nw)
     x(*,*,*)=xamr(0:nx(0)-1,0:nx(1)-1,*,igrid) 
     w(*,*,*)=wamr(0:nx(0)-1,0:nx(1)-1,*,igrid)
  end
  ;-------------- 3D ----------------------
  3: begin
     x=dblarr(nx(0),nx(1),nx(2),ndim)
     w=dblarr(nx(0),nx(1),nx(2),nw)
     x(*,*,*,*)=xamr(0:nx(0)-1,0:nx(1)-1,0:nx(2)-1,*,igrid) 
     w(*,*,*,*)=wamr(0:nx(0)-1,0:nx(1)-1,0:nx(2)-1,*,igrid)
  end
  endcase
end

;==========================================
pro getpict_asc,unit,npict,$
   x,w,headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,$
   nxamr,corneramr
;==========================================
  on_error,2                      
  ;----------------------------------------
  ;Read header information
  ;----------------------------------------
  gethead,unit,'ascii',physics,$
     headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,$
     nxamr,corneramr
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
    corneramr(0,0)=min(x(*,*))
    corneramr(1,0)=max(x(*,*))
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
    corneramr(0,0)=min(x(*,*,0))
    corneramr(1,0)=max(x(*,*,0))
    corneramr(2,0)=min(x(*,*,1))
    corneramr(3,0)=max(x(*,*,1))
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
    corneramr(0,0)=min(x(*,*,*,0))
    corneramr(1,0)=max(x(*,*,*,0))
    corneramr(2,0)=min(x(*,*,*,1))
    corneramr(3,0)=max(x(*,*,*,1))
    corneramr(4,0)=min(x(*,*,*,2))
    corneramr(5,0)=max(x(*,*,*,2))
  end
  endcase
end

;==========================================
pro getpict_bin,unit,npict,$
   x,w,headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,$
   nxamr,corneramr
;==========================================
  on_error,2                      
  ;----------------------------------------
  ;Read header information
  ;----------------------------------------
  gethead,unit,'binary',physics,$
     headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,$
     nxamr,corneramr
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
    corneramr(0,0)=min(x(*,*))
    corneramr(1,0)=max(x(*,*))
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
    corneramr(0,0)=min(x(*,*,0))
    corneramr(1,0)=max(x(*,*,0))
    corneramr(2,0)=min(x(*,*,1))
    corneramr(3,0)=max(x(*,*,1))
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
    corneramr(0,0)=min(x(*,*,*,0))
    corneramr(1,0)=max(x(*,*,*,0))
    corneramr(2,0)=min(x(*,*,*,1))
    corneramr(3,0)=max(x(*,*,*,1))
    corneramr(4,0)=min(x(*,*,*,2))
    corneramr(5,0)=max(x(*,*,*,2))
  end
  endcase
end

;==========================================
pro getpict_real,unit,npict,$
   x,w,headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,$
   nxamr,corneramr
;==========================================
  on_error,2                      
  ;----------------------------------------
  ;Read header information
  ;----------------------------------------
  gethead,unit,'real4',physics,$
     headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,$
     nxamr,corneramr
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
    corneramr(0,0)=min(x(*,*))
    corneramr(1,0)=max(x(*,*))
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
    corneramr(0,0)=min(x(*,*,0))
    corneramr(1,0)=max(x(*,*,0))
    corneramr(2,0)=min(x(*,*,1))
    corneramr(3,0)=max(x(*,*,1))
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
    corneramr(0,0)=min(x(*,*,*,0))
    corneramr(1,0)=max(x(*,*,*,0))
    corneramr(2,0)=min(x(*,*,*,1))
    corneramr(3,0)=max(x(*,*,*,1))
    corneramr(4,0)=min(x(*,*,*,2))
    corneramr(5,0)=max(x(*,*,*,2))
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
         read,PROMPT=prompt+'='+strtrim(string(var),2)+' ?? ',tmp
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
      endif else print,prompt,'= ',var
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
   plotmode,plotmodes,plottitle,plottitles,autorange,autoranges,doask, $
   transform,noautodomain,domainset
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
      if siz(0) eq 3 and siz(2) eq 1 then begin
          cut0=lonarr(siz(1),siz(3))
          cut0(*,*)=cut(*,0,*)
      endif
   endelse

   askstr,'func(s) (e.g. rho p m1;m2 r*m1 -T) ',func,doask
   if plotdim eq 1 then begin
      ; ***AMR: choose plot/plotamr
      if transform eq 'amr' then begin
      print,'1D plotmodes: plot_amr/plots_amr/plotn_amr/plotf_amr'
      askstr,'plotmode(s)                ',plotmode,doask
      print,'1D plotmode: ',plotmode
      if plotmode eq 'plotf_amr' then loadct,2
      if plotmode eq 'plots_amr' then loadct,2
      if plotmode eq 'plotn_amr' then loadct,2
      if domainset eq 0 then begin
         asknum,'noautodomain (0/1)         ',noautodomain,doask
      endif
      endif else begin
      plotmode='plot'
      print,'1D plotmode: ',plotmode
      endelse
   endif else begin
      if plotmode eq 'plot' then plotmode=''
      ; ***AMR: choose contfill_amr
      if transform eq 'amr' then begin
       print,'2D plotmodes: contfill_amr/contour_amr'
       askstr,'plotmode(s)                ',plotmode,doask
       print,'2D plotmode: ',plotmode
       if domainset eq 0 then begin
          asknum,'noautodomain (0/1)         ',noautodomain,doask
       endif
      endif else begin
       print,'2D plotmodes: contour/contlabel/contfill/shade/shadeirr/surface/tv/'
       print,'              vel/velovect/ovelovect/contfillB'
       askstr,'plotmode(s)                ',plotmode,doask
      endelse
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

 if transform eq 'amr' then begin
   print,'skipping use of grid in readtransform since amr file'
 endif else begin
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
 endelse
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


   for ifunc=0,nfunc-1 do begin
     if autoranges[ifunc] eq 'n' then begin
      if noautorange then begin
        print,'read f ranges:'
        f_min=0.
        f_max=0.
      endif else begin
        f_min=fmin(ifunc)
        f_max=fmax(ifunc)
      endelse
      asknum,'Min value for '+funcs(ifunc),f_min,doask
      asknum,'Max value for '+funcs(ifunc),f_max,doask
      fmin(ifunc)=f_min
      fmax(ifunc)=f_max
     endif
   endfor
   noautorange=0

end
;===========================================================================
pro getlimits,first,nfunc,funcs,funcs1,funcs2,autoranges,fmax,fmin,$
                x,w,xreg,wreg,usereg,physics,eqpar,wnames,cut,doask
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
function getvalue,xp,funcname,xamr,wamr,ngrids,ndim,nxamr,nw,physics,eqpar,wnames
; to get function value at a given position in amr data
; variables : an array of function names
;===========================================================================
if(n_elements(xp) eq 0) then message,'position is undefined'
if(n_elements(funcname) eq 0) then message,'function is undefined'
for igrid=0,ngrids-1 do begin
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
     x=dblarr(nxamr(0,igrid),ndim)
     x(*,0)=xamr(0:nxamr(0,igrid)-1,*,igrid)
     if( xp le max(x) and xp ge min(x)) then begin
       w=dblarr(nxamr(0,igrid),nw)
       w(*,*)=wamr(0:nxamr(0,igrid)-1,*,igrid)
       f=animfunc(x,w,funcname,physics,eqpar,wnames)
       ;jl=0
       ;jh=nxamr(0,igrid)+1
       ;while(jh-jl gt 1) do begin
       ;  jc=(jh+jl)/2
       ;  if (xp >= x(jc,0)) then begin
       ;   jl=jc
       ;  endif else begin
       ;   jh=jc
       ;  endelse
       ;endwhile
       ;value=f(jl)+(xp-x(jl,0))*(f(jh)-f(jl))/(x(jh,0)-x(jl,0))
       np=(xp-x(0,0))/(x(1,0)-x(0,0))
       value=interpolate(f,np)
     endif
  end
  ;-------------- 2D ----------------------
  2: begin
     x=dblarr(nx(0),nx(1),ndim)
     x(*,*,*)=xamr(0:nx(0)-1,0:nx(1)-1,*,igrid)
     if(xp(0) le max(x(*,*,0)) and xp(0) ge min(x(*,*,0)) $
       and xp(1) le max(x(*,*,1)) and xp(1) ge min(x(*,*,1))) then begin
       w=dblarr(nx(0),nx(1),nw)
       w(*,*,*)=wamr(0:nx(0)-1,0:nx(1)-1,*,igrid)
       f=animfunc(x,w,funcname,physics,eqpar,wnames)
       np1=(xp(0)-x(0,0,0))/(x(1,0,0)-x(0,0,0))
       np2=(xp(1)-x(0,0,1))/(x(0,1,1)-x(0,0,1))
       value=interpolate(f,np1,np2)
     endif
  end
  ;-------------- 3D ----------------------
  3: begin
     x=dblarr(nx(0),nx(1),nx(2),ndim)
     x(*,*,*,*)=xamr(0:nx(0)-1,0:nx(1)-1,0:nx(2)-1,*,igrid)
     if(xp(0) le max(x(*,*,*,0)) and xp(0) ge min(x(*,*,*,0)) $
       and xp(1) le max(x(*,*,*,1)) and xp(1) ge min(x(*,*,*,1))$
       and xp(2) le max(x(*,*,*,2)) and xp(2) ge min(x(*,*,*,2))) then begin
       w=dblarr(nx(0),nx(1),nx(2),nw)
       w(*,*,*,*)=wamr(0:nx(0)-1,0:nx(1)-1,0:nx(2)-1,*,igrid)
       f=animfunc(x,w,funcname,physics,eqpar,wnames)
       np1=(xp(0)-x(0,0,0,0))/(x(1,0,0,0)-x(0,0,0,0))
       np2=(xp(1)-x(0,0,0,1))/(x(0,1,0,1)-x(0,0,0,1))
       np3=(xp(2)-x(0,0,0,2))/(x(0,0,1,2)-x(0,0,0,2))
       value=interpolate(f,np1,np2)
     endif
  end
  endcase
endfor
return,value
end

;===========================================================================
pro getfunc,f,f1,f2,func1,func2,x,w,physics,eqpar,wnames,cut
;===========================================================================
on_error,2

; Added following statement for trailing blank removal
physics=strcompress(physics,/remove_all)

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
  nfunc,multix,multiy,plotix,plotiy,funcs,funcs1,funcs2,fmin,fmax,f,$
  domainamr,info1amr,corner1amr,rgb,vertical=vertical
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
   if n_elements(vertical) eq 0 then vertical=0 ;;if 1 then plot maps in a column
   if vertical then begin
     multix=1
     multiy=nfunc
   endif
   ppos=fltarr(4,multix*multiy)
   for iy=0,multiy-1 do begin
     for ix=0,multix-1 do begin
        ppos(0,iy*multix+ix)=float(ix)/float(multix)
        ppos(1,iy*multix+ix)=1.-float(iy+1)/float(multiy)
        ppos(2,iy*multix+ix)=float(ix+1)/float(multix)
        ppos(3,iy*multix+ix)=1.-float(iy)/float(multiy)
     endfor
   endfor
   for ifunc=0,nfunc-1 do begin
      if ndim eq 2 then begin
         aspect=(domainamr(3)-domainamr(1))/(domainamr(2)-domainamr(0))
      ;   pos=getpos(aspect,position=ppos(*,ifunc),margin=0.05)
      endif
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

      if strmid(plotmod,0,4) eq 'cont' then begin
         if abs(f_max - f_min) lt 1.0d-8 then begin
           f_max=f_max+1.0d0
           f_min=f_min-1.0d0
         endif
         levels=findgen(contourlevel)/contourlevel*(f_max-f_min)+f_min
      endif

      ; print,'contourlevel=',contourlevel
      ; print,'levels=',levels
      ; print,'fmin=',f_min
      ; print,'fmax=',f_max
    

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
	 'contour': contour,f>f_min,LEVELS=levels,/FOLLOW,$
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
         'velovect' :velovect,f1,f2,/NOERASE
         'ovelovect':velovect,f1,f2,/NOERASE,$
	    XRANGE=[0,n_elements(f1(*,0))-1],YRANGE=[0,n_elements(f1(0,*))-1]
         endcase
      'coord': case plotmod of
	 'volume'   :xvolume,f
	 'contour'  :contour,f>f_min,xx,yy,$
			LEVELS=levels,XSTYLE=1,YSTYLE=1,/NOERASE
	 'contlabel':contour,f>f_min,xx,yy,$
			LEVELS=levels,/FOLLOW,XSTYLE=1,YSTYLE=1,/NOERASE
	 'contfill' :contour,f>f_min,xx,yy, $
			LEVELS=levels,/FILL,XSTYLE=1,YSTYLE=1,/NOERASE
         'contfillB' : begin
                        device,decomposed=0
                        loadct,0
                        nlevel=n_elements(levels)
                        ncolor=nlevel+1
                        bottom=1
                        c_color=indgen(ncolor)+bottom
                        loadct,33,ncolors=ncolor,bottom=bottom
                        contour,f>f_min,xx,yy,LEVELS=levels,c_colors=c_color,$
                        position=pos,/FILL,XSTYLE=1,YSTYLE=1,/NOERASE
                        loadct,33
                        tkf=findgen(7)
                        tkf=(f_max-f_min)/6.*tkf+f_min
                        tka=strmid(strtrim(string(tkf),2),0,4)
                        colorbar,nlevels=100,position=[pos[2]+0.005,pos[1],pos[2]+0.02,pos[3]],/vertical,/right,$
                        ticknames=tka,color=180,xtitle='',ytitle='',title=''
                        loadct,0
                        hcc='pth=1 B=!9'+ string( "123b)+'!x2 !4 '+string( "142b)+'!x=0.5 rho!Iin!n=1 rho!Iout!n=10'
                        hcaa='pth=1 B=10 !4'+string( "142b)+'!x=0.01 rho!Iin!n=1 rho!Iout!n=10'
                        hcab='pth=1 B=100 !4'+string( "142b)+'!x=0.0001 rho!Iin!n=1 rho!Iout!n=10'
                        hcac='pth=1 B=0.1 !4'+string( "142b)+'!x=100 rho!Iin!n=1 rho!Iout!n=10'
                        hcay='pth=1 B=1 !4'+string( "142b)+'!x=1 rho!Iin!n=1 rho!Iout!n=10'
                        hcaw='pth=1 B=0.5 !4'+string( "142b)+'!x=4 rho!Iin!n=1 rho!Iout!n=10'
                        ;xyouts,0.34,0.49,hcac,/normal,charsize=2
                        b1=animfunc(x,w,'b1','mhd22',eqpars,wnames)
                        b2=animfunc(x,w,'b2','mhd22',eqpars,wnames)
                        siz=size(b1)
                        nxx=siz[1]
                        nyy=siz[2]
                        Bfield=fltarr(2,nxx,nyy)
                        Bfield[0,*,*]=b1
                        Bfield[1,*,*]=b2
                        nsx=5
                        nsy=nsx
                        seeds=intarr(2,nsx+nsy)
                        seeds[1,0:nsx-1]=0
                        ;seeds[1,0:nsx-1]=nyy-1
                        seeds[0,0:nsx-1]=long(findgen(nsx)*nxx/nsx)
                        seeds[0,nsx:nsx+nsy-1]=0
                        ;seeds[0,nsx:nsx+nsy-1]=nxx-1
                        seeds[1,nsx:nsx+nsy-1]=long(findgen(nsy)*nyy/nsy)
                        ;nsy=8
                        ;seeds=intarr(2,nsy)
                        ;seeds[0,0:nsy-1]=0
                        ;seeds[1,0:nsy-1]=indgen(nsy)*nyy/(nsy-1)                        
                        PARTICLE_TRACE,Bfield, seeds, verts, conn, MAX_ITERATIONS=3500
                        i = 0
                        sz = SIZE(verts, /STRUCTURE)
                        WHILE (i LT sz.dimensions[1]) DO BEGIN
                          nverts = conn[i]
                          PLOTS, x[verts[0, conn[i+1:i+nverts]],0], $
                                 x[verts[1, conn[i+1:i+nverts]],1]
                          i=i+nverts+1
                        ENDWHILE
                       end

         ; ***AMR: need overall domain size held fixed while overplotting
         ;         individual grids
         'contour_amr' : begin
                   ; to ensure the right plot extent for first frame: 
                   ;  use contour without data
                   contour,f>f_min,xx,yy, $
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        /NODATA,XSTYLE=1,YSTYLE=1,/NOERASE
                   ; overplot the grid extent from level lbase upward
                   lbase=2
                   if (info1amr ge lbase) then begin
                   polyfill,[corner1amr(0),corner1amr(0), $
                             corner1amr(2),corner1amr(2),corner1amr(0)],$
                            [corner1amr(1),corner1amr(3), $
                             corner1amr(3),corner1amr(1),corner1amr(1)], $
                            col=(1-info1amr/10.0)*!d.n_colors
                   endif
                   contour,f>f_min,xx,yy, $
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        LEVELS=levels,XSTYLE=1,YSTYLE=1,/NOERASE
                   if (info1amr ge lbase) then begin
                        ; overplot the grid extent
                        plots,corner1amr(0),corner1amr(1)
                        plots,corner1amr(0),corner1amr(3),/continue
                        plots,corner1amr(2),corner1amr(3),/continue
                        plots,corner1amr(2),corner1amr(1),/continue
                        plots,corner1amr(0),corner1amr(1),/continue
                   endif
                   end
         'contour2_amr' : begin
                   ; to ensure the right plot extent for first frame: 
                   ;  use contour without data
                   contour,f>f_min,xx,yy, $
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        /NODATA,XSTYLE=1,YSTYLE=1,/NOERASE
                   ; overplot the grid extent up to level lbase
                   lbase=2
                   if (info1amr ge lbase) then begin
                   polyfill,[corner1amr(0),corner1amr(0), $
                             corner1amr(2),corner1amr(2),corner1amr(0)],$
                            [corner1amr(1),corner1amr(3), $
                             corner1amr(3),corner1amr(1),corner1amr(1)], $
                            col=(info1amr/6.0)*!d.n_colors
                   ;         col=(1-info1amr/10.0)*!d.n_colors
                   endif
                   contour,f>f_min,xx,yy, $
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        LEVELS=levels,XSTYLE=1,YSTYLE=1,/NOERASE
                   end
         'contfilltest_amr' : begin
                   contour,f>f_min,xx,yy, $
                        XRANGE=[0.0,4.0], $
                        YRANGE=[0.0,4.0], $
                        LEVELS=levels,/FILL,XSTYLE=1,YSTYLE=1,/NOERASE
                   end
         'contfill_amr' : begin
                   contour,f>f_min,xx,yy, $
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        LEVELS=levels,/FILL,XSTYLE=1,YSTYLE=1,/NOERASE
                        ; overplot the grid extent
                        plots,corner1amr(0),corner1amr(1)
                        plots,corner1amr(0),corner1amr(3),/continue
                        plots,corner1amr(2),corner1amr(3),/continue
                        plots,corner1amr(2),corner1amr(1),/continue
                        plots,corner1amr(0),corner1amr(1),/continue
                   end
         'contfillT_amr' : begin
                   lth=6
                   contour,f>f_min,xx,yy, $
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        LEVELS=levels,/FILL,XSTYLE=1,YSTYLE=1,/NOERASE
                        ; overplot the grid extent
                        plots,corner1amr(0),corner1amr(1)
                        plots,corner1amr(0),corner1amr(3),/continue,thick=lth
                        plots,corner1amr(2),corner1amr(3),/continue,thick=lth
                        plots,corner1amr(2),corner1amr(1),/continue,thick=lth
                        plots,corner1amr(0),corner1amr(1),/continue,thick=lth
                   end
         'contfilln_amr' : begin
                   contour,f>f_min,xx,yy, $
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        LEVELS=levels,/FILL,XSTYLE=1,YSTYLE=1,/NOERASE
                   end
         'contfilln3_amr' : begin
                   lbase=1
                   lfine=2
                   if not(info1amr eq lbase or info1amr ge lfine) then return
                   contour,f>f_min,xx,yy, $
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        LEVELS=levels,/FILL,XSTYLE=1,YSTYLE=1,/NOERASE
                   end
         'contfilln2_amr' : begin
                   lbase=4
                   contour,f>f_min,xx,yy, $
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        LEVELS=levels,/FILL,XSTYLE=1,YSTYLE=1,/NOERASE
                   if (info1amr eq lbase) then begin
                        ; overplot the grid extent
                        plots,corner1amr(0),corner1amr(1)
                        plots,corner1amr(0),corner1amr(3),/continue
                        plots,corner1amr(2),corner1amr(3),/continue
                        plots,corner1amr(2),corner1amr(1),/continue
                        plots,corner1amr(0),corner1amr(1),/continue
                   endif
                   end
         'contfillN_amr' : begin
                   contour,f>f_min,xx,yy, $
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        LEVELS=levels,/FILL,XSTYLE=1,YSTYLE=1,/NOERASE
                   end
         'contfillNi_amr' : begin
                   contour,f>f_min,xx,yy, $
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        LEVELS=levels,/FILL,XSTYLE=1,YSTYLE=1,/NOERASE,/isotropic
                   end
         'contfillDD_amr' : begin
                   contour,f>f_min,xx,yy, $
                        XRANGE=[-domainamr(2),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        LEVELS=levels,/FILL,XSTYLE=1,YSTYLE=5,/NOERASE
                   contour,f>f_min,-xx,yy, $
                        XRANGE=[-domainamr(2),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        LEVELS=levels,/FILL,XSTYLE=5,YSTYLE=1,/NOERASE
                   end
         'contfill4D_amr' : begin
                   contour,f>f_min,xx,yy, $
                        XRANGE=[-domainamr(2),domainamr(2)], $
                        YRANGE=[-domainamr(3),domainamr(3)], $
                        LEVELS=levels,/FILL,XSTYLE=5,YSTYLE=5,/NOERASE
                   contour,f>f_min,-xx,yy, $
                        XRANGE=[-domainamr(2),domainamr(2)], $
                        YRANGE=[-domainamr(3),domainamr(3)], $
                        LEVELS=levels,/FILL,XSTYLE=5,YSTYLE=1,/NOERASE
                   contour,f>f_min,-xx,-yy, $
                        XRANGE=[-domainamr(2),domainamr(2)], $
                        YRANGE=[-domainamr(3),domainamr(3)], $
                        LEVELS=levels,/FILL,XSTYLE=1,YSTYLE=5,/NOERASE
                   contour,f>f_min,xx,-yy, $
                        XRANGE=[-domainamr(2),domainamr(2)], $
                        YRANGE=[-domainamr(3),domainamr(3)], $
                        LEVELS=levels,/FILL,XSTYLE=5,YSTYLE=5,/NOERASE
                   end
         'contfillNf_amr' : begin
                   lbase=1
                   lfine=4
                   if not(info1amr eq lbase or info1amr ge lfine) then return
                   contour,f>f_min,xx,yy, $
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        LEVELS=levels,/FILL,XSTYLE=5,YSTYLE=5,/NOERASE
                   end
         'contfillNfy_amr' : begin
                   lbase=1
                   lfine=4
                   if not(info1amr eq lbase or info1amr ge lfine) then return
                   contour,f>f_min,xx,yy, $
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        ycharsize=0.8,$
                        LEVELS=levels,/FILL,XSTYLE=5,YSTYLE=9,/NOERASE
                   end
         'contfillNfxy_amr' : begin
                   lbase=1
                   lfine=4
                   if not(info1amr eq lbase or info1amr ge lfine) then return
                   contour,f>f_min,xx,yy, $
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        ycharsize=0.8,$
                        LEVELS=levels,/FILL,XSTYLE=9,YSTYLE=9,/NOERASE
                   end
         'contourfy_amr' : begin
                   lbase=1
                   lfine=4
                   if not(info1amr eq lbase or info1amr ge lfine) then return
                   contour,f>f_min,xx,yy, $
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        /NODATA,LEVELS=levels,XSTYLE=5,YSTYLE=9,/NOERASE
                        ; overplot the grid extent
                        plots,corner1amr(0),corner1amr(1)
                        plots,corner1amr(0),corner1amr(3),/continue
                        plots,corner1amr(2),corner1amr(3),/continue
                        plots,corner1amr(2),corner1amr(1),/continue
                        plots,corner1amr(0),corner1amr(1),/continue
                   if (info1amr ge lfine) then begin
                        polyfill,[corner1amr(0),corner1amr(0), $
                             corner1amr(2),corner1amr(2),corner1amr(0)],$
                            [corner1amr(1),corner1amr(3), $
                             corner1amr(3),corner1amr(1),corner1amr(1)], $
                            col=(info1amr/5.0)*!d.n_colors
                   endif
                   contour,f>f_min,xx,yy, $
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        LEVELS=levels,XSTYLE=5,YSTYLE=9,/NOERASE
                   end
         'contourfxy_amr' : begin
                   lbase=1
                   lfine=4
                   if not(info1amr eq lbase or info1amr ge lfine) then return
                   contour,f>f_min,xx,yy, $
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        /NODATA,LEVELS=levels,XSTYLE=9,YSTYLE=9,/NOERASE
                        ; overplot the grid extent
                        plots,corner1amr(0),corner1amr(1)
                        plots,corner1amr(0),corner1amr(3),/continue
                        plots,corner1amr(2),corner1amr(3),/continue
                        plots,corner1amr(2),corner1amr(1),/continue
                        plots,corner1amr(0),corner1amr(1),/continue
                   if (info1amr ge lfine) then begin
                        polyfill,[corner1amr(0),corner1amr(0), $
                             corner1amr(2),corner1amr(2),corner1amr(0)],$
                            [corner1amr(1),corner1amr(3), $
                             corner1amr(3),corner1amr(1),corner1amr(1)], $
                            col=(info1amr/5.0)*!d.n_colors
                   endif
                   contour,f>f_min,xx,yy, $
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        LEVELS=levels,XSTYLE=9,YSTYLE=9,/NOERASE
                   end
         'contfillfxy_amr' : begin
                   lbase=1
                   lfine=4
                   if not(info1amr eq lbase or info1amr ge lfine) then return
                   contour,f>f_min,xx,yy, $
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        LEVELS=levels,/FILL,XSTYLE=9,YSTYLE=9,/NOERASE
                        ; overplot the grid extent
                        plots,corner1amr(0),corner1amr(1)
                        plots,corner1amr(0),corner1amr(3),/continue
                        plots,corner1amr(2),corner1amr(3),/continue
                        plots,corner1amr(2),corner1amr(1),/continue
                        plots,corner1amr(0),corner1amr(1),/continue
                   end
         'contfillfy_amr' : begin
                   lbase=1
                   lfine=4
                   if not(info1amr eq lbase or info1amr ge lfine) then return
                   contour,f>f_min,xx,yy, $
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        LEVELS=levels,/FILL,XSTYLE=5,YSTYLE=9,/NOERASE
                        ; overplot the grid extent
                        plots,corner1amr(0),corner1amr(1)
                        plots,corner1amr(0),corner1amr(3),/continue
                        plots,corner1amr(2),corner1amr(3),/continue
                        plots,corner1amr(2),corner1amr(1),/continue
                        plots,corner1amr(0),corner1amr(1),/continue
                   end
         'contfillf_amr' : begin
                   lbase=1
                   lfine=4
                   if not(info1amr eq lbase or info1amr ge lfine) then return
                   contour,f>f_min,xx,yy, $
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)], $
                        LEVELS=levels,/FILL,XSTYLE=1,YSTYLE=1,/NOERASE
                        ; overplot the grid extent
                        plots,corner1amr(0),corner1amr(1)
                        plots,corner1amr(0),corner1amr(3),/continue
                        plots,corner1amr(2),corner1amr(3),/continue
                        plots,corner1amr(2),corner1amr(1),/continue
                        plots,corner1amr(0),corner1amr(1),/continue
                   end
	 'plot'    :plot,xx,f,YRANGE=[f_min,f_max],XSTYLE=18,YSTYLE=18,$
						   LINE=linestyle,/NOERASE
         ; ***AMR: need x-range fixed
	 'plot1_amr'    :begin
                   plot,xx,f,XRANGE=[domainamr(0),domainamr(1)], $
                        YRANGE=[f_min,f_max],XSTYLE=18,YSTYLE=18,$
					   LINE=linestyle,/NOERASE
                   end
	 'plot_amr'    :begin
                   ;print,'Domain=',domainamr(*)
                   ;print,'Corner=',corner1amr(*)
                   plot,xx,f,XRANGE=[domainamr(0),domainamr(1)], $
                        YRANGE=[f_min,f_max],XSTYLE=18,YSTYLE=18,$
					   LINE=linestyle,/NOERASE
                   plot,[corner1amr(0),corner1amr(0)],[f_min,f_max],$
                        line=0,psym=0,XRANGE=[domainamr(0),domainamr(1)], $
                        YRANGE=[f_min,f_max],XSTYLE=18,YSTYLE=18,/NOERASE
                   plot,[corner1amr(1),corner1amr(1)],[f_min,f_max],$
                        line=0,psym=0,XRANGE=[domainamr(0),domainamr(1)], $
                        YRANGE=[f_min,f_max],XSTYLE=18,YSTYLE=18,/NOERASE
                   end
	 'plots_amr'    :begin
                   eps=(f_max-f_min)*0.1
                   plot,xx,f,XRANGE=[domainamr(0),domainamr(1)], $
                        YRANGE=[f_min-eps,f_max+eps],XSTYLE=1,YSTYLE=17,$
					   LINE=linestyle,/NOERASE

                   polyfill,[corner1amr(0),corner1amr(0), $
                             corner1amr(1),corner1amr(1)],$
                            [f_min-eps,f_max+eps,f_max+eps,f_min-eps], $
                            col=(1-info1amr/10.0)*!d.n_colors
                   
                   plot,xx,f,XRANGE=[domainamr(0),domainamr(1)], $
                        YRANGE=[f_min-eps,f_max+eps],XSTYLE=1,YSTYLE=17,$
					   LINE=linestyle,/NOERASE

                   plot,[corner1amr(0),corner1amr(0)],[f_min,f_max],$
                        line=0,psym=0,XRANGE=[domainamr(0),domainamr(1)], $
                        YRANGE=[f_min-eps,f_max+eps],XSTYLE=1,YSTYLE=17,$
                        /NOERASE
                   plot,[corner1amr(1),corner1amr(1)],[f_min,f_max],$
                        line=0,psym=0,XRANGE=[domainamr(0),domainamr(1)], $
                        YRANGE=[f_min-eps,f_max+eps],XSTYLE=1,YSTYLE=17,$
                        /NOERASE
                   end
	 'plotn_amr'    :begin
                   eps=(f_max-f_min)*0.1
                   eps=0.0
                   plot,xx,f,XRANGE=[domainamr(0),domainamr(1)], $
                        YRANGE=[f_min-eps,f_max+eps],XSTYLE=1,YSTYLE=17,$
					   LINE=linestyle,/NOERASE

                   ;polyfill,[corner1amr(0),corner1amr(0), $
                   ;          corner1amr(1),corner1amr(1)],$
                   ;         [f_min-eps,f_max+eps,f_max+eps,f_min-eps], $
                   ;         col=(1-info1amr/10.0)*!d.n_colors
                   
                   ;plot,xx,f,XRANGE=[domainamr(0),domainamr(1)], $
                   ;     YRANGE=[f_min-eps,f_max+eps],XSTYLE=1,YSTYLE=17,$
                   ;     LINE=linestyle,/NOERASE

                   end
	 'plotf_amr'    :begin
                   eps=(f_max-f_min)*0.1
                   lmax=6.0
                   plot,xx,f,XRANGE=[domainamr(0),domainamr(1)], $
                        YRANGE=[f_min-2*eps,f_max+eps],XSTYLE=1,YSTYLE=17,$
                        SYMSIZE=0.3,LINE=linestyle,/NOERASE

                   polyfill,[corner1amr(0),corner1amr(0), $
                             corner1amr(1),corner1amr(1)],$
                            [f_min-2*eps,f_max+eps,f_max+eps,f_min-2*eps], $
                            ;col=((4+info1amr)/(lmax+5.0))*!d.n_colors
                            col=((4+lmax)/(lmax+5.0))*!d.n_colors
                   
                   plot,xx,f,XRANGE=[domainamr(0),domainamr(1)], $
                        YRANGE=[f_min-2*eps,f_max+eps],XSTYLE=1,YSTYLE=17,$
                        SYMSIZE=0.3,LINE=linestyle,/NOERASE

                    polyfill,[corner1amr(0),corner1amr(0), $
                             corner1amr(1),corner1amr(1)],$
                   [f_min-2*eps,f_min-2*eps+((info1amr-1)/lmax)*1.9*eps, $
                    f_min-2*eps+((info1amr-1)/lmax)*1.9*eps,f_min-2*eps], col=0 
                   end
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
         'vel_amr'    :begin
                   domainsize=(domainamr(2)-domainamr(0))* $
                              (domainamr(3)-domainamr(1))
                   gridsize=(corner1amr(2)-corner1amr(0))* $
                            (corner1amr(3)-corner1amr(1))
                   velvector=400*gridsize/domainsize
                   if(velvector le 10) then velvector=10
                   print,'gridsize,domainsize,velvector'
                   print,gridsize,domainsize,velvector
                   velamr, $
                   f1,f2,xx,yy,$
                        NVECS=velvector,MAXVAL=f_max,$
                        DYNAMIC=0,SEED=velseed,/NOERASE,$
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)]
                   end
         'velovect' :velovect,f1,f2,xx(*,0),yy(0,*),/NOERASE
         'ovelovect':velovect,f1,f2,xx(*,0),yy(0,*),/NOERASE,$
			XRANGE=[min(xx),max(xx)],YRANGE=[min(yy),max(yy)]
         'ovelovect_amr':velovect,f1,f2,xx(*,0),yy(0,*),/NOERASE,$
                        XRANGE=[domainamr(0),domainamr(2)], $
                        YRANGE=[domainamr(1),domainamr(3)]
         endcase
      else:print,'Unknown axistype:',axistype
      endcase
      !p.multi(0)=multi0
   endfor
end
;===========================================================================
pro putbottom,multix,multiy,ix,iy,ninfo,nx,it,time,timeunit

on_error,2

if ninfo lt 1 then return
info=''
if ninfo gt 2 then info='nx='+string(nx,format='(3(i4,","))')+' '
if ninfo gt 1 then info=info+'it='+string(it,format='(i10)')+', '
if n_elements(timeunit) eq 0 then begin
info=info+'time='+string(time,format='(g12.5)')
endif else begin
info=info+'time='+string(time,format='(g12.5)')+timeunit
endelse
xyouts,5+(ix*!d.x_size)/multix,8+(iy*!d.y_size)/multiy,/DEV,info ;;,FONT=1

end
;===========================================================================
pro putheader,multix,multiy,ix,iy,ninfo,headline,nx

on_error,2
   
if ninfo lt 1 then return
info=strtrim(headline,2)
if ninfo gt 1 then info=info+' (nx='+string(nx,format='(3(i4))')+')'
xyouts,5+(ix*!d.x_size)/multix,-12+((iy+1)*!d.y_size)/multiy,/DEV,info ;;,FONT=1

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

if n_elements(ndir) eq 0 or n_elements(f) eq 0 $
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
pro extrapolateamr,x,w,ndim,nw,nx,nxamr,corneramr,xamr,wamr,ngrids,nexp,islope,doask
;==========================================
  on_error,2
  print,'nx=',nx  
  print,'ndim=',ndim
  print,' give new dimensions (nx(ndim) array)'
  nx=lonarr(ndim)
  if n_elements(nexp) eq 0 then begin
     read,nx  
  endif else begin
     nx=nexp
  endelse
  ;print,'nx=',nx  
  getdomain_amr,ndim,domainamr,corneramr,0,0
  if n_elements(islope) eq 0 then islope=0
  asknum,'copy or linear reconstruction (0/1) (for 1D or 2D only)',islope,doask
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
    x=dblarr(nx(0),ndim)
    w=dblarr(nx(0),nw)
    delxf=(domainamr(1)-domainamr(0))/nx(0)
    slope1C=dblarr(nw)
    slope1R=dblarr(nw)
    slope1L=dblarr(nw)
    sign1C=dblarr(nw)
    slope1=dblarr(nw)
    for igrid=0,ngrids-1 do begin
      n1g=nxamr(0,igrid)
      delxl=(corneramr(1,igrid)-corneramr(0,igrid))/n1g
      cl0=delxl/delxf
      ;print,'cl0=',cl0,delxl,delxf
      icl0=long(cl0+0.0001)
      ;print,'icl0=',icl0
      nfine=n1g*cl0
      icounter=float(0.0)
      for x0=0,n1g-1  do begin
         rndex0=(corneramr(0,igrid)-domainamr(0))/delxf+x0*cl0
         ; cell center for coarse cell
         rccoarse=corneramr(0,igrid)+(x0+0.5)*delxl
         if islope eq 1 then begin
            x0R=min([n1g-1,x0+1])
            x0L=max([0,x0-1])
            slope1C(0:nw-1)=0.5*(wamr(x0R,0:nw-1,igrid)-wamr(x0L,0:nw-1,igrid))
            if (x0R-x0L) lt 1.5 then slope1C(0:nw-1)=slope1C(0:nw-1)*2.0
            slope1L(0:nw-1)=(wamr(x0,0:nw-1,igrid)-wamr(x0L,0:nw-1,igrid))
            if (x0-x0L) lt 0.5 then slope1L(0:nw-1)=slope1C(0:nw-1)
            slope1R(0:nw-1)=(wamr(x0R,0:nw-1,igrid)-wamr(x0,0:nw-1,igrid))
            if (x0R-x0) lt 0.5 then slope1R(0:nw-1)=slope1C(0:nw-1)
            for iw=0,nw-1 do begin
               if slope1C(iw) ge 0.0 then sign1C(iw)=1.0
               if slope1C(iw) lt 0.0 then sign1C(iw)=-1.0
               slope1L(iw)=sign1C(iw)*slope1L(iw)
               slope1R(iw)=sign1C(iw)*slope1R(iw)
               slope1(iw)=sign1C(iw)*max([0.0,min([abs(slope1C(iw)), $
                   slope1L(iw),slope1R(iw)])])
            endfor
         endif else begin
            for iw=0,nw-1 do begin
               slope1(iw)=0.0
            endfor
         endelse
         for x0f=0,icl0-1 do begin
            ; cell center for fine cell
            rcfine=domainamr(0)+(rndex0+x0f+0.5)*delxf
            eta1=(rcfine-rccoarse)/delxl
            index01=long(float(rndex0+x0f+0.00001))
            if abs((index01*delxf+0.5*delxf+domainamr(0))-rcfine) $
                  ge 0.000001 then begin
                  print,'rcfine,eta1,delxf,domainamr(0),index01'
                  print,rcfine,eta1,delxf,domainamr(0),index01
                  stop
            endif
            x(index01,0)= rcfine
            w(index01,0:nw-1)=wamr(x0,0:nw-1,igrid)+slope1(0:nw-1)*eta1
            icounter=icounter+1.0
         endfor
      endfor
      if icounter ne float(nfine) then begin
         print,'Warning:',icounter,nfine, $
               ' while ',float(nfine)
         stop
      endif
    endfor
  end
  ;-------------- 2D ----------------------
  2: begin
    x=dblarr(nx(0),nx(1),ndim)
    w=dblarr(nx(0),nx(1),nw)
    delx1f=(domainamr(2)-domainamr(0))/nx(0)
    ;print,'del_x1_fine=',delx1f
    delx2f=(domainamr(3)-domainamr(1))/nx(1)
    ;print,'del_x2_fine=',delx2f
    slope1C=dblarr(nw)
    slope1R=dblarr(nw)
    slope1L=dblarr(nw)
    sign1C=dblarr(nw)
    slope1=dblarr(nw)
    slope2C=dblarr(nw)
    slope2R=dblarr(nw)
    slope2L=dblarr(nw)
    sign2C=dblarr(nw)
    slope2=dblarr(nw)
    for igrid=0,ngrids-1 do begin
      ;print,'igrid=',igrid
      n1g=nxamr(0,igrid)
      n2g=nxamr(1,igrid)
      ;print,'n1g,n2g=',n1g,n2g
      delx1l=(corneramr(2,igrid)-corneramr(0,igrid))/n1g
      ;print,'del_x1_igrid=',delx1l
      delx2l=(corneramr(3,igrid)-corneramr(1,igrid))/n2g
      ;print,'del_x2_igrid=',delx2l
      cl0=delx1l/delx1f
      ;print,'cl0=',cl0,delx1l,delx1f
      icl0=long(cl0+0.0001)
      ;print,'icl0=',icl0
      nfine1=n1g*cl0
      nfine2=n2g*cl0
      icounter=float(0.0)
      for x01=0,n1g-1  do begin
         rndex01=(corneramr(0,igrid)-domainamr(0))/delx1f+x01*cl0
         ; cell center for coarse cell
         rc1coarse=corneramr(0,igrid)+(x01+0.5)*delx1l
         x01R=min([n1g-1,x01+1])
         x01L=max([0,x01-1])
         for x02=0,n2g-1  do begin
            rndex02=(corneramr(1,igrid)-domainamr(1))/delx2f+x02*cl0
            ; cell center for coarse cell
            rc2coarse=corneramr(1,igrid)+(x02+0.5)*delx2l
            x02R=min([n2g-1,x02+1])
            x02L=max([0,x02-1])
            if islope eq 1 then begin
               slope1C(0:nw-1)=  $
                 0.5*(wamr(x01R,x02,0:nw-1,igrid)-wamr(x01L,x02,0:nw-1,igrid))
               if (x01R-x01L) lt 1.5 then slope1C(0:nw-1)=slope1C(0:nw-1)*2.0
               slope1L(0:nw-1)=  $
                (wamr(x01,x02,0:nw-1,igrid)-wamr(x01L,x02,0:nw-1,igrid))
               if (x01-x01L) lt 0.5 then slope1L(0:nw-1)=slope1C(0:nw-1)
               slope1R(0:nw-1)=  $
                (wamr(x01R,x02,0:nw-1,igrid)-wamr(x01,x02,0:nw-1,igrid))
               if (x01R-x01) lt 0.5 then slope1R(0:nw-1)=slope1C(0:nw-1)
               slope2C(0:nw-1)=  $
                 0.5*(wamr(x01,x02R,0:nw-1,igrid)-wamr(x01,x02L,0:nw-1,igrid))
               if (x02R-x02L) lt 1.5 then slope2C(0:nw-1)=slope2C(0:nw-1)*2.0
               slope2L(0:nw-1)=  $
                (wamr(x01,x02,0:nw-1,igrid)-wamr(x01,x02L,0:nw-1,igrid))
               if (x02-x02L) lt 0.5 then slope2L(0:nw-1)=slope2C(0:nw-1)
               slope2R(0:nw-1)=  $
                (wamr(x01,x02R,0:nw-1,igrid)-wamr(x01,x02,0:nw-1,igrid))
               if (x02R-x02) lt 0.5 then slope2R(0:nw-1)=slope2C(0:nw-1)
               for iw=0,nw-1 do begin
                  if slope1C(iw) ge 0.0 then sign1C(iw)=1.0
                  if slope1C(iw) lt 0.0 then sign1C(iw)=-1.0
                  slope1L(iw)=sign1C(iw)*slope1L(iw)
                  slope1R(iw)=sign1C(iw)*slope1R(iw)
                  slope1(iw)=sign1C(iw)*max([0.0,min([abs(slope1C(iw)), $
                             slope1L(iw),slope1R(iw)])])
                  if slope2C(iw) ge 0.0 then sign2C(iw)=1.0
                  if slope2C(iw) lt 0.0 then sign2C(iw)=-1.0
                  slope2L(iw)=sign2C(iw)*slope2L(iw)
                  slope2R(iw)=sign2C(iw)*slope2R(iw)
                  slope2(iw)=sign2C(iw)*max([0.0,min([abs(slope2C(iw)), $
                             slope2L(iw),slope2R(iw)])])
               endfor
            endif else begin
               for iw=0,nw-1 do begin
                  slope1(iw)=0.0
                  slope2(iw)=0.0
               endfor
            endelse
            for x0f1=0,icl0-1 do begin
               ; cell center for fine cell
               rc1fine=domainamr(0)+(rndex01+x0f1+0.5)*delx1f
               eta1=(rc1fine-rc1coarse)/delx1l
               index01=long(float(rndex01+x0f1+0.00001))
               if abs((index01*delx1f+0.5*delx1f+domainamr(0))-rc1fine) $
                  ge 0.000001 then begin
                  print,'rc1fine,eta1,delx1f,domainamr(0),index01'
                  print,rc1fine,eta1,delx1f,domainamr(0),index01
                  stop
               endif
               x(index01,0:nx(1)-1,0)= rc1fine
            for x0f2=0,icl0-1 do begin
               ; cell center for fine cell
               rc2fine=domainamr(1)+(rndex02+x0f2+0.5)*delx2f
               eta2=(rc2fine-rc2coarse)/delx2l
               index02=long(float(rndex02+x0f2+0.00001))
               if abs((index02*delx2f+0.5*delx2f+domainamr(1))-rc2fine) $
                  ge 0.000001 then begin
                  print,'rc2fine,eta2,delx2f,domainamr(1),index02'
                  print,rc2fine,eta2,delx2f,domainamr(1),index02
                  stop
               endif
               x(0:nx(0)-1,index02,1)= rc2fine
               w(index01,index02,0:nw-1)= $
                 wamr(x01,x02,0:nw-1,igrid)+slope1(0:nw-1)*eta1 $
                                           +slope2(0:nw-1)*eta2
               icounter=icounter+1.0
            endfor
            endfor
         endfor
      endfor
      if icounter ne float(nfine1*nfine2) then begin
         print,'Warning:',icounter,nfine1*nfine2, $
               ' while ',float(nfine1*nfine2)
         stop
      endif
    endfor
  end
  ;-------------- 3D ----------------------
  3: begin
    x=dblarr(nx(0),nx(1),nx(2),ndim)
    w=dblarr(nx(0),nx(1),nx(2),nw)
    delx1f=(domainamr(3)-domainamr(0))/nx(0)
    delx2f=(domainamr(4)-domainamr(1))/nx(1)
    delx3f=(domainamr(5)-domainamr(2))/nx(2)
    for igrid=0,ngrids-1 do begin
      n1g=nxamr(0,igrid)
      n2g=nxamr(1,igrid)
      n3g=nxamr(2,igrid)
      delx1l=(corneramr(3,igrid)-corneramr(0,igrid))/n1g
      delx2l=(corneramr(4,igrid)-corneramr(1,igrid))/n2g
      delx3l=(corneramr(5,igrid)-corneramr(2,igrid))/n3g
      cl0=delx1l/delx1f
      print,'cl0=',cl0,delx1l,delx1f
      icl0=long(cl0+0.0001)
      print,'icl0=',icl0
      nfine1=n1g*cl0
      nfine2=n2g*cl0
      nfine3=n3g*cl0
      icounter=float(0.0)
      for x01=0,n1g-1  do begin
         rndex01=(corneramr(0,igrid)-domainamr(0))/delx1f+x01*cl0
         for x02=0,n2g-1  do begin
            rndex02=(corneramr(1,igrid)-domainamr(1))/delx2f+x02*cl0
            for x03=0,n3g-1  do begin
               rndex03=(corneramr(2,igrid)-domainamr(2))/delx3f+x03*cl0
               for x0f1=0,icl0-1 do begin
               x(rndex01+x0f1,0:nx(1)-1,0:nx(2)-1,0)= $
                   domainamr(0)+(rndex01+x0f1+0.5)*delx1f
               for x0f2=0,icl0-1 do begin
               x(0:nx(0)-1,rndex02+x0f2,0:nx(2)-1,1)= $
                   domainamr(1)+(rndex02+x0f2+0.5)*delx2f
               for x0f3=0,icl0-1 do begin
               x(0:nx(0)-1,0:nx(1)-1,rndex03+x0f3,2)= $
                   domainamr(2)+(rndex03+x0f3+0.5)*delx3f
               w(rndex01+x0f1,rndex02+x0f2,rndex03+x0f3,0:nw-1)= $
                  wamr(x01,x02,x03,0:nw-1,igrid)
               icounter=icounter+1.0
               endfor
               endfor
               endfor
            endfor
         endfor
      endfor
      if icounter ne float(nfine1*nfine2*nfine3) then stop
    endfor
  end
  endcase

end
;==========================================
pro extrapolateamr3d,x,w,ndim,nw,nx,nxamr,corneramr,xamr,wamr,ngrids
;==========================================
  on_error,2
  print,'nx=',nx  
  print,'ndim=',ndim
  if ndim ne 3 then stop
  print,' give new dimensions (nx(ndim) array)'
  nx=lonarr(ndim)
  read,nx  
  print,'nx=',nx  
  getdomain_amr,ndim,domainamr,corneramr,0,0
  
;   x=dblarr(nx(0),nx(1),nx(2),ndim)
    x=dblarr(2,2,2,ndim)
    w=dblarr(nx(0),nx(1),nx(2),1)
    delx1f=(domainamr(3)-domainamr(0))/nx(0)
    delx2f=(domainamr(4)-domainamr(1))/nx(1)
    delx3f=(domainamr(5)-domainamr(2))/nx(2)
    for igrid=0,ngrids-1 do begin
      n1g=nxamr(0,igrid)
      n2g=nxamr(1,igrid)
      n3g=nxamr(2,igrid)
      delx1l=(corneramr(3,igrid)-corneramr(0,igrid))/n1g
      delx2l=(corneramr(4,igrid)-corneramr(1,igrid))/n2g
      delx3l=(corneramr(5,igrid)-corneramr(2,igrid))/n3g
      cl0=delx1l/delx1f
      print,'cl0=',cl0,delx1l,delx1f
      icl0=long(cl0+0.0001)
      print,'icl0=',icl0
      nfine1=n1g*cl0
      nfine2=n2g*cl0
      nfine3=n3g*cl0
      icounter=float(0.0)
      for x01=0,n1g-1  do begin
         rndex01=(corneramr(0,igrid)-domainamr(0))/delx1f+x01*cl0
         for x02=0,n2g-1  do begin
            rndex02=(corneramr(1,igrid)-domainamr(1))/delx2f+x02*cl0
            for x03=0,n3g-1  do begin
               rndex03=(corneramr(2,igrid)-domainamr(2))/delx3f+x03*cl0
               for x0f1=0,icl0-1 do begin
;               x(rndex01+x0f1,0:nx(1)-1,0:nx(2)-1,0)= $
;                   domainamr(0)+(rndex01+x0f1+0.5)*delx1f
               for x0f2=0,icl0-1 do begin
;               x(0:nx(0)-1,rndex02+x0f2,0:nx(2)-1,1)= $
;                   domainamr(1)+(rndex02+x0f2+0.5)*delx2f
               for x0f3=0,icl0-1 do begin
;               x(0:nx(0)-1,0:nx(1)-1,rndex03+x0f3,2)= $
;                   domainamr(2)+(rndex03+x0f3+0.5)*delx3f
;               w(rndex01+x0f1,rndex02+x0f2,rndex03+x0f3,0:nw-1)= $
;                  wamr(x01,x02,x03,0:nw-1,igrid)
               w(rndex01+x0f1,rndex02+x0f2,rndex03+x0f3,0)= $
                  wamr(x01,x02,x03,0,igrid)
               icounter=icounter+1.0
               endfor
               endfor
               endfor
            endfor
         endfor
      endfor
      if icounter ne float(nfine1*nfine2*nfine3) then stop
    endfor

end
;======================================================================================================================|
pro particle_trajectory, field_data, field_seed, verts, X = x, Y = y, Z = z,  $
  MAX_ITERATIONS = max_iterations, INTEGRATION = integration, STEP_SIZE = step_size
;======================================================================================================================|
;  Purpose
;    traces the path of a massless particle through a vector field (uniform or non-uniform) from a given starting point. 
;    (IDL procedure 'particle_trace' is not a good choice)
;
;  Syntax
;    particle_trajectory, data, seeds, verts [, X = x] [, Y = y] [, Z = z] [,  $
;      MAX_ITERATIONS = max_iterations] [, INTEGRATION = integration] [, STEP_SIZE = step_size]
;      
;  Inputs
;    field_data: three- or four-dimensional array that defines the vector field.
;      for a two-dimensional field: [2, nx, ny]
;      for a three-dimensional field: [3, nx, ny, nz]
;    field_seed: starting point
;    x: the x coordinate, if not, the device coordinate will be used
;    y: the y coordinate, if not, the device coordinate will be used
;    z: the z coordinate, if not, the device coordinate will be used
;    max_iterations: the maximum number of line segments of trajectory (default 200)
;    integration: 
;      0: 2nd order Runge-Kutta (fast, low precision, default)
;      1: 4th order Runge-Kutta (fast, high precision)
;    step_size = iteration step size
;
;  Outputs
;    verts: a array that will contain the output trajectory 
;    
;  Examples
;  
;  Notes
;
;  History
;    2010-08-01 written by R. L. Jiang at Kwasan Observatory
;
;  Bug report
;    rljiang@nju.edu.cn
;======================================================================================================================|
  on_error, 2

  if n_elements(field_data) eq 0 then message, 'no input vector field.'
  if n_elements(field_seed) eq 0 then message, 'no input starting point.'
  if not keyword_set(max_iterations) then max_iterations = 200l
  if not keyword_set(integration) then integration = 0
  if not keyword_set(step_size) then step_size = 1.0

  if size(field_data, /n_dimensions) eq 3 then begin
    nx = (size(field_data, /dimensions))[1]
    ny = (size(field_data, /dimensions))[2]
    data = fltarr(3, nx, ny, 2)
    data[0:1, *, *, 0] = field_data
    data[0:1, *, *, 1] = field_data
    data[2, *, *, *] = 0.0
    seed = [field_seed, 0.5]
  endif else begin
    data = field_data
    seed = field_seed
  endelse

;----------------------------------------------------------------------------------------------------------------------|
;  two or three dimensional
;----------------------------------------------------------------------------------------------------------------------|
  error = (machar()).eps
  nx = (size(data, /dimensions))[1]
  ny = (size(data, /dimensions))[2]
  nz = (size(data, /dimensions))[3]
  if n_elements(x) ne nx then x = findgen(nx)
  if n_elements(y) ne ny then y = findgen(ny)
  if n_elements(z) ne nz then z = findgen(nz)
  boundary = float([[min(x), max(x)], [min(y), max(y)], [min(z), max(z)]])
  verts = seed
  seed_x = value_locate(x, seed[0])
  seed_y = value_locate(y, seed[1])
  seed_z = value_locate(z, seed[2])
  if seed_x lt 0 or seed_x gt nx-1 or  $
     seed_y lt 0 or seed_y gt ny-1 or  $
     seed_z lt 0 or seed_z gt nz-1 then return
  step_f = step_size*min((boundary[1, *]-boundary[0, *])/[nx-1, ny-1, nz-1])
  out_of_range = 0

  for nstep = 1, max_iterations do begin
    step_s = step_f/max(abs(data[*, seed_x, seed_y, seed_z]))
    catch, out_of_range
    if out_of_range ne 0 then begin
      if out_of_range eq -161 or out_of_range eq -163 then begin
        catch, /cancel
        break
      endif
    endif

;----------------------------------------------------------------------------------------------------------------------|
;  Runge-Kutta 2nd
;----------------------------------------------------------------------------------------------------------------------|
    if integration eq 0 then begin
      while seed[0] lt x(seed_x) do seed_x -= 1 & while seed[0] ge x(seed_x+1) do seed_x += 1
      while seed[1] lt y(seed_y) do seed_y -= 1 & while seed[1] ge y(seed_y+1) do seed_y += 1
      while seed[2] lt z(seed_z) do seed_z -= 1 & while seed[2] ge z(seed_z+1) do seed_z += 1
      box = data[*, seed_x:seed_x+1, seed_y:seed_y+1, seed_z:seed_z+1]
      len = [[seed[0]-x(seed_x), x(seed_x+1)-seed[0]],  $
             [seed[1]-y(seed_y), y(seed_y+1)-seed[1]],  $
             [seed[2]-z(seed_z), z(seed_z+1)-seed[2]]]
      field_vector1 = fltarr(3) & for i = 0, 1 do for j = 0, 1 do for k = 0, 1 do  $
        field_vector1 += box[*, i, j, k]*len[1-i, 0]*len[1-j, 1]*len[1-k, 2]
      field_vector1 /= float(product(total(len, 1)))
      seed1 = seed + step_s*field_vector1

      while seed1[0] lt x(seed_x) do seed_x -= 1 & while seed1[0] ge x(seed_x+1) do seed_x += 1
      while seed1[1] lt y(seed_y) do seed_y -= 1 & while seed1[1] ge y(seed_y+1) do seed_y += 1
      while seed1[2] lt z(seed_z) do seed_z -= 1 & while seed1[2] ge z(seed_z+1) do seed_z += 1
      box = data[*, seed_x:seed_x+1, seed_y:seed_y+1, seed_z:seed_z+1]
      len = [[seed1[0]-x(seed_x), x(seed_x+1)-seed1[0]],  $
             [seed1[1]-y(seed_y), y(seed_y+1)-seed1[1]],  $
             [seed1[2]-z(seed_z), z(seed_z+1)-seed1[2]]]
      field_vector2 = fltarr(3) & for i = 0, 1 do for j = 0, 1 do for k = 0, 1 do  $
        field_vector2 += box[*, i, j, k]*len[1-i, 0]*len[1-j, 1]*len[1-k, 2]
      field_vector2 /= float(product(total(len, 1)))
      seed = seed + step_s*(field_vector1+field_vector2)/2.0
    endif else begin

;----------------------------------------------------------------------------------------------------------------------|
;  Runge-Kutta 4th
;----------------------------------------------------------------------------------------------------------------------|
      while seed[0] lt x(seed_x) do seed_x -= 1 & while seed[0] ge x(seed_x+1) do seed_x += 1
      while seed[1] lt y(seed_y) do seed_y -= 1 & while seed[1] ge y(seed_y+1) do seed_y += 1
      while seed[2] lt z(seed_z) do seed_z -= 1 & while seed[2] ge z(seed_z+1) do seed_z += 1
      box = data[*, seed_x:seed_x+1, seed_y:seed_y+1, seed_z:seed_z+1]
      len = [[seed[0]-x(seed_x), x(seed_x+1)-seed[0]],  $
             [seed[1]-y(seed_y), y(seed_y+1)-seed[1]],  $
             [seed[2]-z(seed_z), z(seed_z+1)-seed[2]]]
      field_vector1 = fltarr(3) & for i = 0, 1 do for j = 0, 1 do for k = 0, 1 do  $
        field_vector1 += box[*, i, j, k]*len[1-i, 0]*len[1-j, 1]*len[1-k, 2]
      field_vector1 /= float(product(total(len, 1)))
      seed1 = seed + step_s*field_vector1/2.0

      while seed1[0] lt x(seed_x) do seed_x -= 1 & while seed1[0] ge x(seed_x+1) do seed_x += 1
      while seed1[1] lt y(seed_y) do seed_y -= 1 & while seed1[1] ge y(seed_y+1) do seed_y += 1
      while seed1[2] lt z(seed_z) do seed_z -= 1 & while seed1[2] ge z(seed_z+1) do seed_z += 1
      box = data[*, seed_x:seed_x+1, seed_y:seed_y+1, seed_z:seed_z+1]
      len = [[seed1[0]-x(seed_x), x(seed_x+1)-seed1[0]],  $
             [seed1[1]-y(seed_y), y(seed_y+1)-seed1[1]],  $
             [seed1[2]-z(seed_z), z(seed_z+1)-seed1[2]]]
      field_vector2 = fltarr(3) & for i = 0, 1 do for j = 0, 1 do for k = 0, 1 do  $
        field_vector2 += box[*, i, j, k]*len[1-i, 0]*len[1-j, 1]*len[1-k, 2]
      field_vector2 /= float(product(total(len, 1)))
      seed2 = seed + step_s*field_vector2/2.0

      while seed2[0] lt x(seed_x) do seed_x -= 1 & while seed2[0] ge x(seed_x+1) do seed_x += 1
      while seed2[1] lt y(seed_y) do seed_y -= 1 & while seed2[1] ge y(seed_y+1) do seed_y += 1
      while seed2[2] lt z(seed_z) do seed_z -= 1 & while seed2[2] ge z(seed_z+1) do seed_z += 1
      box = data[*, seed_x:seed_x+1, seed_y:seed_y+1, seed_z:seed_z+1]
      len = [[seed2[0]-x(seed_x), x(seed_x+1)-seed2[0]],  $
             [seed2[1]-y(seed_y), y(seed_y+1)-seed2[1]],  $
             [seed2[2]-z(seed_z), z(seed_z+1)-seed2[2]]]
      field_vector3 = fltarr(3) & for i = 0, 1 do for j = 0, 1 do for k = 0, 1 do  $
        field_vector3 += box[*, i, j, k]*len[1-i, 0]*len[1-j, 1]*len[1-k, 2]
      field_vector3 /= float(product(total(len, 1)))
      seed3 = seed + step_s*field_vector3

      while seed3[0] lt x(seed_x) do seed_x -= 1 & while seed3[0] ge x(seed_x+1) do seed_x += 1
      while seed3[1] lt y(seed_y) do seed_y -= 1 & while seed3[1] ge y(seed_y+1) do seed_y += 1
      while seed3[2] lt z(seed_z) do seed_z -= 1 & while seed3[2] ge z(seed_z+1) do seed_z += 1
      box = data[*, seed_x:seed_x+1, seed_y:seed_y+1, seed_z:seed_z+1]
      len = [[seed3[0]-x(seed_x), x(seed_x+1)-seed3[0]],  $
             [seed3[1]-y(seed_y), y(seed_y+1)-seed3[1]],  $
             [seed3[2]-z(seed_z), z(seed_z+1)-seed3[2]]]
      field_vector4 = fltarr(3) & for i = 0, 1 do for j = 0, 1 do for k = 0, 1 do  $
        field_vector4 += box[*, i, j, k]*len[1-i, 0]*len[1-j, 1]*len[1-k, 2]
      field_vector4 /= float(product(total(len, 1)))
      seed = seed + step_s*(field_vector1+2.0*field_vector2+2.0*field_vector3+field_vector4)/6.0
    endelse

    verts = [[verts], [seed]]
  endfor

  if size(field_data, /n_dimensions) eq 3 then begin
    verts = verts[0:1, *]
  endif

end
;===================================================================================
FUNCTION GETPOS,ASPECT,POSITION=POSITION,MARGIN=MARGIN
;===================================================================================
;purpose: to give position keywords of a fixed ratio of vertical and horizontal axis
; in a given window

;Check arguments
if(n_params() ne 1) then message,'Usage:RESULT=GETPOS(ASPECT)'
if(n_elements(aspect) eq 0) then message,'Argument ASPECT is undefined'

;Check keywords
if(n_elements(position) eq 0) then position=[0.0,0.0,1.0,1.0]
if(n_elements(margin) eq 0) then margin=0.1

;Get range limited aspect ratio and margin input values
aspect_val=(float(aspect[0])>0.01)<100.0
margin_val=(float(margin[0])>0.0)<0.49
;Compute aspect ratio of position vector in this window
xsize=(position[2]-position[0])*!d.x_vsize
ysize=(position[3]-position[1])*!d.y_vsize
cur_aspect=ysize/xsize

;Compute aspect ratio of this window
win_aspect=float(!d.y_vsize)/float(!d.x_vsize)

;Compute height and width in normalized units
if(aspect_val ge cur_aspect) then begin
   height=(position[3]-position[1])-2.0*margin
   width=height*(win_aspect/aspect_val)
endif else begin
   width=(position[2]-position[0])-2.0*margin
   height=width*(aspect_val/win_aspect)
endelse

;Compute and return position vector
xcenter=0.5*(position[0]+position[2])
ycenter=0.5*(position[1]+position[3])
x0=xcenter-0.5*width
y0=ycenter-0.5*height
x1=xcenter+0.5*width
y1=ycenter+0.5*height
return,[x0,y0,x1,y1]

END

