;pro evol,color=coloron,ps=psout
coloron=1
psout=1
rebinset=1
fname='evolnT'
nexp=38400
iextrapol=1
func='T'
plotmode='plotn_amr'
   print,'======= CURRENT ANIMATION PARAMETERS ================'
   print,'firstpict=',firstpict,', dpict=',dpict,$
        ', npictmax=',npictmax,', savemovie (y/n)=',savemovie,$
        FORMAT='(a,i4,a,i4,a,i4,a,a)'
   print,'ax,az=',ax,',',az,', contourlevel=',contourlevel,$
         ', velvector=',velvector,', velspeed (0..5)=',velspeed,$
        FORMAT='(a,i4,a,i3,a,i3,a,i4,a,i2)'
   if keyword_set(multiplot) then begin
        siz=size(multiplot)
        if siz(0) eq 0 then multiplot=[multiplot,1,1]
        print,'multiplot= ',multiplot,', axistype (coord/cells)=',axistype,$
        FORMAT='(a,"[",i2,",",i2,",",i2,"]",a,a)'
   endif else $
        print,'multiplot= 0 (default), axistype (coord/cells)=',axistype,$
              FORMAT='(a,a)'
   print,'bottomline=',bottomline,', headerline=',headerline,$
        FORMAT='(a,i1,a,i1)'

   if keyword_set(cut) then help,cut
   if keyword_set(wsubtract) then help,wsubtract
   if keyword_set(velpos) then help,velpos
   velpos0 = velpos

   print,'======= FILE DESCRIPTION ============================'
   nfile=0
   if filename eq '' and logfilename ne '' then begin
      filename=logfilename
      while strpos(filename,'.log') ge 0 $
         do strput,filename,'.out',strpos(filename,'.log')
      askstr,'filename(s)   ',filename,1
   endif else $
      askstr,'filename(s)   ',filename,doask

   str2arr,filename,filenames,nfile
   gettype,filenames,filetypes,npictinfiles
   print,'filetype(s)   =',filetypes
   print,'npictinfile(s)=',npictinfiles

   ;====== OPEN FILE(S) AND READ AND PRINT HEADER(S)

   str2arr,physics,physicss,nfile
   physics=''
   anygencoord=0
   for ifile=0,nfile-1 do begin
      openfile,1,filenames(ifile),filetypes(ifile)
      phys=physicss(ifile)
      gethead,1,filetypes(ifile),phys, $
         headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables
      anygencoord=anygencoord or gencoord
      print,         'headline                  =',strtrim(headline,2)
      print,FORMAT='("variables                 =",100(a," "),$)',variables
      print,FORMAT='(" (ndim=",i2,", nw=",i2,")")',ndim,nw
      askstr,'physics (e.g. mhd12)      ',phys,doask
      physicss(ifile)=phys
      physics=physics + phys + ' '
   endfor

   print,'======= PLOTTING PARAMETERS ========================='
   ; ***AMR: skip transformation
   if (filetypes(0) eq 'ascii_amr') or  $
      (filetypes(0) eq 'bin_amr') then begin
      print,'No transformations allowed for AMR files'
      print,'setting transform to amr'
      transform='amr'
   endif
   if transform eq 'amr' then begin
     if n_elements(iextrapol) eq 0 then iextrapol=0
     asknum,'extrapolate to single grid (0/1/2)       ',iextrapol,doask
   endif
   help,noautodomain
   ; ***AMR: added argument to readplotpar
   readplotpar,ndim,cut,cut0,plotdim,nfunc,func,funcs,funcs1,funcs2,$
      plotmode,plotmodes,plottitle,plottitles,autorange,autoranges,doask,$
      transform,noautodomain,domainset
   print,nx
   readtransform,ndim,nx,anygencoord,transform,nxreg,wregpad,$
                 physicss(nfile-1),nvector,vectors,grid,doask

   print,'======= DETERMINE PLOTTING RANGES ==================='

   readlimits,nfunc,funcs,autoranges,noautorange,fmax,fmin,doask

   if noautorange then begin
      npict=(min(npictinfiles)-firstpict)/dpict+1
      if npict gt npictmax then npict=npictmax
      if npict lt 0 then npict=0
   endif else begin
      npict=0
      for ifile=0,nfile-1 do $
         openfile,ifile+1,filenames(ifile),filetypes(ifile)
      error=0
      while npict lt npictmax and not error do begin
         if npict eq 0 then nextpict=firstpict else nextpict=dpict
         for ifile=0,nfile-1 do begin
            ; ***AMR: extra arguments for getpict
            getpict,ifile+1,filetypes(ifile),nextpict,x,w,$
              headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,err,$ 
              nxamr,corneramr,infoamr,xamr,wamr,ngrids,0
            if keyword_set(wsubtract) then w=w-wsubtract
            wnames=variables(ndim:ndim+nw-1)
            usereg=(not gencoord and transform eq 'unpolar') or $
               (gencoord and (transform eq 'polar' or transform eq 'regular' $
                              or transform eq 'sphere'))
            error=err or error

            if not error then begin
               if usereg then case transform of
	          'regular':regulargrid,x_old,nxreg_old,x,xreg,nxreg,dxreg, $
                                     w,wreg,nw,w(0,0,*),triangles
		  'polar'  :polargrid  ,nvector,vectors,x,w,xreg,wreg
		  'sphere' :spheregrid  ,nvector,vectors,x,w,xreg,wreg
		  'unpolar':unpolargrid,nvector,vectors,x,w,xreg,wreg
	       endcase

               ; ***AMR: limts determined by looping over all grids
               if transform eq 'amr' then begin
                 erase
                 noerase=1
                 ;print,'animate: ngrids=',ngrids
                 for igrid=0,ngrids-1 do begin
                   getgrid_amr,ndim,nw,nx,x,w,nxamr,xamr,wamr,igrid
                   first= npict eq 0 and ifile eq 0 and igrid eq 0
                   getlimits,first,nfunc,funcs,funcs1,funcs2,autoranges, $
                      fmax,fmin,x,w,xreg,wreg,usereg,physicss(ifile), $
                      eqpar,wnames,cut0
                 endfor
               endif else begin
                 ; case of no amr
                 first= npict eq 0 and ifile eq 0
                 getlimits,first,nfunc,funcs,funcs1,funcs2,autoranges,$
                         fmax,fmin,x,w,xreg,wreg,usereg,physicss(ifile),$
                         eqpar,wnames,cut0
               endelse

               if ifile eq nfile-1 then begin
                  if npict eq 0 then print,FORMAT='("ipict:    ",$)'
                  print,FORMAT='(i4,$)',firstpict+npict*dpict
                  npict=npict+1
               endif
            endif
         endfor
      endwhile
      print
      for ifunc=0,nfunc-1 do $
      print,'Min and max value for ',funcs(ifunc),':',fmin(ifunc),fmax(ifunc)

   endelse
   ;print,'npict=',npict
   ;if transform eq 'amr' and npict eq 1 then print,'ngrids=',ngrids
   if npict eq 0 then begin
      print,'There are no frames to animate!'
      print,'   Check firstpict=',firstpict,' and dpict=',dpict
      if min(npictinfiles) lt firstpict then $
      print,'   The value of firstpict is larger than the minimum of' 
      print,'   npictinfiles=',npictinfiles
      retall
   endif

   help,corneramr
   getdomain_amr,ndim,domainamr,corneramr,noautodomain,domainset

   ;===== DO ANIMATION IN MULTIX * MULTIY MULTIPLE WINDOWS

   if keyword_set(multiplot) then begin
      multix=multiplot(0)
      multiy=multiplot(1)
      multidir=multiplot(2)
      npict1=(multix*multiy)/(nfunc*nfile)
      if npict1 eq 0 then npict1=1
   endif else if nfile eq 1 then begin
      multix=long(sqrt(nfunc-1)+1)
      multiy=long((nfunc-1)/multix+1)
      multidir=0
      npict1=1
   endif else begin
      multix=nfile
      multiy=nfunc
      multidir=1
      npict1=1
   endelse


   ipict=0
   ipict1=0
   iplot=0
   for ifile=0,nfile-1 do openfile,ifile+1,filenames(ifile),filetypes(ifile)
   error=0
   nHall=fltarr(nexp,npict)
   Teall=fltarr(nexp,npict)
   times=fltarr(npict)
   while ipict lt npict and not error do begin
      if ipict eq 0 then print,FORMAT='("ipict:    ",$)'
      print,FORMAT='(i4,$)',firstpict+ipict*dpict
      if ipict eq 0 then nextpict=firstpict else nextpict=dpict
      for ifile=0,nfile-1 do begin

         if npict gt 1 or nfile gt 1 or noautorange then begin
            ; ***AMR: added arguments to getpict
            getpict,ifile+1,filetypes(ifile),nextpict,x,w,$
               headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,err,$
               nxamr,corneramr,infoamr,xamr,wamr,ngrids,0
            if keyword_set(wsubtract) then w=w-wsubtract
            wnames=variables(ndim:ndim+nw-1)
            usereg=(not gencoord and transform eq 'unpolar') or $
               (gencoord and (transform eq 'polar' or transform eq 'regular' $
                              or transform eq 'sphere'))
            error=error or err
         endif

         if not error then begin
            if usereg then case transform of
		'regular':regulargrid,x_old,nxreg_old,x,xreg,nxreg,dxreg, $
				w,wreg,nw,wregpad,triangles
		'polar'  :begin
                            polargrid,nvector,vectors,x,w,xreg,wreg
                            variables(0:1)=['r','phi']
                          end
		'sphere' :begin
			    spheregrid  ,nvector,vectors,x,w,xreg,wreg
			    variables(0:2)=['r','theta','phi']
			  end
		'unpolar':begin
                            unpolargrid,nvector,vectors,x,w,xreg,wreg
                            variables(0:1)=['x','y']
                          end
            endcase

	    linestyle=0
            if multix*multiy lt nfunc*nfile then linestyle=ifile

            ; ***AMR: loop needed for plotting all grids
            if transform eq 'amr' then begin
               if iextrapol eq 1 then begin
                  extrapolateamr,x,w,ndim,nw,nx,nxamr,corneramr,xamr,wamr,$
                  ngrids,nexp,islope,doask
                 ; plot one extrapolated big grid
                 ;plotfunc,x,w,xreg,wreg,usereg,ndim,physicss(ifile),eqpar,$
                 ;variables,wnames,axistype,plotmodes,plottitles,$
                 ;ax,az,contourlevel,linestyle,$
                 ;velvector,velspeed,velseed,velpos,velx,vely,veltri,$
                 ;cut,cut0,plotdim,$
                 ;nfunc,multix,multiy,plotix,plotiy,funcs,funcs1,funcs2,fmin,fmax,f,domainamr
                 nH=animfunc(x,w,'n',physics,eqpar,wnames)
                 Te=animfunc(x,w,'T',physics,eqpar,wnames)
               endif else begin
                 for igrid=0,ngrids-1 do begin
                  getgrid_amr,ndim,nw,nx,x,w,nxamr,xamr,wamr,igrid
                  plotfunc,x,w,xreg,wreg,usereg,ndim,physicss(ifile),eqpar,$
                  variables,wnames,axistype,plotmodes,plottitles,$
                  ax,az,contourlevel,0,$
                  velvector,0,velseed,velpos,velx,vely,veltri,$
                  cut,cut0,plotdim,nfunc,multix,multiy,plotix,plotiy,$
                  funcs,funcs1,funcs2,fmin,fmax,f,$
                  domainamr,infoamr(igrid),corneramr(*,igrid)
                 endfor 
                endelse
                noerase=0
            endif else begin
              ; normal case of no amr
              plotfunc,x,w,xreg,wreg,usereg,ndim,physicss(ifile),eqpar,$
              variables,wnames,axistype,plotmodes,plottitles,$
              ax,az,contourlevel,linestyle,$
	      velvector,velspeed,velseed,velpos,velx,vely,veltri,$
              cut,cut0,plotdim,$
              nfunc,multix,multiy,plotix,plotiy,funcs,funcs1,funcs2,fmin,fmax,f
            endelse
            nHall[*,ipict]=nH
            Teall[*,ipict]=Te/1.e6
            times[ipict]=time/3600. 
         endif
      endfor

      ipict=ipict+1

   endwhile
   for ifile=1,nfile do close,ifile
; rebin 
if rebinset eq 1 then begin
  nx=nexp/200
  ny=npict
  dims=size(nHall,/dimensions)
  xloc=(findgen(nx)+0.5)*(dims[0]/float(nx))-0.5
  yloc=(findgen(ny)+0.5)*(dims[1]/float(ny))-0.5
  nHall=interpolate(nHall,xloc,yloc,/grid)
  Teall=interpolate(Teall,xloc,yloc,/grid)
  x=interpolate(x,xloc,/grid)
  times=interpolate(times,yloc,/grid)
endif
x=x*10.
; contour results
if psout eq 1 then begin
  xsize=18.0
  ysize=xsize
  charsizecb=1.
  charsizepl=1.
  symsz=0.12
  thic=5
  !P.FONT = 0
  entry_dev=!d.name
  set_plot,'ps'
  if coloron then begin
  device,SET_FONT='Times-Roman',filename=fname+'.eps',/encapsulated,xsize=xsize,ysize=ysize,/color,bits_per_pixel=8
  endif else begin
  device,SET_FONT='Times-Roman',filename=fname+'bw.eps',/encapsulated,xsize=xsize,ysize=ysize,bits_per_pixel=8
  endelse
endif else begin
  set_plot,'X'
  window,0,xsize=1400,ysize=1000
  charsizecb=1.3
  charsizepl=1.5
  symsz=0.4
  thic=5.
endelse
Temax=max(Teall)
Temin=min(Teall)
nHmax=max(nHall)
nHmin=min(nHall)
contourlevel=100
levels1=findgen(contourlevel)/contourlevel*(alog10(nHmax)-alog10(nHmin))+alog10(nHmin)
levels2=findgen(contourlevel)/contourlevel*(Temax-Temin)+Temin
nlevel=n_elements(levels1)
ncolor=nlevel+1
bottom=55
c_color=indgen(ncolor)*2.+bottom
; plot logn
if coloron then begin
  loadct,33,bottom=1
endif else begin
  loadct,0
  c_color=255-c_color
endelse
contour,alog10(nHall),x,times,position=[0.1,0.1,0.46,0.9],xtitle='!8s!x (Mm)',ytitle='Time (hr)',$
levels=levels1,c_colors=c_color,/FILL,XSTYLE=1,YSTYLE=1,charsize=charsizepl
;plot colorbar for logn
tkf=findgen(6)
tkf=(alog10(nHmax)-alog10(nHmin))/5.*tkf+alog10(nHmin)
tka=strmid(strtrim(string(tkf),2),0,4)
if coloron then begin
  colorbar,nlevels=100,ncolors=200,position=[0.1,0.91,0.46,0.92],ticknames=tka,color=0,/top,divisions=5,$
  charsize=charsizecb,bottom=bottom
endif else begin
  colorbar,nlevels=100,ncolors=200,position=[0.1,0.91,0.46,0.92],ticknames=tka,color=0,/top,divisions=5,$
  charsize=charsizecb,/invertcolors
endelse
; plot T
if coloron then begin
  loadct,33,bottom=1
endif else begin
  loadct,0
  c_color=255-c_color
endelse
contour,Teall,x,times,position=[0.54,0.1,0.9,0.9],xtitle='!8s!x (Mm)',ytitle='Time (hr)',$
levels=levels2,c_colors=c_color,/FILL,XSTYLE=1,YSTYLE=1,/NOERASE,charsize=charsizepl
;plot colorbar for T
;if coloron then loadct,33
tkf=findgen(6)
tkf=(Temax-Temin)/5.*tkf+Temin
tka=strmid(strtrim(string(tkf),2),0,4)
colorbar,nlevels=100,ncolors=255-bottom,position=[0.54,0.91,0.9,0.92],ticknames=tka,color=0,/top,divisions=5,$
charsize=charsizecb,bottom=bottom

if coloron then loadct,0
xyouts,0.28,0.96,'log !16n!x!dH!n (cm!u-3!n)',/normal,alignment=0.5,charsize=charsizepl+0.5
xyouts,0.72,0.96,'!16T!x (MK)',/normal,alignment=0.5,charsize=charsizepl+0.5
xyouts,0.5,0.975,'Case S2',/normal,alignment=0.5,charsize=charsizepl+0.5
if psout eq 0 then write_png,fname+'.png',tvrd(true=1)
if psout eq 1 then begin
  device,/close
  set_plot,entry_dev
endif

end
