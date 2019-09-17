;pro Tcevol ;Thermal instability criterion and Temperature evolution
func='T criteria'
iextrapol=-1
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
   ;print,'doask',doask
   if transform eq 'amr' then begin
     if n_elements(iextrapol) eq 0 then iextrapol=0
     asknum,'extrapolate to single grid (0/1/2)       ',iextrapol,doask
   endif
   help,noautodomain
   ; ***AMR: added argument to readplotpar
   readplotpar,ndim,cut,cut0,plotdim,nfunc,func,funcs,funcs1,funcs2,$
      plotmode,plotmodes,plottitle,plottitles,autorange,autoranges,doask,$
      transform,noautodomain,domainset
   print,'nx=',nx
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
                      eqpar,wnames,cut0,doask
                 endfor
               endif else begin
                 ; case of no amr
                 first= npict eq 0 and ifile eq 0
                 getlimits,first,nfunc,funcs,funcs1,funcs2,autoranges,$
                         fmax,fmin,x,w,xreg,wreg,usereg,physicss(ifile),$
                         eqpar,wnames,cut0,doask
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

   ;===== DO calculation

   ipict=0
   for ifile=0,nfile-1 do openfile,ifile+1,filenames(ifile),filetypes(ifile)

   error=0
   times=fltarr(npict)
   Tecc=fltarr(npict,2)
   zerol=fltarr(npict)
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
         case ndim of
         ;-------------- 1D ----------------------
         1: begin
          cp=13.05
          radius=0.
          for j=0,1 do begin
           if(radius eq 0.) then begin
            xp=cp
            av=getvalue(xp,funcs[j],xamr,wamr,ngrids,ndim,nxamr,nw,physics,eqpar,wnames)
           endif else begin
           if(cp-radius lt domainamr[0] or cp+radius gt domainamr[1]) then message,'The radius is too large!'
           npt=10
           dx=2.*radius/float(npt-1)
           av=0.
           for i=0,npt-1 do begin
             xp=cp-radius+i*dx
             av=av+getvalue(xp,funcs[j],xamr,wamr,ngrids,ndim,nxamr,nw,physics,eqpar,wnames)
           endfor
           av=av/float(npt)
           endelse
           Tecc[ipict,j]=av
          endfor
          times[ipict]=time
         end
         endcase
      endfor
      ipict=ipict+1
   endwhile

   Tmax=max(Tecc[*,0],maxid)
   timemax=times[maxid]
   print,'Tmax=',Tmax,'at time=',timemax
   times=times/3600.
   Tecc[*,1]=Tecc[*,1]*1.e7
   ytit='!16C!x (10!u-7!n erg cm!u-3!n s!u-1!n K!u-1!n)'
   yrangc=[-2,10]
   tzoom=[2.1,3.1]
   index=where(times ge tzoom[0] and times le tzoom[1])
   timesz=times[index]
   ccz=Tecc[index,1]
   Tez=Tecc[index,0]

   for i=0, n_elements(times)-1 do begin
        if Tecc[i,1] le 0. then begin
           cindx=i
           break
        endif 
   endfor
   print,'ctime=',times[cindx],' cc=',Tecc[cindx,1],' cT=',Tecc[cindx,0]
   xsize=18.0
   ysize=xsize
   charsiz=1.4
   chthick=1.5
   entry_dev=!d.name
   !P.FONT = 0
   set_plot,'ps'
   device,SET_FONT='Times-Roman',filename='Tcevol.eps',/encapsulated,xsize=xsize,ysize=ysize

   plot,times,Tecc[*,1],xtitle='Time (hr)',ystyle=4,xstyle=1,/nodata,background=255,color=0,$
        xthick=2,ythick=2,thick=2,charsize=charsiz,position=[0.12,0.55,0.88,0.9],xrange=[0,5],charthick=chthick
   axis,yaxis=1,ystyle=1,yrange=yrangc,/save,color=0,xthick=2,ythick=2,charsize=charsiz,charthick=chthick
   oplot,times,Tecc[*,1],linestyle=2,color=0,thick=5
   oplot,times,zerol,linestyle=1,color=0,thick=2
   axis,yaxis=0,ystyle=1,yrange=[1.e2,6.e6],/ylog,/save,ytitle='!16T!x (K)',color=0,xthick=2,ythick=2,charsize=charsiz,charthick=chthick
   oplot,times,Tecc[*,0],linestyle=0,color=0,thick=5
   xyouts,0.98,0.568,ytit,/normal,charsize=charsz,charthick=chthick,orientation=90

   x12=0.18
   dx=0.06
   y00=0.75
   dy=0.031
   mt=0.01
   plots,[x12,x12+dx],[y00,y00],/normal,linestyle=0,color=0,thick=5
   plots,[x12,x12+dx],[y00-1.*dy,y00-1.*dy],/normal,linestyle=2,color=0,thick=5
   xyouts,x12+dx+0.01,y00-mt,'!16T!x',color=0,/normal,charsize=charsiz,charthick=chthick
   xyouts,x12+dx+0.01,y00-1.*dy-mt,'!16C!x',color=0,/normal,charsize=charsiz,charthick=chthick

   plot,timesz,ccz,ystyle=4,xstyle=1,xrange=tzoom,/nodata,background=255,color=0,$
        xthick=2,ythick=2,thick=2,charsize=charsiz,position=[0.12,0.1,0.88,0.45],/noerase,charthick=chthick
   axis,yaxis=1,ystyle=1,yrange=yrangc,/save,color=0,xthick=2,ythick=2,charsize=charsiz,charthick=chthick
   oplot,timesz,ccz,linestyle=2,color=0,thick=5
   oplot,timesz,zerol,linestyle=1,color=0,thick=2
   axis,yaxis=0,ystyle=1,yrange=[1.e2,6.e6],/ylog,/save,ytitle='!16T!x (K)',color=0,xthick=2,ythick=2,charsize=charsiz,charthick=chthick
   oplot,timesz,Tez,linestyle=0,color=0,thick=5
   xyouts,0.435,0.012,'Time (hr)',/normal,charsize=charsz,charthick=chthick
   xyouts,0.98,0.118,ytit,/normal,charsize=charsz,charthick=chthick,orientation=90

   x12=0.18
   dx=0.06
   y00=0.36
   dy=0.031
   mt=0.01
   plots,[x12,x12+dx],[y00,y00],/normal,linestyle=0,color=0,thick=5
   plots,[x12,x12+dx],[y00-1.*dy,y00-1.*dy],/normal,linestyle=2,color=0,thick=5
   xyouts,x12+dx+0.01,y00-mt,'!16T!x',color=0,/normal,charsize=charsiz,charthick=chthick
   xyouts,x12+dx+0.01,y00-1.*dy-mt,'!16C!X',color=0,/normal,charsize=charsiz,charthick=chthick

   plots,[0.12,0.439],[0.45,0.55],/normal,linestyle=0,color=0,thick=2
   plots,[0.591,0.88],[0.55,0.45],/normal,linestyle=0,color=0,thick=2

   xar=0.623
   arrow,xar,0.32,xar,0.37,/normalized,thick=2
   arrow,xar,0.21,xar,0.165,/normalized,thick=2
   device,/close
   set_plot,entry_dev

   
end
