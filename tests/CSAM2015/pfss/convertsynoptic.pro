;pro convertsynoptic

file=dialog_pickfile()
data=readfits(file,hdr)
ss=size(data)
nlon=ss[1]
nlat=ss[2]
nloni=nlon
nlati=long(ss[2]*!dpi/2.)
print,'nlon,nlat',nlon,nlat
print,'nloni,nlati',nloni,nlati



;  set up grid
;  uniformly spaced grid
lat=linrange(nlati+1,-90.,90.)
lat=(lat(0:nlati-1)+lat(1:nlati))/2.
theta=(90-lat)*!dpi/180
  
lon=linrange(nloni+1,0.,360.)
lon=(lon(0:nloni-1)+lon(1:nloni))/2.
phi=lon*!dtor

dlatix=asin(linrange(nlat,nlat/2.-0.5,-nlat/2.+0.5)/nlat*2)*180/!dpi
dlonix=linrange(nlon+1,0.,360.) 
dlonix=(dlonix(0:nlon-1)+dlonix(1:nlon))/2.
dlatinterp=get_interpolation_index(reverse(dlatix),lat)
dloninterp=get_interpolation_index(dlonix,lon)
magout=interpolate(data,dloninterp,dlatinterp,/grid)

openw,lun,'hmisynopticmap2111.dat',/get_lun
writeu,lun,nloni
writeu,lun,nlati
writeu,lun,theta
writeu,lun,phi
writeu,lun,magout
free_lun,lun


contourlevel=200
levels1=findgen(contourlevel)/contourlevel*(max(magout)*0.2-min(magout)*0.2)+min(magout)*0.2
nlevel=n_elements(levels1)
ncolor=nlevel+1
bottom=0
c_color=indgen(ncolor)+bottom
charsizecb=1.3
charsizepl=2
symsz=0.4
thic=5.
loadct,0
window,0,xsize=1280,ysize=600
contour,magout,lon,lat,$
levels=levels1,c_colors=c_color,/FILL,XSTYLE=1,YSTYLE=1,charsize=charsizepl,/nodata
contour,magout,lon,lat,xtitle='longitude',ytitle='latitude',title='synoptic magnetogram',$
levels=levels1,c_colors=c_color,/FILL,XSTYLE=1,YSTYLE=1,charsize=charsizepl
write_png,'synopticmapraw.png',tvrd()
end
