;pro converthmi

;fn1 = 'hmi.sharp_cea_720s.3894.20140328_140000_TAI.Bp.fits'    ;B_phi = Bx
;read_sdo,fn1,index1,data1
;index2map,index1,data1,mapbx
;fn2 = 'hmi.sharp_cea_720s.3894.20140328_140000_TAI.Br.fits'
fn2=dialog_pickfile()
read_sdo,fn2,index2,data2
index2map,index2,data2,mapbz
;fn3 = 'hmi.sharp_cea_720s.3894.20140328_140000_TAI.Bt.fits'      ;B_theta equals -By
;read_sdo,fn3,index3,data3
;index2map,index3,data3,mapby
;mapby.data = -mapby.data

intimei = utc2int(mapbz.time)
mjdi  = intimei.mjd + fix((intimei.time + 120000.0)/(24.0*3600.0*1000.0))
timei = (intimei.time + 120000.0) mod (24.0*3600.0*1000.0)
intimei.mjd = mjdi
intimei.time= timei
utctime = anytim2utc(intimei,/vms)
mapbz.time = utctime
;mapbx.time = utctime
;mapby.time = utctime

carlon = tim2carr(utctime)   ;The Carrington longitude of the central meridian of the sun
lon = mapbz.xc - carlon      ; mapbz.xc records the Carrington longitude of the center of the field of view
lat = mapbz.yc
xcyc = lonlat2xy([lon[0],lat],utctime)
;mapbx.xc = xcyc[0]
;mapbx.yc = xcyc[1]
;mapby.xc = xcyc[0]
;mapby.yc = xcyc[1]
mapbz.xc = xcyc[0]
mapbz.yc = xcyc[1]
rsun = mapbz.rsun    ;in arcsec
dx = 2.0*!pi*rsun*0.03/360.    ;in arcsec. The CEA vector magnetic has a spatial resolution of 0.03 degree in longitude and latitude.
;mapbx.dx = dx
;mapbx.dy = dx
;mapby.dx = dx
;mapby.dy = dx
mapbz.dx = dx
mapbz.dy = dx

loadct,0
wdef,1,800,800
!p.background = 255
plot_map,mapbz,dmax=2000,dmin=-2000,color=0,charsize=1.8
write_png,'bz.png',tvrd()

;prepare input data file for lfff extrapolation in AMRVAC
arcsec2cm=7.25e7 ; in cm
dx=dx*arcsec2cm ; in cm
dy=dx
xc=mapbz.xc*arcsec2cm
yc=mapbz.yc*arcsec2cm
Bz=mapbz.data
sizebz=size(Bz)
nx1=sizebz[1]
nx2=sizebz[2]
filename='hmi20140328_140000.dat'
openw,lun,filename,/get_lun
writeu,lun,nx1
writeu,lun,nx2
writeu,lun,double(xc)
writeu,lun,double(yc)
writeu,lun,double(dx)
writeu,lun,double(dy)
writeu,lun,double(Bz)
free_lun,lun
print,'Bz range (Gauss):', min(Bz),max(Bz)
;wdef,2,800,800
;tvscl,Bz
print,'nx1,nx2',nx1,nx2
print,'xc,yc (cm)',xc,yc
print,'dx,dy (cm)',dx,dy
x1=xc-nx1*dx/2
x2=xc+nx1*dx/2
y1=yc-nx1*dx/2
y2=yc+nx1*dx/2
; output 
print,'x range (10 Mm):',' xprobmin1:',x1*1.e-9,'  xprobmax1:',x2*1.e-9
print,'y range (10 Mm):',' xprobmin1:',y1*1.e-9,'  xprobmax2:',y2*1.e-9
;sub_map,mapbz,smapbz,xrange=[0,779],yrange=[0,571],/pixel    ;,xrange=[500,700.0],yrange=[-300,-160.0]
;print,'Flux balance coefficient:',total(smapbz.data)/total(abs(smapbz.data))
;sub_map,mapbx,smapbx,ref_map=smapbz                                 
;sub_map,mapby,smapby,ref_map=smapbz
;wdef,1,800,800
;plot_map,smapbz,dmax=2000,dmin=-2000,color=0,charsize=1.8
;write_png,'sbz.png',tvrd()
;wdef,2,800,800
;plot_map,smapbz,dmax=2000,dmin=-2000,color=0,charsize=1.8
;plot_vmap,/over,smapbx,smapby,mapbz=smapbz,limit=150,scale=0.005,iskip=10,jskip=10,$
;          v_color=255,axis_color=0,/Nolabels,v_thick=2.0   ;,/No_arrow_head  ;,/Noaxis
;write_png,'sbxyz.png',tvrd()
;save,smapbx,smapby,smapbz,filename='bxyz_submap.sav'

end
