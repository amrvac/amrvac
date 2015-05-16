set_plot,'ps'
device,/encapsulated
device,filename='amr.eps'
device,xsize=13,ysize=12,xoffset=0,yoffset=0

nn=5
ss=fltarr(nn)
np=fltarr(nn)
ttime=fltarr(nn)
bctime=fltarr(nn)
amrtime=fltarr(nn)
rtime=fltarr(nn)
np=[1,2,4,8,16]
;ss=[1,1.715,2.910,3.526,4.05]
;ttime=[2852.883,1663.038,980.268,809.104,704.335]
;bctime=[846.287,510.514,314.914,320.026,276.548]
;amrtime=[228.853,165.627,129.04,141.114,119.799]
ttime=[2852.883,1663.038,980.268,728.603,581.299]
bctime=[846.287,510.514,314.914,262.231,180.716]
amrtime=[228.853,165.627,129.04,127.044,114.529]
rtime=ttime-bctime-amrtime
ss=ttime(0)/ttime
print,'speedup'
print,ss
print,'advance speedup'
print,rtime(0)/rtime
XYouts,0.1,0.1,/normal,'!17 '
lch=1.2
lth=2.0
plot_oo,np,ss,psym=4,xcharsize=lch,ycharsize=lch, $
     xthick=lth,ythick=lth,xrange=[0.9,20],yrange=[0.9,20],xstyle=1,ystyle=1,$
     xtitle='!17 Number of processors',ytitle='!17 Speedup', $
     title='!17 Full AMR scaling'
oplot,np,ss,thick=lth
oplot,np,np,line=1,thick=lth
oplot,np,rtime(0)/rtime,line=2,thick=lth
oplot,np,rtime(0)/rtime,psym=2
oplot,np,bctime(0)/bctime,line=3,thick=lth
oplot,np,bctime(0)/bctime,psym=6
oplot,np,amrtime(0)/amrtime,line=4,thick=lth
oplot,np,amrtime(0)/amrtime,psym=5
XYouts,1,15,'!17 2D HD Rayleigh-Taylor'
oplot,[1,1.5],[10,10]
Xyouts,1.6,10,'!17total'
oplot,[1,1.5],[8,8],line=2
Xyouts,1.6,8,'!17advance'
oplot,[1,1.5],[6,6],line=3
Xyouts,1.6,6,'!17BC'
oplot,[1,1.5],[5,5],line=4
Xyouts,1.6,5,'!17AMR'
;end
device,/close
