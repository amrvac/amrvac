set_plot,'ps'
device,/encapsulated
device,filename='amrvacDD.2.eps'
device,xsize=13,ysize=12,xoffset=0,yoffset=0

nn=5
ss=fltarr(nn)
np=fltarr(nn)
np=[1,2,4,8,16]
ss=[1,2.12,3.85,4.76,6.34]
ssapo=[1,72.066/72.254,72.066/71.876,72.066/72.921,72.066/92.370]
ss24=[1,49.785/25.288,49.785/13.987,49.785/8.898,49.785/5.599]
ss24apo=[1,48.529/27.625,48.529/17.531,48.529/14.511,48.529/14.424]
XYouts,0.1,0.1,/normal,'!17 '
lch=1.2
lth=2.0
plot_oo,np,ssapo,psym=2,xcharsize=lch,ycharsize=lch, $
     xthick=lth,ythick=lth,xrange=[0.9,20],yrange=[0.5,20],xstyle=1,ystyle=1,$
     xtitle='!17 Number of processors',ytitle='!17 Speedup', $
     title='!17 DD usage of AMRVAC'
oplot,np,ssapo,thick=lth
oplot,np,ss,thick=lth
oplot,np,ss,psym=2,thick=lth
oplot,np,ss24,thick=lth
oplot,np,ss24,psym=4,thick=lth
oplot,np,ss24apo,thick=lth
oplot,np,ss24apo,psym=4,thick=lth
oplot,np,np,line=1,thick=lth
lch=0.7
oplot,[1,1],[15,15],psym=2
XYouts,1,15,'!17 3D MHD jet 512 blocks (10!u3!n)' 
oplot,[1,1],[12,12],psym=4
XYouts,1,12,'!17 3D MHD jet 64 blocks (20!u3!n)' 
XYouts,14,0.65,'!17Auto'
XYouts,14,2.7,'!17Auto'
XYouts,12,7,'!17Manual'
XYouts,14,14,'!17Ideal'
;end
device,/close
