set_plot,'ps'
device,/encapsulated
device,filename='amrvacDD.1.eps'
device,xsize=13,ysize=12,xoffset=0,yoffset=0
nn=6
ss=fltarr(nn)
np=fltarr(nn)
np=[1,2,4,8,16,32]
;ss=[1,2.12,4.08,11.24,13.8,12.7]
ssapo=[1,2.96,6.19,16.27,19.84,20.05]
ssman=[1,2.18,4.62,5.19,18.36,26.17]
tapo=[368.79,124.621,59.546,22.662,18.589,18.392]
tman=[325.119,149.007,70.352,62.624,17.712,12.425]
XYouts,0.1,0.1,/normal,'!17 '
lch=1.2
lth=2.0
plot_oo,np,ssapo,psym=4,xcharsize=lch,ycharsize=lch, $
;plot,np,tapo,psym=4,xcharsize=lch,ycharsize=lch, $
     xthick=lth,ythick=lth,xrange=[0.9,64],yrange=[0.9,64],xstyle=1,ystyle=1,$
;     xthick=lth,ythick=lth,xrange=[0.9,64],yrange=[0,380],xstyle=1,ystyle=1,$
     xtitle='!17 Number of processors',ytitle='!17 Speedup', $
;     xtitle='!17 Number of processors',ytitle='!17 Execution time', $
     title='!17 DD usage of AMRVAC'
;oplot,np,tapo,thick=lth
;oplot,np,tman,thick=lth
;oplot,np,tman,psym=2,thick=lth
;oplot,np,tapo(0)/np,line=1,thick=lth
;oplot,np,tman(0)/np,line=1,thick=lth
oplot,np,ssapo,thick=lth
oplot,np,ssman,thick=lth
oplot,np,ssman,psym=2,thick=lth
oplot,np,np,line=1,thick=lth
oplot,[32,32],[2,2],psym=2
Xyouts,34,2,'!17Manual'
oplot,[32,32],[4,4],psym=4
Xyouts,34,4,'!17Auto'
lch=0.7
XYouts,1.1,1,'52x104x104',charsize=lch
XYouts,2.2,2,'26x104x104',charsize=lch
XYouts,4.2,4,'52x52x52',charsize=lch
XYouts,6,9.5,'26x52x52',charsize=lch
XYouts,16,15,'52x26x26',charsize=lch
XYouts,32.2,21.5,'26x26x26',charsize=lch
XYouts,1,40,'!17 3D MHD jet simulation' 
device,/close
;end
