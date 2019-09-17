nn=6
ss=fltarr(nn)
np=fltarr(nn)
np=[1,2,4,8,16,32]
ss=[1,2.12,4.08,11.24,13.8,12.7]
XYouts,0.1,0.1,/normal,'!17 '
lch=1.2
lth=2.0
plot_oo,np,ss,psym=4,xcharsize=lch,ycharsize=lch, $
     xthick=lth,ythick=lth,xrange=[0.9,64],yrange=[0.9,64],xstyle=1,ystyle=1,$
     xtitle='!17 Number of processors',ytitle='!17 Speedup', $
     title='!17 DD usage of AMRVAC'
oplot,np,ss,thick=lth
oplot,np,np,line=1,thick=lth
lch=0.7
XYouts,1.1,1,'52x104x104',charsize=lch
XYouts,2.2,2,'26x104x104',charsize=lch
XYouts,4.2,4,'52x52x52',charsize=lch
XYouts,6,12,'26x52x52',charsize=lch
XYouts,16,15,'52x26x26',charsize=lch
XYouts,32.2,13,'26x26x26',charsize=lch
XYouts,1,40,'!17 3D MHD jet simulation' 
end
