import numpy as np
import read,amrplot,synchroSun,emissSun
import matplotlib.colors as colors

# usage: 
# $ ipython
# In [1]: import aia
# In [2]: aia.euv('promfluxroperelaxb',80,80,171,0,0,0,24,18,12,100,100,100)
def euv(filenameout,nstart,nend,ion,theta,phiy,phi,xra,yra,los,xn,yn,losn,Teunit=1.e6,nunit=1.e9,Lunit=1.e9,\
    function=None,vmin=None,vmax=None,xlabel='x',ylabel='y',xcenter=0,ycenter=0):
    '''Synthetic view along z-axis of ion wavelength from snapshots nstart to nend 
    after clock-wisely rotating data cube theta degree around x-axis, phiy degree 
    around y-axis, and then phi around z-axis. Note that x, y, and z axes are 
    corotating with the data cube. The output image has a size of xra by yra with
    the resolution of xn by yn, the line-of-sight depth is los and its resolution
    is losn. Function can be 'log' or 'sqrt', vmin/vmax set the lower/upper value
    limit of the image.'''
    if nend < nstart:
        nend=nstart
    nend=nend+1
    for i in range(nstart,nend):
        # import data
        data=synchroSun.get_pointdata(i,filenameout=filenameout,type='vtu')
        # get emissivity
        ed=emissSun.emiss(data,ion,Teunit,nunit,Lunit)
        j=emissSun.viewangle(ed,theta,phiy,phi)
        map=synchroSun.shotgun(j,[xn,yn,losn],[xra,yra,los])
        # create color table for the waveband
        cdict=emissSun.aiact(ion)
        my_cmap = colors.LinearSegmentedColormap('my_colormap',cdict,256)
        if nstart == nend-1:
        # show an interactive window
            synchroSun.show_map(map,function=function,vmax=vmax,cmap=my_cmap,xlabel=xlabel,ylabel=ylabel,xcenter=xcenter,ycenter=ycenter)
        else:
        # write png and pdf image files
            filen='synAIA'+str(ion)+'n'+str(i)
            synchroSun.show_map(map,filename=filen,function=function,vmax=vmax,cmap=my_cmap,xlabel=xlabel,ylabel=ylabel,xcenter=xcenter,ycenter=ycenter)
