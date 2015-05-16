# -------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import griddata
import sys, time
from streamplot import streamplot
import read as read
import tools as tb
from matplotlib.ticker import MaxNLocator
from scipy import ndimage
import copy
import colorsys
from helper import scale_to_unit_interval as scale_to_unit_interval
import lic_internal

if sys.platform == "win32":
# On Windows, the best timer is time.clock()
    default_timer = time.clock
else:
# On most other platforms the best timer is time.time()
    default_timer = time.time     


def lic(u1,u2,d,nmax=600,fig=None,kernellen=31,color=None,saturation=0.8,lmin=0.05,lmax=0.95,vmin=None,vmax=None,ncolors=256):

    if fig==None:
        ax=plt.gca()
    else:
        ax=fig.figure.gca()

    xrange=[ax.get_xlim()[0],ax.get_xlim()[1]]
    yrange=[ax.get_ylim()[0],ax.get_ylim()[1]]

    if nmax == 600 and fig != None:
        nmax = fig.dpi * max([fig.fig_w,fig.fig_h])


    # Aspect ratio:
    r=(xrange[1]-xrange[0])/(yrange[1]-yrange[0])
    if r<1:
        ny=nmax
        nx=int(r*nmax)
    else:
        nx=nmax
        ny=int(nmax/r)
    nregrid = [nx,ny]
    
    CC=d.getCenterPoints()
    tmp0=np.complex(0,nregrid[0])
    tmp1=np.complex(0,nregrid[1])
    x=np.linspace(xrange[0],xrange[1],nregrid[0])
    y=np.linspace(yrange[0],yrange[1],nregrid[1])
    grid_x, grid_y = np.mgrid[xrange[0]:xrange[1]:tmp0, yrange[0]:yrange[1]:tmp1]

    u = griddata(CC, u1, (grid_x, grid_y), method='linear')
    v = griddata(CC, u2, (grid_x, grid_y), method='linear')
    uisnan=np.isnan(u)
    visnan=np.isnan(v)
    un = np.empty(np.shape(u))
    vn = np.empty(np.shape(v))
    un[uisnan] = griddata(CC, u1, (grid_x[uisnan], grid_y[uisnan]), method='nearest')
    vn[visnan] = griddata(CC, u2, (grid_x[visnan], grid_y[visnan]), method='nearest')
    u[uisnan]=un[uisnan]
    v[visnan]=vn[visnan]        
        

    texture = np.random.rand(nx,ny).astype(np.float32)
    vectors = np.zeros((ny,nx,2),dtype=np.float32)
    vectors[...,0] = u.transpose().astype(np.float32)
    vectors[...,1] = v.transpose().astype(np.float32)

    kernel = np.sin(np.arange(kernellen)*np.pi/kernellen)
    kernel = kernel.astype(np.float32)

    image = lic_internal.line_integral_convolution(vectors, texture, kernel)

    if color == None:
        plt.gray()
        plt.figure()
        plt.imshow(image,origin='lower',extent=(xrange[0],xrange[1],yrange[0],yrange[1]))
    else:
        col = griddata(CC, color, (grid_x, grid_y), method='linear')
        colisnan=np.isnan(col)
        col[colisnan] = griddata(CC, color, (grid_x[colisnan], grid_y[colisnan]), method='nearest')
        col = col.transpose()
        if vmin==None:
            vmin=col.min()
        if vmax==None:
            vmax=col.max()
        col = np.clip(col,vmin,vmax)
        color_scaled = (1.-scale_to_unit_interval(col))*0.7
        luminance_scaled = scale_to_unit_interval(image)*(lmax-lmin) + lmin
        image_rgb = np.empty((ny,nx,3))
        for i in range(image.shape[0]):
            for j in range(image.shape[1]):
                image_rgb[i,j,:] = colorsys.hls_to_rgb(color_scaled[i,j],luminance_scaled[i,j],saturation)
        # Create the colormap for the colorbar:
        colors_h = (1.-np.linspace(0,1,ncolors))*0.7
        colormap_rgb = np.empty((ncolors,3))
        for i in range(ncolors):
            colormap_rgb[i,:] = colorsys.hls_to_rgb(colors_h[i],0.5,saturation)
        cmap_instance = colors.ListedColormap(colormap_rgb,'custom_map',ncolors)

        figure=plt.figure()
        image=plt.imshow(image_rgb,origin='lower',extent=(xrange[0],xrange[1],yrange[0],yrange[1]))
        plt.set_cmap(cmap_instance)
        plt.clim(vmin,vmax)
        plt.colorbar()
