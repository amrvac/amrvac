"""2D plotting routines for vtu AMR data"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import griddata
import sys, time
from matplotlib.ticker import MaxNLocator
from scipy import ndimage
import copy
from scipy import interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
from copy import deepcopy

from amrvac_tools.vtkfiles import streamplot, read

if sys.platform == "win32":
# On Windows, the best timer is time.clock()
    default_timer = time.clock
else:
# On most other platforms the best timer is time.time()
    default_timer = time.time     

#=============================================================================
class polyplot():

    """Simple 2D plotter class using matplotlib, plots every cell in one patch"""
    
    def __init__(self,value,data,nlevels=256, grid=None, blocks=None, blockWidth = 8, blockHeight = 8, nlevel1=0,cmap='jet', min=None, max=None,
                 xrange=None, yrange=None, orientation='vertical', right=True, fixzoom=None, fixrange=None, fig=None, axis=None,
                 filenameout=None, clear=True,
                 edgecolor='k',nancolor='magenta',smooth=0,
                 swap=0,**kwargs):

        self.swap=swap
        self.nlevels=nlevels
        self.grid = grid
        self.blocks = blocks
        self.blockWidth = blockWidth
        self.blockHeight = blockHeight
        self.nlevel1 = nlevel1
        self.cmap=cmap
        self.orientation=orientation
        self.right=right # If True, colorbar is at the right, if False it is at the left
        self.cbarwidth=0.15
        self.cbarpad=0.70
        self.fixzoom=fixzoom
        self.fixrange=fixrange
        self.filenameout=filenameout
        self.clear=clear # If True, the figure is cleared when initializing amrplot, preventing us e.g. for drawing on a different subplot.

        self.fontsize=10
        self.fig_w=2.5
        self.fig_h=4
        self.dpi=300
        self.maxXticks=None
        self.maxYticks=None
        self.cbarticks=None
        self.edgecolor=edgecolor
        self.nancolor=nancolor
        self.smooth=smooth
            
        self.xrange=xrange
        self.yrange=yrange
        if xrange is None:
            self.xrange=[data.getBounds()[0],data.getBounds()[1]]
        if yrange is None:
            self.yrange=[data.getBounds()[2],data.getBounds()[3]]

        self.setValue(value,min=min,max=max)

        # If a figure and axis were not given, create new ones
        if fig is None:
            self.figure=plt.figure(figsize=(self.fig_w,self.fig_h),dpi=100)
        else:
            self.figure=fig

        if axis is None:
            self.ax = self.figure.gca()
        else:
            self.ax = axis

        self.show(var=value,data=data,min=min,max=max)
        if self.filenameout is None:
            self.figure.canvas.mpl_connect('button_press_event', self.onkey)


    def setValue(self,value,min=None,max=None):
        """Sets the min and max values of the data to saturate the display"""
        self.value=value
        self.min=min
        self.max=max
        if self.min==None:
            self.min=value.min()
        if self.max==None:
            self.max=value.max()

    def update(self,var=None,data=None,min=None,max=None,reset=None,fixrange=None,filenameout=None):
        """Prepare to re-draw the window, check if data was updated"""
        if var is not None:
            newdata = np.any(var!=self.value)
        else:
            newdata = False

        if newdata:
            self.value=var
            if fixrange is None:
                if min is None:
                    self.min=self.value.min()
                if max is None:
                    self.max=self.value.max()
        if data is not None:
            self.data=data
        if reset is not None:
            self.min=self.value.min()
            self.max=self.value.max()
        if min is not None:
            self.min = min
        if max is not None:
            self.max = max
        
        if filenameout is not None:
            self.filenameout=filenameout
            self.figure.set_size_inches( (self.fig_w,self.fig_h) )

        self.ax.set_rasterization_zorder(-9)

        # save the view for later:
        self.viewXrange=deepcopy( self.ax.xaxis.get_view_interval() )
        self.viewYrange=deepcopy( self.ax.yaxis.get_view_interval() )

        
    def info(self):
        """Print info to the console"""
        print('=======================================================')
        print('plotting range between %e and %e' % (self.min,self.max))
        if self.fixzoom is None:
            print('xrange = [%e,%e]     yrange = [%e,%e]' % (self.xrange[0],self.xrange[1],self.yrange[0],self.yrange[1]))
        else:
            print("""Fixing zoomlevel to
xrange = [%e,%e]     yrange = [%e,%e]""" % (
                self.viewXrange[0],self.viewXrange[1],self.viewYrange[0],self.viewYrange[1]))
        if self.nlevels<=1:
            print('Need more than one color-level, resetting nlevels')
            self.nlevels=256
        print('colormap = %s; nlevels=%d; orientation=%s' % (self.cmap,self.nlevels,self.orientation))
        if self.grid is not None:
            print('Also showing gridlines')
        if self.blocks is not None:
            print('Also showing blocks')
        print('=======================================================')
            
    def show(self,var=None,data=None,min=None,max=None,reset=None,fixrange=None,filenameout=None):
        """Draw the plotting-window"""
        t0 = default_timer()

        self.update(var=var,data=data,min=min,max=max,reset=reset,fixrange=fixrange,filenameout=filenameout)

        if self.clear:
            self.figure.clf()
            self.ax = self.figure.gca()     
        self.info()
     
        colormap = plt.cm.get_cmap(self.cmap, self.nlevels)

        valueRange=self.max-self.min
        if valueRange == 0.:
            valueRange = 1

        self.valueClip = np.clip(self.value,self.min,self.max)

        tdata0= default_timer()
        if (self.swap == 0):
            [myxlist,myylist]=self.data.getPointList()
        else:
            [myylist,myxlist]=self.data.getPointList()
        # List for regular cells with finite values (no infinitys or NaNs)
        self.xlist = [[] for i in range(self.nlevels)]
        self.ylist = [[] for i in range(self.nlevels)]
        ilevel = ((self.nlevels-1)*(self.valueClip-self.min)/valueRange).astype(int)
        # Special list to contain NaNs or other problematic cells
        self.xlistspecial = []
        self.ylistspecial = []
        for icell in range(self.data.ncells):
             if ilevel[icell] > -1 and ilevel[icell]<self.nlevels:
                  self.xlist[ilevel[icell]].extend(myxlist[5*icell:5*(icell+1)])
                  self.ylist[ilevel[icell]].extend(myylist[5*icell:5*(icell+1)])
             else:
                  self.xlistspecial.extend(myxlist[5*icell:5*(icell+1)])
                  self.ylistspecial.extend(myylist[5*icell:5*(icell+1)])
        tdata1=default_timer()

        # Fill cells with finite values
        for ilevel in range(self.nlevels):
            self.ax.fill(self.xlist[ilevel],self.ylist[ilevel],
                     facecolor=colormap(ilevel/(self.nlevels-1.)), closed=False, edgecolor='none',antialiased=False,zorder=-10)

        # Fill cells with special values
        if self.xlistspecial: # If the list of special cells is not empty
            print('WARNING: There are NaNs or Inftys, NaN color:',self.nancolor)
            self.ax.fill(self.xlistspecial,self.ylistspecial,
                     facecolor=self.nancolor, closed=False, edgecolor='none',antialiased=False,zorder=-10)

        if self.grid is not None:
            self.ax.fill(myxlist,myylist, 
                     facecolor='none', edgecolor=self.edgecolor,aa=True,linewidth=0.2,alpha=0.8)

        if self.blocks is not None:
            [myxBlockList,myyBlockList] = self.data.getPieces(self.blockWidth,self.blockHeight,self.nlevel1)
            self.ax.fill(myxBlockList,myyBlockList, 
                    facecolor='none', edgecolor=self.edgecolor,aa=True,linewidth=0.2,alpha=0.8)

        if self.orientation is not None:
            self.colorbar()
        
        if self.fixzoom is None:
            self.ax.set_xlim(self.xrange[0],self.xrange[1])
            self.ax.set_ylim(self.yrange[0],self.yrange[1])
        else:
            self.ax.set_xlim(self.viewXrange)
            self.ax.set_ylim(self.viewYrange)
        self.ax.set_aspect('equal')
# ticks:
        if self.maxXticks is not None:
            self.ax.xaxis.set_major_locator(MaxNLocator(self.maxXticks-1))
        if self.maxYticks is not None:
            self.ax.yaxis.set_major_locator(MaxNLocator(self.maxYticks-1))

        for tick in self.ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(self.fontsize)
        for tick in self.ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(self.fontsize)
        self.ax.xaxis.get_offset_text().set_fontsize(self.fontsize-2)
        self.ax.yaxis.get_offset_text().set_fontsize(self.fontsize-2)

        tend = default_timer()
        print('time for arranging the data= %f sec' % (tdata1-tdata0))
        print('Execution time = %f sec' % (tend-t0))
        print('=======================================================')
        if self.filenameout is None:
            plt.draw()

# Make the main axis active:
        plt.sca(self.ax)



    def colorbar(self,cax=None):
        """Draw the colorbar.
        """
        colormap = plt.cm.get_cmap(self.cmap, self.nlevels)
        m = plt.cm.ScalarMappable(cmap=colormap)
        m.set_array(self.valueClip)
        m.set_clim(vmin=self.min,vmax=self.max)
        
        if cax is None:
            divider = make_axes_locatable(self.ax)
            if self.orientation == 'horizontal':
                self.cax = divider.append_axes("bottom", self.cbarwidth, pad=self.cbarpad)
            else:
                if self.right:
                    self.cax = divider.append_axes("right", "4%", pad="6%")
                else:
                    self.cax = divider.append_axes("left", "4%", pad="6%")
        else:
            self.cax=cax
                    
        self.cbar=self.figure.colorbar(m, orientation=self.orientation,cax=self.cax)
                    
        self.cbar.solids.set_rasterized(True)

        if self.cbarticks is not None:
            self.cbar.locator=MaxNLocator(nbins=self.cbarticks-1)
            self.cbar.update_ticks()


        for tick in self.cbar.ax.get_xticklabels():
            tick.set_fontsize(self.fontsize)
        for tick in self.cbar.ax.get_yticklabels():
            tick.set_fontsize(self.fontsize)
        self.cbar.ax.xaxis.get_offset_text().set_fontsize(self.fontsize-2)
        self.cbar.ax.yaxis.get_offset_text().set_fontsize(self.fontsize-2)

    def save(self,filenameout=None):
        """Save the figure"""
        if filenameout is not None:
            self.filenameout=filenameout
        print('saving plot to file %s' % (self.filenameout))
        self.figure.set_size_inches( (self.fig_w,self.fig_h) )
        self.figure.savefig(self.filenameout, transparent=False,aa=True,dpi=self.dpi,interpolation='bicubic',bbox_inches='tight')
        self.filenameout=None

        
    def onkey(self,event):
        """
        Get data at mousepoint-location.  Press middle button at location to activate,
        press outside of plotting range to remove cross-hair.
        """
        try:
            if event.button !=2:
                return True
        except AttributeError:
                return True
        if event.xdata is None:
            try:
                self.selection.pop(0).remove()
                plt.show()
            except:
                return True
            return True
       
        icell=self.data.getIcellByPoint(event.xdata,event.ydata)
        try:
            self.selection.pop(0).remove()
        except AttributeError:
            print('first time defining selection, please wait for centerpoints...')
        except ValueError:
            pass
        except IndexError:
            pass
        self.selection=self.ax.plot(self.data.getCenterPoints()[icell,0],self.data.getCenterPoints()[icell,1],'m+', markersize=20,markeredgewidth=3)
        plt.show()
        self.data.showValues(icell)
        print('value=%e' %(self.value[icell]))
        if self.value[icell] != self.valueClip[icell]:
            print('exceeding plot range')
        print('=======================================================')


#=============================================================================
class rgplot(polyplot):
    """
    As polyplot, but use regridded data to display
    """
    
    def show(self,var=None,data=None,min=None,max=None,reset=None,fixrange=None,filenameout=None):
        """Draw the plotting-window"""

        t0 = default_timer()

        self.update(var=var,data=data,min=min,max=max,reset=reset,fixrange=fixrange,filenameout=filenameout)

        self.viewXrange=self.ax.xaxis.get_view_interval()
        self.viewYrange=self.ax.yaxis.get_view_interval()
        self.figure.clf()
        self.ax = self.figure.gca()        
        self.ax.set_rasterization_zorder(-9)
        self.info()

        if self.fixzoom is None:
            self.ax.set_xlim(self.xrange[0],self.xrange[1])
            self.ax.set_ylim(self.yrange[0],self.yrange[1])
        else:
            self.ax.set_xlim(self.viewXrange)
            self.ax.set_ylim(self.viewYrange)
        self.ax.set_aspect('equal')


        valueRange=self.max-self.min
        if valueRange == 0.:
            valueRange = 1

        self.valueClip = np.clip(self.value,self.min,self.max)

# Now regrid and imshow:
        tdata0= default_timer()
        xrange=[self.ax.get_xlim()[0],self.ax.get_xlim()[1]]
        yrange=[self.ax.get_ylim()[0],self.ax.get_ylim()[1]]

        nregrid = [self.dpi*self.fig_w,self.dpi*self.fig_h]

        if (self.swap == 0):
            CC=self.data.getCenterPoints()
        else:
            CC=self.data.getCenterPoints()[:,[1,0]]
        tmp0=np.complex(0,nregrid[0])
        tmp1=np.complex(0,nregrid[1])
        grid_x, grid_y = np.mgrid[xrange[0]:xrange[1]:tmp0, yrange[0]:yrange[1]:tmp1]

# regrid with linear interpolation and fall back to nearest neighbor at NaN, which should just be the borders.
        gridvar = griddata(CC, self.valueClip, (grid_x, grid_y), method='cubic',fill_value=float(np.NaN))
        isnan=np.isnan(gridvar)
        gridvar[isnan] = griddata(CC, self.value, (grid_x[isnan], grid_y[isnan]), method='nearest',fill_value=float(np.NaN))

# smooth the data slightly:
        gridvar = ndimage.gaussian_filter(gridvar, sigma=self.smooth)

        extent = xrange[0], xrange[1], yrange[0], yrange[1]
        gridvarClip = np.clip(gridvar,self.min,self.max)
        self.image=gridvarClip
        tdata1=default_timer()


        im2 = self.ax.imshow(gridvarClip.transpose(), cmap=plt.cm.get_cmap(self.cmap, self.nlevels), interpolation='bicubic',extent=extent,origin='lower',zorder=-10)
        norm = colors.Normalize(vmin=self.min, vmax=self.max)
        im2.set_norm(norm) 
#        Q=velovect(self.data.u1,self.data.u2,self.data,nvect=[20,20],scale=30,fig=self)

        if self.grid is not None:
            if (self.swap == 0):
                [myxlist,myylist]=self.data.getPointList()
            else:
                [myylist,myxlist]=self.data.getPointList()
            self.ax.fill(myxlist,myylist, 
                     facecolor='none', edgecolor=self.edgecolor,aa=True,linewidth=0.2,alpha=0.8)

        if self.blocks is not None:
            [myxBlockList,myyBlockList] = self.data.getPieces(self.blockWidth,self.blockHeight,self.nlevel1)
            self.ax.fill(myxBlockList,myyBlockList, 
                    facecolor='none', edgecolor=self.edgecolor,aa=True,linewidth=0.4,alpha=0.5)

        if self.orientation is not None:
            self.colorbar() 

# ticks:
        if self.maxXticks is not None:
            self.ax.xaxis.set_major_locator(MaxNLocator(self.maxXticks-1))
        if self.maxYticks is not None:
            self.ax.yaxis.set_major_locator(MaxNLocator(self.maxYticks-1))

        for tick in self.ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(self.fontsize)
        for tick in self.ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(self.fontsize)
        self.ax.xaxis.get_offset_text().set_fontsize(self.fontsize-2)
        self.ax.yaxis.get_offset_text().set_fontsize(self.fontsize-2)

        tend = default_timer()
        print('time for arranging the data= %f sec' % (tdata1-tdata0))
        print('Execution time = %f sec' % (tend-t0))
        print('=======================================================')
        if self.filenameout is None:
            plt.draw()

# Make the main axis active:
        plt.sca(self.ax)


#=============================================================================
class polyanim():
    """Animator class with polyplot"""
    
    def __init__(self,offsets=None,filename='data',type='pvtu',function=None,
                 nlevels=256, grid=None, cmap='jet',
                 min=None, max=None,
                 xrange=None, yrange=None, orientation='horizontal', filenameout='anim', **kwargs):
        
        self.nlevels=nlevels
        self.grid = grid
        self.cmap=cmap
        self.orientation=orientation
        self.offsets=offsets
        self.filename=filename
        self.type=type
        self.filenameout=filenameout
        if function is None:
            self.function = lambda x: x.rho
        else:
            self.function=function
            
        self.xrange=xrange
        self.yrange=yrange
        self.min=min
        self.max=max

    def setup(self):
        data=read.load(self.offsets[0],file=self.filename,type=self.type)

        if self.xrange is None:
            self.xrange=[data.getBounds()[0],data.getBounds()[1]]
        if self.yrange is None:
            self.yrange=[data.getBounds()[2],data.getBounds()[3]]
            
    def run(self):
        self.setup()
        
        offsets=np.arange(self.offsets[1]+1-self.offsets[0])+self.offsets[0]
        for offset in offsets:
            fo=''.join([self.filenameout,str(offset).zfill(4),'.png'])
            data=read.load(offset,file=self.filename,type=self.type)
            
#            exec ''.join(['var=',self.function])            
            var = self.function(data)
            self.polyplot = polyplot(var,data=data,
                                     nlevels=self.nlevels,grid=self.grid,cmap=self.cmap,
                                     orientation=self.orientation,
                                     xrange=self.xrange,yrange=self.yrange,
                                     min=self.min,max=self.max,
                                     filenameout=fo
                                     
                )
            self.polyplot.save(fo)
            plt.close()

            
class line():
    """returns array of cells on a given line"""
    
    def __init__(self,data,alice,bob):

        self.data=data
        self.alice=np.array(alice)
        self.bob=np.array(bob)
        self.icells=[]
        self.n=np.empty(2)
        self.n[0] = (self.bob[0]-self.alice[0])/np.sqrt((self.bob[0]-self.alice[0])**2+(self.bob[1]-self.alice[1])**2)
        self.n[1] = (self.bob[1]-self.alice[1])/np.sqrt((self.bob[0]-self.alice[0])**2+(self.bob[1]-self.alice[1])**2)
        self.epsilon = 2.e-1
        self.stepmax=10000

    def run(self):

        self.x=copy.deepcopy(self.alice)
        i=0
        while (self.n[0]*(self.bob[0]-self.x[0])+self.n[1]*(self.bob[1]-self.x[1]))>0 and i<self.stepmax:
            myIcell = self.data.getIcellByPoint(self.x[0],self.x[1])
            self.icells.append(myIcell)
            self.step(myIcell)
            i=i+1

    def step(self,icell):
        centerpoints=self.data.getCenterPoints()
        vert=self.data.getVert(icell)
        delta = np.min([np.abs(vert[0][0]-vert[1][0]),np.abs(vert[0][1]-vert[3][1])])/2.
        delta = delta * (1.-self.epsilon)
        myIcell=icell
        while (myIcell == icell and
               self.n[0]*(self.bob[0]-self.x[0])+self.n[1]*(self.bob[1]-self.x[1]))>0:
            myIcell=self.data.getIcellByPoint(self.x[0],self.x[1])
            self.x[0] = self.x[0] + delta * self.n[0]
            self.x[1] = self.x[1] + delta * self.n[1]

def plotoverline(var,data,alice,bob):
    fig=plt.figure()
    l=line(data,alice,bob)
    l.run()

    x=[]
    y=[]
    myvar=[]
    linecoords=data.getCenterPoints()[l.icells]
    for i in range(len(l.icells)):
        x.append(linecoords[i,0])
        y.append(linecoords[i,1])
    exec('myvar=%s[l.icells]' % var)
    x = np.array(x)
    y = np.array(y)
    s = np.sqrt(x**2+y**2)
    plt.plot(s,myvar,marker='+',linestyle='-')

    return l


#=============================================================================
def velovect(u1,u2,d,nvect=None,scalevar=None,scale=100,color='k',ax=None,alpha=1.):
    """Plots normalized velocity vectors"""

    minvel=1e-40

    if ax is None:
        ax=plt.gca()

    CC=d.getCenterPoints()
    n=np.sqrt(u1**2+u2**2)
    # remove zero velocity:
    m=n<minvel
    vr=np.ma.filled(np.ma.masked_array(u1/n,m),0.)
    vz=np.ma.filled(np.ma.masked_array(u2/n,m),0.)
    if scalevar is not None:
        vr = vr*scalevar
        vz = vz*scalevar
    if nvect is None:
        Q=ax.quiver(CC[:,0],CC[:,1],vr,vz,pivot='middle',width=1e-3,minlength=0.,scale=scale,
                    headwidth=6,alpha=alpha)
    else:
        # regrid the data:
        tmp0=np.complex(0,nvect[0])
        tmp1=np.complex(0,nvect[1])
        grid_r, grid_z = np.mgrid[ax.get_xlim()[0]:ax.get_xlim()[1]:tmp0, ax.get_ylim()[0]:ax.get_ylim()[1]:tmp1]
        grid_vr = griddata(CC, vr, (grid_r, grid_z), method='nearest')
        grid_vz = griddata(CC, vz, (grid_r, grid_z), method='nearest')
        Q=ax.quiver(grid_r,grid_z,grid_vr,grid_vz,pivot='middle',width=2e-3,minlength=minvel,scale=scale,
                    headwidth=10,headlength=10,color=color,edgecolor=color,rasterized=True,alpha=alpha)

    plt.draw()
    return Q     

#=============================================================================
def contour(var,d,levels=None,nmax=600,colors='k',linestyles='solid',ax=None,linewidths=1,smooth=1.,alpha=1.):
    if ax is None:
        ax=plt.gca()

    xrange=[ax.get_xlim()[0],ax.get_xlim()[1]]
    yrange=[ax.get_ylim()[0],ax.get_ylim()[1]]

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
    grid_x, grid_y = np.mgrid[xrange[0]:xrange[1]:tmp0, yrange[0]:yrange[1]:tmp1]
    gridvar = griddata(CC, var, (grid_x, grid_y), method='linear')
    isnan=np.isnan(gridvar)
    gridvar[isnan] = griddata(CC, var, (grid_x[isnan], grid_y[isnan]), method='nearest')

# smooth the data slightly:
    blurred_gridvar = ndimage.gaussian_filter(gridvar, sigma=smooth)


    if levels is None:
        cs = ax.contour(grid_x,grid_y,blurred_gridvar,16,alpha=alpha)
    else:
        cs = ax.contour(grid_x,grid_y,blurred_gridvar,levels=levels,colors=colors,linestyles=linestyles,linewidths=linewidths,alpha=alpha)

    plt.draw()
    return cs

#=============================================================================
def streamlines(u1,u2,d,x0=None,y0=None,nmax=600,density=1,fig=None,color='b',linewidth=1,arrowsize=1,alpha=1.,smooth=0):
    """plots streamlines from a vector field.  Use density=[densx,densy] to control how close streamlines are allowed to get."""
    if fig is None:
        ax=plt.gca()
    else:
        ax=fig.ax

    xrange=[ax.get_xlim()[0],ax.get_xlim()[1]]
    yrange=[ax.get_ylim()[0],ax.get_ylim()[1]]

    if nmax == 600 and fig is not None:
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
    print(nregrid)
    
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

    if (smooth != 0):
        u = ndimage.gaussian_filter(u, sigma=smooth)
        v = ndimage.gaussian_filter(v, sigma=smooth)

    if (x0 is not None and y0 is not None):
        for myx in zip(x0,y0):
            streamplot.streamplot(x, y, u.transpose(), v.transpose(), x_0=myx[0],
                       y_0=myx[1], density=density, linewidth=linewidth,
                       INTEGRATOR='RK4', color=color, arrowsize=arrowsize,alpha=alpha)
    else:
        streamplot.streamplot(x, y, u.transpose(), v.transpose(),
                   density=density, linewidth=linewidth,
                   INTEGRATOR='RK4', color=color, arrowsize=arrowsize,alpha=alpha)

#=============================================================================
def streamline(u1,u2,d,x1_0,x2_0,dl=1.,fig=None,nmax=600):

    if fig is None:
        ax=plt.gca()
    else:
        ax=fig.figure.gca()

    xrange=[ax.get_xlim()[0],ax.get_xlim()[1]]
    yrange=[ax.get_ylim()[0],ax.get_ylim()[1]]

    if nmax == 600 and fig is not None:
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
# Normalize:
    mag=np.sqrt(u**2+v**2)
    u=u/mag
    v=v/mag

    fu = interpolate.interp2d(grid_x, grid_y, u, kind='cubic')
    fv = interpolate.interp2d(grid_x, grid_y, v, kind='cubic')
    x1=[x1_0]
    x2=[x2_0]
    while(x1[-1] >= xrange[0] and x1[-1] <= xrange[1] and
    x2[-1] >= yrange[0] and x2[-1] <= yrange[1]):
        dt=dl
        x1.append(x1[-1]+fu(x1[-1],x2[-1])*dt)
        x2.append(x2[-1]+fv(x1[-1],x2[-1])*dt)
    return [np.array(x1),np.array(x2)]

#=============================================================================
def gridvar(var,d,xrange,yrange,nmax=600,smooth=0):

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
    grid_x, grid_y = np.mgrid[xrange[0]:xrange[1]:tmp0, yrange[0]:yrange[1]:tmp1]
    gridvar = griddata(CC, var, (grid_x, grid_y), method='linear')
    isnan=np.isnan(gridvar)
    gridvar[isnan] = griddata(CC, var, (grid_x[isnan], grid_y[isnan]), method='nearest')

# smooth the data slightly:
    blurred_gridvar = ndimage.gaussian_filter(gridvar, sigma=smooth)

    return blurred_gridvar,grid_x,grid_y
