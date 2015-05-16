import read
import numpy as np
import matplotlib.pyplot as plt
import pickle
from amrplot import line
from helper import dash
from parmap import parmapLimitThreads as parmap
from scipy.signal import lfilter, lfilter_zi, filtfilt, butter

filenamermax='rmax.pdf'
filenamezmax='zmax.pdf'
filenamermin='rmin.pdf'
dpi = 300

onlysmooth = False

#    filenamermax='rmaxResolution.pdf'
#    filenamezmax='zmaxResolution.pdf'
#    filenamermin='rminResolution.pdf'

label=[r'$\sigma_0=0.01$', r'$\sigma_0=1$', r'$\sigma_0=3$', r'$\alpha=10^\circ, \sigma_0=3$','_None','_None','_None','_None']
linestyles=['-',':','--','-.','-',':','--','-.']
colors=['r','g','b','m','r','g','b','m']
linewidth=[0.75,0.75,0.75,0.75,2,2,2,2]
markers=['o','o','s','s','o','o','s','s']

def get(data):
    '''
Get the radii of the shock and nebula
    '''
    lfacpw=9.
    trcut=0.001
    ri=1.e18
    re=5.e18
    c=2.9979e10
    v=2.493e-2*c*ri/re

    
    CC=data.getCenterPoints()
    rmax = (CC[data.lfac>lfacpw][:,0].max() - CC[data.lfac>lfacpw][:,0].min())
    if CC[:,0].min()<0:
        rmax = rmax/2.
    zmax = (CC[data.lfac>lfacpw][:,1].max() - CC[data.lfac>lfacpw][:,1].min())
    if CC[:,1].min()<0:
        zmax = zmax/2.

# minimal extension of shock:
    rmin = (CC[data.lfac<lfacpw][:,0]**2 + CC[data.lfac<lfacpw][:,1]**2).min()
    rmin = np.sqrt(rmin)

# Volume averaged nebula radius:
    dV = (data.getSurface() * 2.*np.abs(np.pi*CC[:,0]) )[data.tr2>trcut]
    V  = dV.sum()
    if CC[:,0].min()<0:
        V = V/2.
    if CC[:,1].min()>0:
        V = V*2.
    rn = np.power(3.*V/(4.*np.pi),1./3.)

# rmax from Hydro expectation:
    rmaxhydro = 0.82*np.sqrt(2)*rn*np.sqrt(v/c)/np.sqrt((1-np.square(ri/rn)))

# along various radial angles:

# investigate theta = 0,5,10,20,30,45,60,90 for starters.
    thetas = np.array((0.,5.,10.,20.,30.,45.,60.,90.))
    smins  = np.zeros(len(thetas))
    R=2e18
    alice=[0,0]
    i=0
    for theta in thetas:
        bob    = [R*np.sin(theta*np.pi/180.),R*np.cos(theta*np.pi/180.)]
        x=[]
        y=[]
        myline = line(data,alice,bob)
        myline.run()
        
        linecoords=data.getCenterPoints()[myline.icells]
        for icell in range(len(myline.icells)):
            x.append(linecoords[icell,0])
            y.append(linecoords[icell,1])
        myvar=data.lfac[myline.icells]
        x = np.array(x)
        y = np.array(y)
        s = np.sqrt(x**2+y**2)
        smins[i] = (s[myvar<lfacpw]).min()
        i=i+1
        

    return {'rmax':rmax,'zmax': zmax,'rmin':rmin,'rn':rn,'rmaxhydro':rmaxhydro,'thetas': thetas, 'smin': smins}

def getFromOffset(offset,filenameout='data',type='pvtu'):
    data = read.load(offset,file=filenameout,type=type)
    myradii=get(data)
    myradii['offset'] = offset
    return myradii

def collect(beg,end,filenameout='data',type='pvtu',dumpfile='radii.dat'):

    offsets=np.arange(end-beg+1)+beg

    radii=[]
    radii=parmap(lambda offset, filenameout=filenameout,type=type:
              getFromOffset(offset,filenameout=filenameout,type=type),
                 offsets,np=4)
        #    for offset in offsets:
        #        data = read.load(offset,file=filenameout,type=type)
        #        myradii=get(data)
        #        myradii['offset'] = offset
        #        radii.append(myradii)

    file=open(dumpfile,'wb') 
    pickle.dump(radii,file)
    file.close()

    return radii

def restore(dumpfile='radii.dat'):

    file=open(dumpfile,'r')
    data = pickle.load(file)
    file.close()
    return data

def plot(radii):

    timescale = 1.05772480719860
    fontsize = 18

    # All for reordering the data:
    nelem=len(radii)
    rmax=[]
    rmaxhydro=[]
    zmax=[]
    rmin=[]
    offset=[]
    rn=[]
    smin=[]
    thetas=[]
    for i in range(nelem):
        offset.append(radii[i].get('offset'))
        rn.append(radii[i].get('rn'))
        rmax.append(radii[i].get('rmax'))
        zmax.append(radii[i].get('zmax'))
        rmin.append(radii[i].get('rmin'))
        rmaxhydro.append(radii[i].get('rmaxhydro'))
        smin.append(radii[i].get('smin'))
        thetas.append(radii[i].get('thetas'))
        
    offset=np.array(offset)
    rn=np.array(rn)
    rmax=np.array(rmax)
    zmax=np.array(zmax)
    rmin=np.array(rmin)
    rmaxhydro=np.array(rmaxhydro)
    smin=np.array(smin)
    thetas=np.array(thetas)

# Now ready to plot:
    plt.figure()
    plt.plot(offset*timescale,rmax/rn,'k+',linestyle='solid')
    plt.plot(offset*timescale,rmaxhydro/rn,'k--')
    plt.xlim((10,100))
    plt.ylim((0,0.3))

    plt.grid(True)
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

    plt.ylabel(r'$r_{max}$/$r_n$',size=fontsize)
    plt.xlabel(r'$t [years]$',size=fontsize)


    plt.figure()
    plt.plot(offset*timescale,zmax/rn,'k+',linestyle='solid')
    plt.plot(offset*timescale,0.23/0.82*rmaxhydro/rn,'k--')
    plt.ylim((0,0.3))

    plt.grid(True)
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

    plt.ylabel(r'$z_{max}$/$r_n$',size=fontsize)
    plt.xlabel('t [years]',size=fontsize)

# Plot the shock size on different angular directions: 
    f3=plt.figure()
    ax3 = f3.add_subplot(111)
    
    for iplot in range(len(smin[0,:])):
        ax3.plot(offset*timescale,smin[:,iplot],'+',markersize=10,
                 linestyle='dashed',dashes=dash(iplot),
                 label=''.join(['theta=',repr(thetas[0,iplot]).zfill(2),' deg']))

    ax3.grid(True)
    for tick in ax3.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax3.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

    ax3.set_ylabel(r'$R_{max}$(theta)',size=fontsize)
    ax3.set_xlabel('t [years]',size=fontsize)
    ax3.legend(handlelength=4,loc=2)

    plt.show()

def multiplotError(rad):

    timescale = 1.05772480719860
    ri=1.e18
    re=5.e18
    c=2.9979e10
    v=2.493e-2*c*ri/re
    alpha=6./5.
    t0=alpha*ri/v

    fontsize = 10
    params = {'legend.fontsize': fontsize}
    plt.rcParams.update(params)

    f1=plt.figure(figsize=(4,4*0.618))
    ax1 = f1.add_subplot(111)
    f2=plt.figure(figsize=(4,4*0.618))
    ax2 = f2.add_subplot(111)
    f3=plt.figure(figsize=(4,4*0.618))
    ax3 = f3.add_subplot(111)

    nplots=len(rad)


    for j in range(nplots):
        radii=rad[j]
    # All for reordering the data:
        nelem=len(radii)
        rmax=[]
        rmaxhydro=[]
        zmax=[]
        rmin=[]
        offset=[]
        rn=[]
        smin=[]
        thetas=[]
        for i in range(nelem):
            offset.append(radii[i].get('offset'))
            rn.append(radii[i].get('rn'))
            rmax.append(radii[i].get('rmax'))
            zmax.append(radii[i].get('zmax'))
            rmin.append(radii[i].get('rmin'))
            rmaxhydro.append(radii[i].get('rmaxhydro'))
            smin.append(radii[i].get('smin'))
            thetas.append(radii[i].get('thetas'))
            
        offset=np.array(offset)
        rn=np.array(rn)
        rmax=np.array(rmax)
        zmax=np.array(zmax)
        rmin=np.array(rmin)
        rmaxhydro=np.array(rmaxhydro)
        smin=np.array(smin)
        thetas=np.array(thetas)
        t=offset*timescale*(365.*24.*3600.)+t0
        rnalpha=ri*(t/t0)**alpha
        rmaxalpha=0.82*(rnalpha**(3./2.)*c**(-0.5)*(1.+alpha)**(0.5)*t**(-0.5)
                        *(1.-(t0/t)**(alpha+1.))**(-0.5))

# Filter the data:
        b, a = butter(3, 0.05, 'low')
        rmaxSmooth = filtfilt(b, a, rmax, padtype='even')
        zmaxSmooth = filtfilt(b, a, zmax, padtype='even')

    # Now ready to plot:
# Plot rmax:

        if onlysmooth == False :
            ax1.plot(t/(365.*24.*3600.),rmax/rn,''.join(['k',markers[j]]),linewidth=0.5,linestyle='-', markersize=2.5,label='_',markerfacecolor=colors[j],alpha=0.3)
        ax1.plot(t/(365.*24.*3600.),rmaxSmooth/rn,''.join([colors[j],'']),linewidth=linewidth[j],linestyle=linestyles[j], markersize=3.0,label=label[j])

#        if j==0:
#            ax1.plot(offset*timescale,rmaxhydro/rn,''.join(['k','']))
#            ax1.plot(offset*timescale,rmaxalpha/rnalpha,''.join(['k','']))
#            ax1.text(700, 0.075, r'H', fontsize=fontsize, ha='center', va='top')
            
#        ax1.grid(True)
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        
        ax1.set_xlim((250,1060))
        ax1.set_ylim((0.0,0.3))
                    
        ax1.set_ylabel(r'$r_{\rm max}/r_{\rm n}$',size=fontsize)
        ax1.set_xlabel(r'$\rm Nebula~age~[years]$',size=fontsize)

# Plot zmax:
        if onlysmooth == False :
            ax2.plot(offset*timescale,zmax/rn,''.join(['k',markers[j]]),linewidth=0.5,linestyle='-', markersize=2.5,label='_',markerfacecolor=colors[j],alpha=0.3)
        ax2.plot(offset*timescale,zmaxSmooth/rn,''.join([colors[j],'']),linewidth=linewidth[j],linestyle=linestyles[j], markersize=3.0,label=label[j])

        if j==0:
            ax2.plot(offset*timescale,0.23/0.82*rmaxalpha/rnalpha,''.join(['k','']))
            ax2.text(700, 0.039, r'H', fontsize=fontsize, ha='center', va='top')
            
        for tick in ax2.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        for tick in ax2.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        
        ax2.set_xlim((10,800))
        ax2.set_ylim((0,0.1))
                    
        ax2.set_ylabel(r'$z_{\rm max}/r_{\rm n}$',size=fontsize)
        ax2.set_xlabel(r'$t~\rm [years]$',size=fontsize)

# Plot zmin:
        ax3.plot(offset*timescale,rmin,''.join([colors[j],'']),linewidth=linewidth[j],linestyle=linestyles[j],label=label[j])
        ax3.set_ylim((1e+14,1e+16))
        ax3.set_yscale('log')
            
        ax3.grid(True)
        for tick in ax3.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        for tick in ax3.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        
        ax3.set_xlim((10,100))
                    
        ax3.set_ylabel(r'$z_{\rm min}~\rm [cm]$',size=fontsize)
        ax3.set_xlabel(r't~\rm [years]',size=fontsize)

        leg = ax1.legend(loc='best', fancybox=True)
        leg.get_frame().set_alpha(0.8)
#        ax2.legend()
        ax3.legend(loc=4)
        plt.draw()

    f1.savefig(filenamermax, transparent=False,aa=True,dpi=dpi,bbox_inches='tight')
    f2.savefig(filenamezmax, transparent=False,aa=True,dpi=dpi,bbox_inches='tight')
    f3.savefig(filenamermin, transparent=False,aa=True,dpi=dpi,bbox_inches='tight')



def multiplot(rad):

    timescale = 1.05772480719860
    fontsize = 10
    params = {'legend.fontsize': fontsize}
    plt.rcParams.update(params)

    f1=plt.figure(figsize=(4,4*0.618))
    ax1 = f1.add_subplot(111)
    f2=plt.figure(figsize=(4,4*0.618))
    ax2 = f2.add_subplot(111)
    f3=plt.figure(figsize=(4,4*0.618))
    ax3 = f3.add_subplot(111)

    nplots=len(rad)


    for j in range(nplots):
        radii=rad[j]
    # All for reordering the data:
        nelem=len(radii)
        rmax=[]
        rmaxhydro=[]
        zmax=[]
        rmin=[]
        offset=[]
        rn=[]
        smin=[]
        thetas=[]
        for i in range(nelem):
            offset.append(radii[i].get('offset'))
            rn.append(radii[i].get('rn'))
            rmax.append(radii[i].get('rmax'))
            zmax.append(radii[i].get('zmax'))
            rmin.append(radii[i].get('rmin'))
            rmaxhydro.append(radii[i].get('rmaxhydro'))
            smin.append(radii[i].get('smin'))
            thetas.append(radii[i].get('thetas'))
            
        offset=np.array(offset)
        rn=np.array(rn)
        rmax=np.array(rmax)
        zmax=np.array(zmax)
        rmin=np.array(rmin)
        rmaxhydro=np.array(rmaxhydro)
        smin=np.array(smin)
        thetas=np.array(thetas)

    # Now ready to plot:
# Plot rmax:
        ax1.plot(offset*timescale,rmax/rn,''.join([colors[j],'']),linewidth=linewidth[j],linestyle=linestyles[j],label=label[j])
        if j==1:
            ax1.plot(offset*timescale,rmaxhydro/rn,''.join(['k','']))
            ax1.text(95, 0.13, r'H', fontsize=fontsize, ha='center', va='top')
            
#        ax1.grid(True)
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        
        ax1.set_xlim((10,100))
        ax1.set_ylim((0.0,0.3))
                    
        ax1.set_ylabel(r'$r_{\rm max}/r_{\rm n}$',size=fontsize)
        ax1.set_xlabel(r'$t~\rm [years]$',size=fontsize)

# Plot zmax:
        ax2.plot(offset*timescale,zmax/rn,''.join([colors[j],'']),linewidth=linewidth[j],linestyle=linestyles[j],label=label[j])
        if j==1:
            ax2.plot(offset*timescale,0.23/0.82*rmaxhydro/rn,''.join(['k','']))
            ax2.text(95, 0.036, r'H', fontsize=fontsize, ha='center', va='top')
            
        for tick in ax2.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        for tick in ax2.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        
        ax2.set_xlim((10,100))
        ax2.set_ylim((0,0.1))
                    
        ax2.set_ylabel(r'$z_{\rm max}/r_{\rm n}$',size=fontsize)
        ax2.set_xlabel(r'$t~\rm [years]$',size=fontsize)

# Plot zmin:
        ax3.plot(offset*timescale,rmin,''.join([colors[j],'']),linewidth=linewidth[j],linestyle=linestyles[j],label=label[j])
        ax3.set_ylim((1e+14,1e+16))
        ax3.set_yscale('log')
            
        ax3.grid(True)
        for tick in ax3.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        for tick in ax3.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        
        ax3.set_xlim((10,100))
                    
        ax3.set_ylabel(r'$z_{\rm min}~\rm [cm]$',size=fontsize)
        ax3.set_xlabel(r't~\rm [years]',size=fontsize)

        leg = ax1.legend(loc='best', fancybox=True)
        leg.get_frame().set_alpha(0.8)
#        ax2.legend()
        ax3.legend(loc=4)
        plt.draw()

    f1.savefig(filenamermax, transparent=False,aa=True,dpi=dpi,bbox_inches='tight')
    f2.savefig(filenamezmax, transparent=False,aa=True,dpi=dpi,bbox_inches='tight')
    f3.savefig(filenamermin, transparent=False,aa=True,dpi=dpi,bbox_inches='tight')




def aspectratios(rad):

    timescale = 1.05772480719860
    fontsize = 18
    f1=plt.figure()
    ax1 = f1.add_subplot(111)

    nplots=len(rad)
    color=['r','g','b']
    label=['sigma=0.01','sigma=1','sigma=3']
    for j in range(nplots):
        radii=rad[j]
    # All for reordering the data:
        nelem=len(radii)
        rmax=[]
        rmaxhydro=[]
        zmax=[]
        rmin=[]
        offset=[]
        rn=[]
        for i in range(nelem):
            offset.append(radii[i].get('offset'))
            rn.append(radii[i].get('rn'))
            rmax.append(radii[i].get('rmax'))
            zmax.append(radii[i].get('zmax'))
            rmin.append(radii[i].get('rmin'))
            rmaxhydro.append(radii[i].get('rmaxhydro'))
            
        offset=np.array(offset)
        rn=np.array(rn)
        rmax=np.array(rmax)
        zmax=np.array(zmax)
        rmin=np.array(rmin)
        rmaxhydro=np.array(rmaxhydro)

    # Now ready to plot:
# Plot aspect ratio:
        ax1.plot(offset*timescale,rmax/zmax,''.join([color[j],'+']),linestyle='solid',label=label[j])
        ax1.axhline(y=0.82/0.23,linestyle='--',color='k')
            
        ax1.grid(True)
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        
        ax1.set_xlim((10,100))
        ax1.set_ylim((2,14))

        ax1.set_ylabel(r'$r_{\rm max}$/$z_{\rm max}$',size=fontsize)
        ax1.set_xlabel('t~\rm [years]',size=fontsize)

        ax1.legend()
