import read as read
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import griddata
import vtk as v
import vtktonumpy as ah
import pickle
from IPython.parallel import Client
from matplotlib.ticker import MaxNLocator
from numba import autojit
from parmap import parmap as parmap
from numba import double
from numba import int32 as int32numba
from numba.decorators import jit

#from pylab import rcParams
#rcParams['xtick.direction'] = 'out'
#rcParams['ytick.direction'] = 'out'

#rc=Client()
#dview=rc[:]

path='movie/'
fileheader='aniso_av_'

def get(pointdata):

    print '=== Calculating local anisotropy ==='

    bcyl = get_bcyl(pointdata)

    m      = np.square(bcyl[0])+np.square(bcyl[2])+np.square(bcyl[1]) <= 1.e-14
 
    lfac=read.extract(pointdata,'lfac',attribute_mode='point')
    m = m | (lfac >= 9.)

    aniso = np.square(bcyl[0])+np.square(bcyl[2])/( np.square(bcyl[0])+np.square(bcyl[2])+np.square(bcyl[1]))

    m = m | (np.isnan(aniso))
    aniso = np.ma.masked_array(aniso,m)

    print '=== Done with local anisotropy ==='

    return [aniso,m]

def fieldstrength(pointdata):
    bx=read.extract(pointdata,'b1',attribute_mode='point')
    by=read.extract(pointdata,'b2',attribute_mode='point')
    bz=read.extract(pointdata,'b3',attribute_mode='point')

    lfac=read.extract(pointdata,'lfac',attribute_mode='point')

    b= np.sqrt(np.square(bx)+np.square(by)+np.square(bz))

    m = b <=1.e-7
    m = m | (lfac >= 9.)

    b = np.ma.masked_array(b,m)

    print '=== Done with field strength ==='

    return [b,m]

def get_bcyl(pointdata):

    bx=read.extract(pointdata,'b1',attribute_mode='point')
    by=read.extract(pointdata,'b2',attribute_mode='point')
    bz=read.extract(pointdata,'b3',attribute_mode='point')
    
    points=ah.vtk2array(pointdata.GetPoints().GetData())

    r      = np.sqrt(np.square(points[:,0])+np.square(points[:,1]))
    m      = r <= 0.

    cosphi = np.ma.array(points[:,0]/r,mask=m,fill_value=0.)
    sinphi = np.ma.array(points[:,1]/r,mask=m,fill_value=0.)

    br     = bx*cosphi + by*sinphi
    bphi   = by*cosphi - bx*sinphi

    return np.array([br,bphi,bz])

def load(offset,filenameout='data',type='pvtu'):
    if type == 'vtu':
        filename=''.join([filenameout,repr(offset).zfill(4),'.vtu'])
        datareader = v.vtkXMLUnstructuredGridReader()
    elif type == 'pvtu':
        filename=''.join([filenameout,repr(offset).zfill(4),'.pvtu'])
        datareader = v.vtkXMLPUnstructuredGridReader()
        
    print '=== loading and converting to pointdata, file:', filename,' ==='

    datareader.SetFileName(filename)
    datareader.Update()
    c2p = v.vtkCellDataToPointData()
    c2p.SetInput(datareader.GetOutput())
    pointdata=c2p.GetOutput()
    pointdata.Update()

    print '=== Done with loading pointdata ==='
    return pointdata

def plot(offset,nregrid,xrange=[-1.8e18,1.8e18],yrange=[-1.8e18,1.8e18],filenameout='data',type='pvtu'):

    nlevels=256

    path='movie'

    timescale = 1.05772480719860e-1
    time_start = 212.16
    time = offset*timescale + time_start

    pointdata=load(offset,filenameout=filenameout,type=type)
    aniso  = get(pointdata)
    points=ah.vtk2array(pointdata.GetPoints().GetData())

    x=points[:,0]
    y=points[:,1]
    z=points[:,2]

    # this is needed because griddata freaks out otherwise:
    smalldouble=1.e-8

    tmp=(xrange[1]-xrange[0])*smalldouble
    xi = np.linspace(xrange[0]-tmp,xrange[1]+tmp,nregrid[0])
    tmp=(yrange[1]-yrange[0])*smalldouble
    yi = np.linspace(yrange[0]-tmp,yrange[1]+tmp,nregrid[1])
    data = griddata((x,y), aniso[0], (xi[None,:],yi[:,None]), method='linear')
    mask = griddata((x,y), aniso[1], (xi[None,:],yi[:,None]), method='linear')
    data = np.ma.masked_array(data,mask)


    plt.clf()
    cs = plt.contourf(xi,yi,data,levels=np.arange(nlevels)/float(nlevels-1))
    cb = plt.colorbar()
    cb.set_ticks(np.arange(11)/10.)

    ax = plt.gca()
    ax.set_aspect('equal')

    plt.ylabel('z [cm]')
    plt.xlabel('y [cm]')
    plt.title(''.join(['Departure from tor. field,  t=',str(time)[0:5],' yrs']))

    plt.savefig(''.join([path,'/aniso_',str(offset).zfill(4),'.png']), format="png", transparent=False)

#    plt.show()
    plt.close()

def plot_range(beg,end,nregrid,xrange=[-1.8e18,1.8e18],yrange=[-1.8e18,1.8e18],filenameout='data',type='pvtu'):

    offsets=np.add(beg,range(end-beg+1))
    for i in offsets:
        plot(i,nregrid,filenameout=filenameout,type=type)
        print 'Offset', i, '; ', str(float(i-beg+1)/float(len(offsets))*100.)[0:4], '% done'

def average(aniso,points,nregrid,rrange=[0.,1.5e18],phirange=[0.,6.283185307179586],zrange=[-1.5e18,1.5e18]):
    
    x=points[:,0]
    y=points[:,1]
    z=points[:,2]

    # create r_i, phi_j and z_k arrays:
    ri   = np.linspace(rrange[0],rrange[1],nregrid[0])
    phij = np.linspace(phirange[0],phirange[1],nregrid[1])
    zk   = np.linspace(zrange[0],zrange[1],nregrid[2])

    # allocate xijk arrays:
    nelements = nregrid[0]*nregrid[1]*nregrid[2]
    xijk=np.empty(nelements)
    yijk=np.empty(nelements)
    zijk=np.empty(nelements)

# create single indexed points: 
    for i in range(nregrid[0]):
        for j in range(nregrid[1]):
            for k in range(nregrid[2]):
                xijk[i+k*nregrid[0]+j*nregrid[0]*nregrid[2]] = ri[i]*np.cos(phij[j])
                yijk[i+k*nregrid[0]+j*nregrid[0]*nregrid[2]] = ri[i]*np.sin(phij[j])
                zijk[i+k*nregrid[0]+j*nregrid[0]*nregrid[2]] = zk[k]

# griddata to points:
    data = griddata((x,y,z), aniso[0], (xijk,yijk,zijk), method='nearest')
    data = np.ma.masked_array(data,np.isnan(data))

# Now create the phi average:
#    aniso_av = avkernel(aniso[0],x,y,z,rij,phij,zk,np.array(nregrid,dtype=np.int16))
    aniso_av = np.zeros((nregrid[2],nregrid[0]))
    for i in range(nregrid[0]):
        for j in range(nregrid[1]):
            for k in range(nregrid[2]):
                aniso_av[k,i] = aniso_av[k,i] + data[i+k*nregrid[0]+j*nregrid[0]*nregrid[2]]
    aniso_av = aniso_av/np.double(nregrid[1])

    return np.ma.masked_array(aniso_av,np.isnan(aniso_av))



def averageNumba(aniso,points,nregrid,rrange=[0.,1.5e18],phirange=[0.,6.283185307179586],zrange=[-1.5e18,1.5e18]):

    print '=== Obtaining average of quantity ==='
    # create r_i, phi_j and z_k arrays:
    ri   = np.linspace(rrange[0],rrange[1],nregrid[0])
    phij = np.linspace(phirange[0],phirange[1],nregrid[1])
    zk   = np.linspace(zrange[0],zrange[1],nregrid[2])

    (xijk,yijk,zijk) = cylKernel(ri,phij,zk,np.array(nregrid,dtype=np.int32))
# griddata to points:
    data = griddata((points[:,0],points[:,1],points[:,2]), 
                    np.array(aniso[0].filled(float('nan')), dtype=np.float64),
                    (xijk,yijk,zijk), method='nearest')

    aniso_av = avKernel(data,np.array(nregrid,dtype=np.int32))

    print '=== Done with average of quantity ==='
    return np.ma.masked_array(aniso_av,np.isnan(aniso_av))

def stddevNumba(aniso,points,nregrid,rrange=[0.,1.5e18],phirange=[0.,6.283185307179586],zrange=[-1.5e18,1.5e18]):

    print '=== Obtaining standard deviation of anisotropy ==='
    # create r_i, phi_j and z_k arrays:
    ri   = np.linspace(rrange[0],rrange[1],nregrid[0])
    phij = np.linspace(phirange[0],phirange[1],nregrid[1])
    zk   = np.linspace(zrange[0],zrange[1],nregrid[2])

    (xijk,yijk,zijk) = cylKernel(ri,phij,zk,np.array(nregrid,dtype=np.int32))
# griddata to points:
    data = griddata((points[:,0],points[:,1],points[:,2]), 
                    np.array(aniso[0].filled(float('nan')), dtype=np.float64),
                    (xijk,yijk,zijk), method='nearest')

    aniso_std = stdKernel(data,np.array(nregrid,dtype=np.int32))

    print '=== Done with standard deviation of anisotropy ==='
    return np.ma.masked_array(aniso_std,np.isnan(aniso_std))

def correlationNumba(a,b,meanA,meanB,stddevA,stddevB,points,nregrid,rrange=[0.,1.5e18],phirange=[0.,6.283185307179586],zrange=[-1.5e18,1.5e18]):

    print '=== Obtaining correlation of anisotropy ==='
    # create r_i, phi_j and z_k arrays:
    ri   = np.linspace(rrange[0],rrange[1],nregrid[0])
    phij = np.linspace(phirange[0],phirange[1],nregrid[1])
    zk   = np.linspace(zrange[0],zrange[1],nregrid[2])

    (xijk,yijk,zijk) = cylKernel(ri,phij,zk,np.array(nregrid,dtype=np.int32))
# griddata to points:
    dataA = griddata((points[:,0],points[:,1],points[:,2]), 
                    np.array(a[0].filled(float('nan')), dtype=np.float64),
                    (xijk,yijk,zijk), method='nearest')
    dataB = griddata((points[:,0],points[:,1],points[:,2]), 
                    np.array(b[0].filled(float('nan')), dtype=np.float64),
                    (xijk,yijk,zijk), method='nearest')

    correlation = correlationKernel(dataA,dataB,meanA,meanB,stddevA,stddevB,np.array(nregrid,dtype=np.int32))

    print '=== Done with correlation of anisotropy ==='
    return np.ma.masked_array(correlation,np.isnan(correlation))

@autojit
def correlationKernel(dataA,dataB,meanA,meanB,stddevA,stddevB,nregrid):
# Now create the phi average:
    correlation = np.zeros((nregrid[2],nregrid[0]))
    for i in range(nregrid[0]):
        for j in range(nregrid[1]):
            for k in range(nregrid[2]):
                correlation[k,i] = correlation[k,i] + (
                    (dataA[i+k*nregrid[0]+j*nregrid[0]*nregrid[2]]-meanA[k,i]) *
                    (dataB[i+k*nregrid[0]+j*nregrid[0]*nregrid[2]]-meanB[k,i]) /
                    ((np.double(nregrid[1])-1.)*stddevA[k,i]*stddevB[k,i])
                    )

    return correlation

#@jit(argtypes=[double[:],double[:],double[:],int32numba[:]])
@autojit
def cylKernel(ri,phij,zk,nregrid):
    # allocate xijk arrays:
    nelements = nregrid[0]*nregrid[1]*nregrid[2]
    xijk=np.empty(nelements)
    yijk=np.empty(nelements)
    zijk=np.empty(nelements)

# create single indexed points: 
    for i in range(nregrid[0]):
        for j in range(nregrid[1]):
            for k in range(nregrid[2]):
                xijk[i+k*nregrid[0]+j*nregrid[0]*nregrid[2]] = ri[i]*np.cos(phij[j])
                yijk[i+k*nregrid[0]+j*nregrid[0]*nregrid[2]] = ri[i]*np.sin(phij[j])
                zijk[i+k*nregrid[0]+j*nregrid[0]*nregrid[2]] = zk[k]

    return (xijk,yijk,zijk)

#@jit(argtypes=[double[:],int32numba[:]])
@autojit
def avKernel(data,nregrid):
# Now create the phi average:
    aniso_av = np.zeros((nregrid[2],nregrid[0]))
    for i in range(nregrid[0]):
        for j in range(nregrid[1]):
            for k in range(nregrid[2]):
                aniso_av[k,i] = aniso_av[k,i] + data[i+k*nregrid[0]+j*nregrid[0]*nregrid[2]]
    aniso_av = aniso_av/np.double(nregrid[1])

    return aniso_av


#@jit(argtypes=[double[:],int32numba[:]])
@autojit
def stdKernel(data,nregrid):

    aniso_av = avKernel(data,nregrid)
# Now create the phi stddev:
    aniso_std = np.zeros((nregrid[2],nregrid[0]))
    for i in range(nregrid[0]):
        for j in range(nregrid[1]):
            for k in range(nregrid[2]):
                aniso_std[k,i] = aniso_std[k,i] + np.square(data[i+k*nregrid[0]+j*nregrid[0]*nregrid[2]]-aniso_av[k,i])
    aniso_std = np.sqrt(aniso_std/np.double(nregrid[1]))

    return aniso_std


def get_average(offset,nregrid,rrange=[0.,1.5e18],phirange=[0.,6.283185307179586],zrange=[-1.5e18,1.5e18],function=get,filenameout='data',type='pvtu'):
    # load the data:
    pointdata=load(offset,filenameout=filenameout,type=type)
    myaniso  = function(pointdata)
    points=ah.vtk2array(pointdata.GetPoints().GetData())

# Gathering data:
    data = averageNumba(myaniso,points,nregrid,rrange=rrange,phirange=phirange,zrange=zrange)
    
    return {'data':data,'nregrid':nregrid,'rrange':rrange,'phirange':phirange,'zrange':zrange,'offset':offset,'filenameout':filenameout,'type':type}

def get_correlation(offset,nregrid,rrange=[0.,1.5e18],phirange=[0.,6.283185307179586],zrange=[-1.5e18,1.5e18],functiona=get,functionb=fieldstrength,filenameout='data',type='pvtu'):
    # load the data:
    pointdata=load(offset,filenameout=filenameout,type=type)
    a  = functiona(pointdata)
    b  = functionb(pointdata)
    points=ah.vtk2array(pointdata.GetPoints().GetData())

# Gathering data:
    meanA = averageNumba(a,points,nregrid,rrange=rrange,phirange=phirange,zrange=zrange)
    meanB = averageNumba(b,points,nregrid,rrange=rrange,phirange=phirange,zrange=zrange)
    
    stddevA = stddevNumba(a,points,nregrid,rrange=rrange,phirange=phirange,zrange=zrange)
    stddevB = stddevNumba(b,points,nregrid,rrange=rrange,phirange=phirange,zrange=zrange)

    data = correlationNumba(
        a,b,meanA,meanB,stddevA,stddevB
,points,nregrid,rrange=rrange,phirange=phirange,zrange=zrange)

    return {'data':data,'nregrid':nregrid,'rrange':rrange,'phirange':phirange,'zrange':zrange,'offset':offset,'filenameout':filenameout,'type':type}

def plot_average(data,nlevels=256,plotrange=None,filename=None):

    offset=data.get('offset')
    rrange=data.get('rrange')
    zrange=data.get('zrange')
    averaged=data.get('data')
    nregrid=data.get('nregrid')

# Preliminaries:
    nlevels=256
    fontsize=10
    fig_w=3
    fig_h=4
    maxXticks=6
    maxYticks=9
    if filename==None:
        filename=''.join([path,fileheader,str(offset).zfill(4)])

    timescale = 0.105772480719860
    time_start = 212.16
    time = offset*timescale + time_start

    if plotrange==None:
        plotrange=np.array([0.,1.])
    plotrange=np.array(plotrange)
    
    averagedClip=averaged
    np.clip(averaged,plotrange[0],plotrange[1],out=averagedClip)

# Now plotting:
    xi = np.linspace(rrange[0],rrange[1],nregrid[0])
    yi = np.linspace(zrange[0],zrange[1],nregrid[2])
    figure=plt.figure(figsize=(fig_w,fig_h))
    ax = figure.gca()
    ax.set_rasterization_zorder(-9)

    cs = plt.contourf(xi,yi,averagedClip,levels=np.arange(nlevels)/float(nlevels-1)*(plotrange[1]-plotrange[0])+plotrange[0],aa=True,zorder=-10)
    cb = plt.colorbar()
    cb.set_ticks(np.arange(11)/10.*(plotrange[1]-plotrange[0])+plotrange[0])

    ax.set_aspect('equal')
    ax.xaxis.set_major_locator(MaxNLocator(maxXticks-1))
    ax.yaxis.set_major_locator(MaxNLocator(maxYticks-1))
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    ax.xaxis.get_offset_text().set_fontsize(fontsize-2)
    ax.yaxis.get_offset_text().set_fontsize(fontsize-2)

    for tick in cb.ax.get_yticklabels():
        tick.set_fontsize(fontsize)


    plt.ylabel(r'$z~\rm [cm]$',size=fontsize)
    plt.xlabel(r'$r~\rm [cm]$',size=fontsize)
    #    plt.title(''.join(['Anisotropy,  t=',str(time)[0:5],' yrs']),size=fontsize)
#    plt.title(''.join(['t=',str(time)[0:5],' yrs']),size=fontsize)

    plt.savefig(''.join([filename,'.pdf']), format="pdf", transparent=False,dpi=300,aa=True,bbox_inches='tight')
    plt.savefig(''.join([filename,'.png']), format="png", transparent=False,dpi=150,aa=True,bbox_inches='tight')

#    plt.show()
#    plt.close()

def plot_average_range(beg,end,nregrid,function=get,rrange=[0.,1.6e18],phirange=[0.,6.283185307179586],zrange=[-1.6e18,1.6e18],plotrange=np.array([0.,1.]),filenameout='data',type='pvtu'):


    offsets=np.add(beg,range(end-beg+1))
    for offset in offsets:
        average=get_average(offset,nregrid,rrange=rrange,phirange=phirange,zrange=zrange,function=function,filenameout=filenameout,type=type)
        plot_average(average,plotrange=plotrange)
#    averages=parmap(lambda offset, nregrid=nregrid,filenameout=filenameout,type=type:
#                            get_average(offset,nregrid,rrange=rrange,phirange=phirange,zrange=zrange,filenameout=filenameout,type=type),
#        offsets,np=2)
#    map(plot_average,
#        averages)
        #    for i in offsets:
        #        plot_average(i,nregrid,filenameout=filenameout,type=type)
        #        print 'Offset', i, '; ', str(float(i-beg+1)/float(len(offsets))*100.)[0:4], '% done'


def get1D(offset,nregrid,rrange=[0.,1.5e18],phirange=[0.,6.283185307179586],zrange=[-1.5e18,1.5e18],filenameout='data',type='pvtu'):

# load the data:
    pointdata=load(offset,filenameout=filenameout,type=type)
    aniso  = get(pointdata)
    points=ah.vtk2array(pointdata.GetPoints().GetData())

# Gathering data:
    datarz = average(aniso,points,nregrid,rrange=rrange,phirange=phirange,zrange=zrange)

    rsrange = [0.,np.sqrt(rrange[1]**2 + np.max(zrange)**2)]


# create r and z arrays:
    r   = np.linspace(rrange[0],rrange[1],nregrid[0])
    z   = np.linspace(zrange[0],zrange[1],nregrid[2])
    rr, zz = np.meshgrid(r,z)
    rr=np.reshape(rr,-1)
    zz=np.reshape(zz,-1)

# create rs_i and theta_j arrays:
    rsi    = np.linspace(rsrange[0],rsrange[1],nregrid[0])
    thetaj = np.linspace(phirange[0],phirange[1]/2.,nregrid[1])

# allocate xij arrays:
    nelements = nregrid[0]*nregrid[1]
    rij=np.empty(nelements)
    zij=np.empty(nelements)

# create single indexed points: 
    for i in range(nregrid[0]):
        for j in range(nregrid[1]):
            zij[j+i*nregrid[1]] = rsi[i]*np.cos(thetaj[j])
            rij[j+i*nregrid[1]] = rsi[i]*np.sin(thetaj[j])

# griddata to points:
    data = griddata((rr,zz), np.reshape(datarz,-1), (rij,zij), method='linear')
    data = np.ma.masked_array(data,np.isnan(data))
    
# Now create the theta average:
    aniso_av = np.zeros(nregrid[0])
    for i in range(nregrid[0]):
        index=np.arange(nregrid[1])+i*nregrid[1]
        aniso_av[i] = data[index].mean()

    aniso_av = np.ma.masked_array(aniso_av,np.isnan(aniso_av))

    return {'aniso_av':aniso_av,
            'nregrid':nregrid,
            'rrange':rrange,
            'phirange':phirange,
            'zrange':zrange,
            'rsph':rsi,
            'offset':offset,
            'filenameout':filenameout,
            'type':type
            }

def collect1D(offsetrange,nregrid,rrange=[0.,1.5e18],phirange=[0.,6.283185307179586],zrange=[-1.5e18,1.5e18],filenameout='data',type='pvtu',dumpfile='aniso1D.dat'):

    offsets=np.arange(offsetrange[1]-offsetrange[0]+1)+offsetrange[0]
    aniso=dview.map_sync(lambda offset,nregrid=nregrid,rrange=rrange,phirange=phirange,zrange=zrange,filenameout=filenameout,type=type:
                         get1D(offset,nregrid,rrange=rrange,phirange=phirange,zrange=zrange,filenameout=filenameout,type=type),
                         offsets)

    file=open(dumpfile,'wb') 
    pickle.dump(aniso,file)
    file.close()

    return aniso


def restore(dumpfile='aniso1D.dat'):

    file=open(dumpfile,'r')
    data = pickle.load(file)
    file.close()
    return data


def plot1D(data,filename=None):

    timescale = 0.105772480719860

# All for reordering the data:
    nelem=len(data)
    aniso=[]
    offsets=[]
    rs=[]
    for i in range(nelem):
        offsets.append(data[i].get('offset'))
        aniso.append(data[i].get('aniso_av'))

    rs=data[0].get('rsph')
    offsets = np.array(offsets)
    aniso   = np.array(aniso)
    rs      = np.array(rs)

    plt.figure()
    plt.subplot(1, 1, 1)
    plt.contourf(offsets*timescale,rs,aniso.transpose(),256,cmap=cm.bone)
    cb=plt.colorbar()
    plt.ylim((0,1.5e18))
    plt.ylabel('r [cm]')
    plt.xlabel('time [years]')

    if filename == None:
        plt.show()
    else:
        plt.savefig(filename, format="png", transparent=False)


def getline(data,theta):

    dr=(data['rrange'][1]-data['rrange'][0])/data['nregrid'][0]
    dz=(data['zrange'][1]-data['zrange'][0])/data['nregrid'][2]
    ds=min(dr,dz)
    send = np.sqrt(data['rrange'][1]**2+data['zrange'][1]**2)
    s = np.linspace(0.,send,send/ds)

    r = s * np.sin(theta*np.pi/180.)
    z = s * np.cos(theta*np.pi/180.)

    ri = np.linspace(data['rrange'][0],data['rrange'][1],data['nregrid'][0])
    zi = np.linspace(data['zrange'][0],data['zrange'][1],data['nregrid'][2])

    tmp0=np.complex(0,data['nregrid'][0])
    tmp1=np.complex(0,data['nregrid'][2])

    grid_x, grid_y = np.mgrid[data['rrange'][0]:data['rrange'][1]:tmp0, 
                              data['zrange'][0]:data['zrange'][1]:tmp1]

    online = griddata(np.array(zip(grid_x.flatten(),grid_y.flatten())),data['data'].data.transpose().flatten(),zip(r,z),method='linear')

    return {'theta':theta, 'r':s, 'data': online}
