import read as read
import numpy as np
import scipy as scipy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from scipy.interpolate import griddata
from congrid import congrid
import vtk as v
#import vtktonumpy as ah
#from tvtk import array_handler as ah
import numpy_support as ah
import img_scale

def get_pointdata(offset,filenameout='data',type='pvtu'):

    if type == 'vtu':
        filename=''.join([filenameout,repr(offset).zfill(4),'.vtu'])
        datareader = v.vtkXMLUnstructuredGridReader()
    elif type == 'pvtu':
        filename=''.join([filenameout,repr(offset).zfill(4),'.pvtu'])
        datareader = v.vtkXMLPUnstructuredGridReader()
        
    print '=== Reading ',filename,' ==='
        
    datareader.SetFileName(filename)
    datareader.Update()
    data = datareader.GetOutput()
    c2p = v.vtkCellDataToPointData()
    c2p.SetInput(data)
    pointdata=c2p.GetOutput()
    pointdata.Update()
 
    vtk_points=data.GetPoints().GetData()
    points=ah.vtk2array(vtk_points)
 
    print '=== Done with reading data! ==='
    return {'pointdata': pointdata, 'points': points}

def rot2D(offset,nphi,filenameout='data',type='pvtu'):
    """
    creates a mock vtu datastructure with rotated values.
    No connectivity is set, vtu is just used as a container for pointdata     since
    the following routines expect to get a vtu unstructured datatype.
    """

    data=get_pointdata(offset,filenameout=filenameout,type='pvtu')
    points2D=data.get('points')
    data3D = v.vtkUnstructuredGrid()

    nelem2D = len(points2D)
    nelem3D = nelem2D*nphi

    print 'preparing ',nelem3D, 'points'

    points3D = np.empty((nelem3D,3))

    phis=np.linspace(0.,2.*np.pi*float(nphi)/float(nphi+1),nphi)

    i3D = 0
    for i2D in range(nelem2D):
# pick next 2D point:
        point = points2D[i2D]
        r = point[0]
        z = point[1]
        for phi in phis:
            sinphi = np.sin(phi)
            cosphi = np.cos(phi)
            points3D[i3D,0] = r*cosphi
            points3D[i3D,1] = r*sinphi
            points3D[i3D,2] = z
            i3D=i3D+1

    vtkpoints=ah.array2vtkPoints(points3D)
    data3D.SetPoints(vtkpoints)

# rotate variables:
    ivar = 0
    while ivar < data.get('pointdata').GetPointData().GetNumberOfArrays():
        var = data.get('pointdata').GetPointData().GetArrayName(ivar)
        print 'rotating variable:',ivar,var
        i3D=0
# Treat vectors:
        if var == 'u1' or var == 'b1' or var == 'v1':
            array2D_vec = np.empty((nelem2D,3))
            array3D_vec = np.empty((nelem3D,3))
            bx = np.empty(nelem3D)
            by = np.empty(nelem3D)
            bz = np.empty(nelem3D)
            array2D_vec[:,0] = read.extract(data.get('pointdata'),data.get('pointdata').GetPointData().GetArrayName(ivar),attribute_mode='point')
            array2D_vec[:,1] = read.extract(data.get('pointdata'),data.get('pointdata').GetPointData().GetArrayName(ivar+1),attribute_mode='point')
            array2D_vec[:,2] = read.extract(data.get('pointdata'),data.get('pointdata').GetPointData().GetArrayName(ivar+2),attribute_mode='point')
            for i2D in range(nelem2D):
                for phi in phis:
# Should have three components (makes no sense otherwise):
# Input in cylindrical coordinates, convention:
# b1==br; b2==bz; b3==bphi  and similar for velocities.
                    array3D_vec[i3D,0] = array2D_vec[i2D,0]
                    array3D_vec[i3D,1] = array2D_vec[i2D,1]
                    array3D_vec[i3D,2] = array2D_vec[i2D,2]
# Now create the rotated vectors:
                    sinphi = np.sin(phi)
                    cosphi = np.cos(phi)
                    bz[i3D] = array3D_vec[i3D,1]  
                    bx[i3D] = array3D_vec[i3D,0] * cosphi - array3D_vec[i3D,2] * sinphi
                    by[i3D] = array3D_vec[i3D,2] * cosphi + array3D_vec[i3D,0] * sinphi
                    i3D=i3D+1

            vtkarray3D = ah.array2vtk(bx)
            vtkarray3D.SetName(data.get('pointdata').GetPointData().GetArrayName(ivar))
            data3D.GetPointData().AddArray(vtkarray3D)

            vtkarray3D = ah.array2vtk(by)
            vtkarray3D.SetName(data.get('pointdata').GetPointData().GetArrayName(ivar+1))
            data3D.GetPointData().AddArray(vtkarray3D)
 
            vtkarray3D = ah.array2vtk(bz)
            vtkarray3D.SetName(data.get('pointdata').GetPointData().GetArrayName(ivar+2))
            data3D.GetPointData().AddArray(vtkarray3D)
  
            ivar = ivar+3
        else:
# Treat scalars:
            array2D = np.empty(nelem2D)
            array3D = np.empty(nelem3D)
            array2D = read.extract(data.get('pointdata'),data.get('pointdata').GetPointData().GetArrayName(ivar),attribute_mode='point')
            for i2D in range(nelem2D):
                for phi in phis:
                    array3D[i3D] = array2D[i2D]
                    i3D=i3D+1
            ivar = ivar+1

            vtkarray3D = ah.array2vtk(array3D)
            vtkarray3D.SetName(var)
            data3D.GetPointData().AddArray(vtkarray3D)
    
    data3D.Update()
    return {'pointdata': data3D, 'points': points3D}


def emiss(data,phi,theta,nu,alpha,recipe=1,polarized=1):

    print '=== Getting emissivities... ==='

# peel data:
    pointdata=data.get('pointdata')

# Synchrotron constants:
    
    c  = 2.99792458e10 # cm s^-1
    c1 = 6.27e18 # cm^-7/2 g^-5/2 s^4
    c2 = 2.368e-3 # g^-2 cm^-1 s^3
    me = 9.1093897e-28 # g

# Euclidian norm:
    norm = lambda x: np.sqrt(
    np.square(x[:,0]) + np.square(x[:,1]) + np.square(x[:,2])
    )
# Scalar product:
    scalar = lambda x,y: x[:,0]*y[:,0] + x[:,1]*y[:,1] + x[:,2]*y[:,2]

# Gamma from alpha:
    Gamma = 2.*alpha + 1.

# points:
    
    lfac=read.extract(pointdata,'lfac',attribute_mode='point')
    ux=read.extract(pointdata,'v1',attribute_mode='point')*lfac
    uy=read.extract(pointdata,'v2',attribute_mode='point')*lfac
    uz=read.extract(pointdata,'v3',attribute_mode='point')*lfac

    nelem = len(ux)


# get line of sight vector:
# phi and theta should be given in deg:
    n = np.empty((nelem,3))
    sintheta = np.sin(np.pi/180.*theta)
    costheta = np.cos(np.pi/180.*theta)
    sinphi   = np.sin(np.pi/180.*phi)
    cosphi   = np.cos(np.pi/180.*phi)
    n[:,0] = sintheta*cosphi
    n[:,1] = sintheta*sinphi
    n[:,2] = costheta


# get doppler factor:
    beta = np.empty((nelem,3))
    beta[:,0] = ux/lfac
    beta[:,1] = uy/lfac
    beta[:,2] = uz/lfac

    D = np.empty(nelem)
    D = 1./(
    lfac*(1.-(scalar(n,beta)))
    )

    print 'min and max Doppler factors: ', D.min(), D.max()

# get ndash:
    ndash = np.empty((nelem,3))
    ndash[:,0] = D*n[:,0] - (D+1.) * lfac/(lfac+1) * beta[:,0]
    ndash[:,1] = D*n[:,1] - (D+1.) * lfac/(lfac+1) * beta[:,1]
    ndash[:,2] = D*n[:,2] - (D+1.) * lfac/(lfac+1) * beta[:,2]

    print 'got ndash min and max:', ndash[:,0].min(), ndash[:,0].max(), ndash[:,1].min(), ndash[:,1].max(), ndash[:,2].min(), ndash[:,2].max()

# get bdash:
    if recipe ==1 or recipe ==2:
        bdash = get_bdash(pointdata)

        print 'got bdash min and max:', bdash[:,0].min(), bdash[:,0].max(), bdash[:,1].min(), bdash[:,1].max(), bdash[:,2].min(), bdash[:,2].max()

# bdash perp:
        bdashp = np.empty(nelem)
        tmp    = np.empty((nelem,3))
        tmp    = np.cross(ndash,bdash)
        bdashp = norm(tmp)

    else:
        bdashp = np.sqrt(read.extract(pointdata,'p',attribute_mode='point'))

    print 'got bdashp min and max: ', bdashp.min(), bdashp.max()

# emissivity:
    if recipe ==1:
        rhoe=read.extract(pointdata,'rhoe',attribute_mode='point')
        rhoe_0=read.extract(pointdata,'rhoe_0',attribute_mode='point')
    elif recipe ==2:
        rhoe=read.extract(pointdata,'n',attribute_mode='point')
        rhoe_0=read.extract(pointdata,'n0',attribute_mode='point')

    if recipe ==1 or recipe ==2:
# electron densities and cutoff:
        eps_inf=read.extract(pointdata,'eps_inf',attribute_mode='point')
        
        emiss = (c2*rhoe_0/me*c2*np.power(c1,(Gamma-3.)/2.)/(8.*np.pi) *
                 np.power(D,2.+(Gamma-1)/2.) *
                 powerorzero(bdashp,(Gamma+1.)/2.)*
        np.power(rhoe_0/rhoe,-((Gamma+2.)/3.)) *
        np.power(nu,(1.-Gamma)/2.) *
        powerorzero((1. - np.sqrt(nu)/(np.sqrt(c1*bdashp)*eps_inf)),Gamma-2.)
        )
# mask out the wind zone:
        m = lfac >9.
        emiss = np.ma.filled(np.ma.masked_array(emiss,m),0.)
    else:
        rho=read.extract(pointdata,'rho',attribute_mode='point')
        emiss = rho*np.power(D,2.+(Gamma-1)/2.)*powerorzero(bdashp,(Gamma+1.)/2.)
            
    if polarized !=1:
        print '=== Done with emissivities! ==='
    return {'emiss': emiss, 'points': data.get('points'),
            'phi': phi, 'theta': theta, 'nu': nu, 'alpha': alpha, 'recipe': recipe}
# The following only if polarized emission is requested:

# direction of emission:
    q       = np.empty((nelem,3))
    absb    = np.empty(nelem)
    b       = np.empty((nelem,3))
    ehat    = np.empty((nelem,3))
    absehat = np.empty((nelem))
    lhat    = np.empty((nelem,3))
    coschie = np.empty(nelem)
    sinchie = np.empty(nelem)
    cos2chie= np.empty(nelem)
    sin2chie= np.empty(nelem)
     
    b[:,0] =read.extract(pointdata,'b1',attribute_mode='point')
    b[:,1] =read.extract(pointdata,'b2',attribute_mode='point')
    b[:,2] =read.extract(pointdata,'b3',attribute_mode='point')

    absb = norm(b)

    q = b + np.cross(n,np.cross(beta,b))
    q[:,0] = q[:,0]/absb
    q[:,1] = q[:,1]/absb
    q[:,2] = q[:,2]/absb

    ehat      = np.cross(n,q) 
    absehat   = norm(ehat)
    ehat[:,0] = ehat[:,0]/absehat
    ehat[:,1] = ehat[:,1]/absehat
    ehat[:,2] = ehat[:,2]/absehat

    lhat[:,0] = -costheta*cosphi
    lhat[:,1] = -costheta*sinphi
    lhat[:,2] = sintheta

    coschie = scalar(lhat,ehat)
    sinchie = scalar(n,np.cross(lhat,ehat))

# complex phase with Doppelwinkelsatz:
    sin2chie = 2.*sinchie*coschie
    cos2chie = np.square(coschie)-np.square(sinchie)

# assemble complex polarized emissivity:
    piem = (Gamma+1)/(Gamma+7./3.)
    emisspol = np.ma.filled(np.ma.masked_array(piem * emiss * (cos2chie + 1j* sin2chie),absb==0),0.)

    print '=== Done with emissivities! ==='
    return {'emiss': emiss, 'emisspol': emisspol, 'points': data.get('points'),
            'phi': phi, 'theta': theta, 'nu': nu, 'alpha': alpha, 'recipe': recipe}



def powerorzero(x,y):
    m = x<=0.
    masked = np.ma.masked_array(np.power(x,y),m)
    return np.ma.filled(masked,0.)



def get_bdash(data):

    c  = 2.99792458e10 # cm s^-1

    bx=read.extract(data,'b1',attribute_mode='point')
    by=read.extract(data,'b2',attribute_mode='point')
    bz=read.extract(data,'b3',attribute_mode='point')
    ux=read.extract(data,'u1',attribute_mode='point')
    uy=read.extract(data,'u2',attribute_mode='point')
    uz=read.extract(data,'u3',attribute_mode='point')
    lfac=read.extract(data,'lfac',attribute_mode='point')
    
    bdash = np.empty((len(bx),3))
    bdash[:,0] = bx/lfac + ux/lfac*(bx*ux+by*uy+bz*uz)/c/c
    bdash[:,1] = by/lfac + uy/lfac*(bx*ux+by*uy+bz*uz)/c/c
    bdash[:,2] = bz/lfac + uz/lfac*(bx*ux+by*uy+bz*uz)/c/c

    return bdash


def show_emiss(inputdata,zi,nregrid,xrange=[-1.8e18,1.8e18],yrange=[-1.8e18,1.8e18]):
   
    fontsize=18
    nlevels=256

    points = inputdata.get('points')
    emiss = inputdata.get('emiss')

    x=points[:,0]
    y=points[:,1]
    z=points[:,2]

   # this is needed because griddata freaks out otherwise:
    smalldouble=1.e-8

    tmp=(xrange[1]-xrange[0])*smalldouble
    xi = np.linspace(xrange[0]-tmp,xrange[1]+tmp,nregrid[0])
    tmp=(yrange[1]-yrange[0])*smalldouble
    yi = np.linspace(yrange[0]-tmp,yrange[1]+tmp,nregrid[1])

    nelem = nregrid[0]*nregrid[1]
    xij=np.empty(nelem)
    yij=np.empty(nelem)
    zij=np.empty(nelem)
    zij[:]=zi
    for i in range(nregrid[0]):
        for j in range(nregrid[1]):
            xij[i+nregrid[0]*j] = xi[i]
            yij[i+nregrid[0]*j] = yi[j]

    data = griddata((x,y,z), emiss, (xij,yij,zij), method='nearest')

    data_2D = np.empty((nregrid[0],nregrid[1]))
    for i in range(nregrid[0]):
        for j in range(nregrid[1]):
            data_2D[j,i] = data[i+nregrid[0]*j]

    plotrange=np.array([data_2D.min(),data_2D.max()])

    plt.figure()
    plt.clf()
    cs = plt.contourf(xi,yi,data_2D,levels=np.arange(nlevels)/float(nlevels-1)*(plotrange[1]-plotrange[0])+plotrange[0])
#    cs = plt.contourf(xi,yi,data_2D,nlevels)
    cb = plt.colorbar()

    ax = plt.gca()
    ax.set_aspect('equal')
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    ax.xaxis.get_offset_text().set_fontsize(fontsize-2)
    ax.yaxis.get_offset_text().set_fontsize(fontsize-2)

    for tick in cb.ax.get_yticklabels():
        tick.set_fontsize(fontsize)


    plt.ylabel('y [cm]',size=fontsize)
    plt.xlabel('x [cm]',size=fontsize)
 
    plt.show()


def get_ray(data,nregrid,beg,end):


    points = data.get('points')
    emiss = data.get('emiss')

    x=points[:,0]
    y=points[:,1]
    z=points[:,2]

    xi = np.linspace(beg[0],end[0],nregrid)
    yi = np.linspace(beg[1],end[1],nregrid)
    zi = np.linspace(beg[2],end[2],nregrid)

    data = griddata((x,y,z), emiss, (xi,yi,zi), method='nearest')    

    return data

def shotgun(inputdata,npixel,L):

    print '=== Calculating radiation map... ==='
    points   = inputdata.get('points')
    emiss    = inputdata.get('emiss')
    if 'emisspol' in inputdata:
        emisspol = inputdata.get('emisspol')
    phi      = inputdata.get('phi')
    theta    = inputdata.get('theta')
    x=points[:,0]
    y=points[:,1]
    z=points[:,2]

    grid = make_grid(phi,theta,npixel,L)

    data    = griddata((x,y,z),emiss, (grid[0],grid[1],grid[2]),method='nearest')
    if 'emisspol' in inputdata:
        datap_r = griddata((x,y,z),emisspol.real, (grid[0],grid[1],grid[2]),method='nearest')
        datap_i = griddata((x,y,z),emisspol.imag, (grid[0],grid[1],grid[2]),method='nearest')
    print '=== Done with regridding, now integrating los ==='
    map    = np.zeros((npixel[0],npixel[1]))
    
    if 'emisspol' in inputdata:
        mapp_r = np.zeros((npixel[0],npixel[1]))
        mapp_i = np.zeros((npixel[0],npixel[1]))
        for i in range(npixel[0]):
            for j in range(npixel[1]):
                map[i,j] = map[i,j] + data[np.arange(npixel[2])+j*npixel[2]+i*npixel[2]*npixel[1]].sum()
                mapp_r[i,j] = mapp_r[i,j] + datap_r[np.arange(npixel[2])+j*npixel[2]+i*npixel[2]*npixel[1]].sum()
                mapp_i[i,j] = mapp_i[i,j] + datap_i[np.arange(npixel[2])+j*npixel[2]+i*npixel[2]*npixel[1]].sum()
    else:
        for i in range(npixel[0]):
            for j in range(npixel[1]):
                map[i,j] = map[i,j] + data[np.arange(npixel[2])+j*npixel[2]+i*npixel[2]*npixel[1]].sum()

    print '=== Done with map! ==='
    if 'emisspol' in inputdata:
        return {'Inu': map*L[2]/npixel[2],'Inupol': (mapp_r + 1j * mapp_i) * L[2]/npixel[2],
                'phi': inputdata.get('phi'), 'theta': inputdata.get('theta'),
            'nu': inputdata.get('nu'), 'alpha': inputdata.get('alpha'), 'recipe': inputdata.get('recipe'), 
            'npixel': npixel, 'L': L}
    else:
        return {'Inu': map*L[2]/npixel[2],
                'phi': inputdata.get('phi'), 'theta': inputdata.get('theta'),
            'nu': inputdata.get('nu'), 'alpha': inputdata.get('alpha'), 'recipe': inputdata.get('recipe'), 
            'npixel': npixel, 'L': L}

def shoot(inputdata,phi,theta,npixel,xgrid,ygrid):

# Euclidian norm:
    norm = lambda x: np.sqrt(np.square(x).sum())

    points = inputdata.get('points')
    emiss  = inputdata.get('emiss')
    x = np.linspace(xgrid[0],xgrid[1],npixel[0])
    y = np.linspace(ygrid[0],ygrid[1],npixel[1])


    dy = (xgrid[1]-xgrid[0])/npixel[0]
    dz = dy*np.sqrt(2.)

    beg = np.array([-1.e18,-1.e18,-2.e18])
    end = np.array([1.e18,-1.e18,0.])

    map = np.empty((npixel[0],npixel[1]))

    for i in range(npixel[0]):
        for j in range(npixel[1]):
            map[i,j] = get_ray(emiss,points,512,beg,end).sum()*norm(end-beg)/512.
            print 'got i,j:', i,j, i*npixel[0]+j,'map:',map[i,j]
            beg[1] = beg[1]+dy
            end[1] = end[1]+dy
        beg[1] = beg[1]-npixel[1]*dy
        end[1] = end[1]-npixel[1]*dy
        beg[2] = beg[2]+dz
        end[2] = end[2]+dz
        

    return map

def phase(map):
    '''map should be a 2D array of complex numbers'''
    phase = np.empty((map.shape[0],map.shape[1]))

    x = map.real
    y = map.imag

    for i in range(map[0].shape[0]):
        for j in range(map[0].shape[1]):
            if x[i,j] > 0.:
                phase[i,j] = np.arctan(y[i,j]/x[i,j])
            elif x[i,j] < 0. and y[i,j] >= 0.:
                phase[i,j] = np.arctan(y[i,j]/x[i,j]) + np.pi
            elif x[i,j] < 0. and y[i,j] < 0.:
                phase[i,j] = np.arctan(y[i,j]/x[i,j]) - np.pi
            elif x[i,j] == 0. and y[i,j] > 0.:
                 phase[i,j] = np.pi/2.
            elif x[i,j] == 0. and y[i,j] < 0.:
                phase[i,j] = -np.pi/2.
            else:
                phase[i,j] = float('nan') 

    return phase

def make_grid(phi,theta,npixel,L):

    xi,yi,zi=np.mgrid[0:npixel[0],0:npixel[1],0:npixel[2]]

    xi=xi/(npixel[0]-1.)
    yi=yi/(npixel[1]-1.)
    zi=zi/(npixel[2]-1.)

    xi = xi*L[0] - L[0]/2.
    yi = yi*L[1] - L[1]/2.
    zi = zi*L[2] - L[2]/2.

    xi = xi.flatten()
    yi = yi.flatten()
    zi = zi.flatten()

    nelem = len(xi)

# Now rotate the points around theta and then phi:

    cosphi = np.cos(phi*np.pi/180.)
    sinphi = np.sin(phi*np.pi/180.)
    costheta = np.cos(theta*np.pi/180.)
    sintheta = np.sin(theta*np.pi/180.)

    xid = xi
    yid = np.empty(nelem)
    zid = np.empty(nelem)
# Rotate by theta around the x-axis:
    yid = costheta*yi - sintheta*zi
    zid = sintheta*yi + costheta*zi

    xidd = np.empty(nelem)
    yidd = np.empty(nelem)
    zidd = zid
# Rotate by phi around the (old) z-axis:
    xidd =  cosphi*xid - sinphi*yid
    yidd =  sinphi*xid + cosphi*yid
    
    return [xidd,yidd,zidd]

def show_map(inputdata,npol='none',filename='none',background='w',plotrange='none',colpol='w',cmap=cm.gist_heat):

    nlevels=256
    fontsize=18

    L   = inputdata.get('L')
    map = inputdata.get('Inu')
    if 'Inupol' in inputdata:
        mappol = inputdata.get('Inupuol')


    xi = np.linspace(-L[0]/2.,L[0]/2.,map.shape[0])
    yi = np.linspace(-L[1]/2.,L[1]/2.,map.shape[1])
    xi2D=np.empty((map.shape[0],map.shape[1]))
    yi2D=np.empty((map.shape[0],map.shape[1]))
    for i in range(map.shape[1]):
        xi2D[:,i] = xi
    for i in range(map.shape[0]):
        yi2D[i,:] = yi
        
    plt.figure()
    plt.clf()

    if plotrange != 'none':
        levels = np.arange(nlevels)/float(nlevels-1)*(plotrange[1]-plotrange[0])+plotrange[0]
    	cs = plt.contourf(xi2D,yi2D,map,levels=levels,cmap=cmap)
    	cb = plt.colorbar()
    	cb.set_ticks(plotrange[1]*np.arange(10)/10.)
    else:
        cs = plt.contourf(xi2D,yi2D,np.log10(map),nlevels,cmap=cmap)
    	cb = plt.colorbar()

    ax = plt.gca()
    ax.set_aspect('equal')
    ax.patch.set_facecolor(background)
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    ax.xaxis.get_offset_text().set_fontsize(fontsize-2)
    ax.yaxis.get_offset_text().set_fontsize(fontsize-2)

    for tick in cb.ax.get_yticklabels():
        tick.set_fontsize(fontsize)

    plt.plot(0,0,'k+')

    plt.ylabel('b [cm]',size=fontsize)
    plt.xlabel('a [cm]',size=fontsize)

    plt.title(''.join((['phi=',str(inputdata.get('phi'))[0:4],' ; theta=',str(inputdata.get('theta'))[0:4]])))

    if npol != 'none' and 'Inupol' in inputdata:
        xbeg=-0.5*L[0]*float(npol[0]-1)/float(npol[0])
        ybeg=-0.5*L[1]*float(npol[1]-1)/float(npol[1])
        x=np.linspace(xbeg,-xbeg,npol[0])
        y=np.linspace(ybeg,-ybeg,npol[1])
        pi=abs(mappol)/map
        u=congrid(pi*np.sqrt(mappol).real/np.sqrt(abs(mappol)),npol,method='linear',centre=True)
        v=congrid(pi*np.sqrt(mappol).imag/np.sqrt(abs(mappol)),npol,method='linear',centre=True)
        plt.quiver(x,y,u,v,width=0.002,headlength=0,headwidth=0,headaxislength=0,pivot='middle',color=colpol)


    if filename == 'none':
        plt.show()
    else: 
        plt.savefig(''.join([filename,'.pdf']), format="pdf", transparent=False)
        plt.close()


def sweep_theta(offset,phi,theta_beg,theta_end,ntheta,nu,alpha,filenameout='data',type='pvtu'):

    npixel=[256,256,1024]
    L=[3e18,3e18,3e18]

    theta=np.linspace(theta_beg,theta_end,ntheta)

    data = get_pointdata(offset,filenameout=filenameout,type=type)
    for theta_i in theta:
        em = emiss(data,phi,theta_i,nu,alpha,polarized=0)
        map = shotgun(em,npixel,L)
        show_map(map,filename=''.join(['movie/t',str(offset).zfill(4),'theta',str(theta_i)[0:4]]))


def sweep_phi(offset,phi_beg,phi_end,theta,nphi,nu,alpha,filenameout='data',type='pvtu'):

    npixel=[256,256,512]
    L=[1e18,1e18,2e18]

    phi=np.linspace(phi_beg,phi_end,nphi)

    data = get_pointdata(offset,filenameout=filenameout,type=type)
    i = 0
    for phi_i in phi:
        em = emiss(data,phi_i,theta,nu,alpha)
        map = shotgun(em,npixel,L)
        show_map(map,filename=''.join(['movie/t',str(offset).zfill(4),'frame',str(i).zfill(4),'phi',str(phi_i)[0:4],'theta',str(theta)[0:4]]))
        show_map(map,filename=''.join(['movie/t',str(offset).zfill(4),'frame',str(i).zfill(4),'_pol_phi',str(phi_i)[0:4],'theta',str(theta)[0:4]]),npol=[32,32])
        i = i+1

def sweep_time(offset_beg,offset_end,phi,theta,nu,alpha,filenameout='data',type='pvtu'):

    npixel=[256,256,512]
    L=[1e18,1e18,2e18]

    offset=np.arange(offset_end-offset_beg+1)+offset_beg

    i = 0
    for offset_i in offset:
        data = get_pointdata(offset_i,filenameout=filenameout,type=type)
        em = emiss(data,phi,theta,nu,alpha)
        map = shotgun(em,npixel,L)
        show_map(map,filename=''.join(['movie_time/t',str(offset_i).zfill(4),'frame',str(i).zfill(4),'phi',str(phi)[0:4],'theta',str(theta)[0:4]]))
        show_map(map,filename=''.join(['movie_time/t',str(offset_i).zfill(4),'frame',str(i).zfill(4),'_pol_phi',str(phi)[0:4],'theta',str(theta)[0:4]]),npol=[32,32])
        i = i+1


def sweep_composite(offset_beg,offset_end,phi,theta,nu_r,nu_g,nu_b,alpha,rotate=0.,filenameout='data',type='pvtu'):

    npixel=[256,256,512]
    L=[3e18,3e18,3e18]

    offset=np.arange(offset_end-offset_beg+1)+offset_beg

    i = 0
    for offset_i in offset:
        data = get_pointdata(offset_i,filenameout=filenameout,type=type)
        em_r = emiss(data,phi,theta,nu_r,alpha,polarized=0)
        em_g = emiss(data,phi,theta,nu_g,alpha,polarized=0)
        em_b = emiss(data,phi,theta,nu_b,alpha,polarized=0)

        map_r = shotgun(em_r,npixel,L)
        map_g = shotgun(em_g,npixel,L)
        map_b = shotgun(em_b,npixel,L)

        composite(map_r,map_g,map_b,rotate=rotate,filename=''.join(['movie_composite/t',str(offset_i).zfill(4),'frame',str(i).zfill(4),'phi',str(phi)[0:4],'theta',str(theta)[0:4]]))

        i = i+1

def composite(map_r,map_g,map_b,rotate=0.0,filename='none'):
    
    img=np.zeros((map_r.get('Inu').shape[0],map_r.get('Inu').shape[1],3))

    img[:,:,0] = img_scale.sqrt(map_r.get('Inu'))
    img[:,:,1] = img_scale.sqrt(map_g.get('Inu'))
    img[:,:,2] = img_scale.sqrt(map_b.get('Inu'))

    img = scipy.ndimage.rotate(img, rotate, reshape=False,order=1) 

    fig=plt.figure(frameon=False)
    plt.imshow(img, aspect='equal',origin='lower')
    ax=fig.add_subplot(1,1,1)
    ax.set_axis_off()
    ax.set_aspect('equal')
    
    if filename == 'none':
        plt.show()
    else: 
        fig.set_size_inches(4,4)
        extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        plt.savefig(filename,bbox_inches=extent)
        plt.close()

def getFlux(offset,phi,theta,nu,alpha,nphi=None,filenameout='data',type='pvtu'):

    npixel=[128,128,128]
    L=[2e18,2e18,2e18]

    if nphi==None:
        data = get_pointdata(offset,filenameout=filenameout,type=type)
    else:
        data = rot2D(offset,nphi,filenameout=filenameout,type=type)
    em   = emiss(data,phi,theta,nu,alpha,polarized=0)
    map  = shotgun(em,npixel,L)

    flux = map.get('Inu').sum()

    return flux
