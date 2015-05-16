import numpy as np
from scipy import interpolate
from scipy.integrate import odeint

# Euclidian norm:
norm = lambda x: np.sqrt(
np.square(x[:,0]) + np.square(x[:,1]) + np.square(x[:,2])
)

c  = 2.99792458e10 # cm s^-1
fopi=4.*np.pi


def sigma(d):

    b = np.empty((d.ncells,3))
    e = np.empty((d.ncells,3))
    v = np.empty((d.ncells,3))
    s = np.empty((d.ncells,3))
    
    b[:,0] = d.b1
    b[:,1] = d.b2
    b[:,2] = d.b3

    v[:,0] = d.u1/d.lfac
    v[:,1] = d.u2/d.lfac
    v[:,2] = d.u3/d.lfac

    e=np.cross(b,v)
    s=1./fopi*np.cross(e,b)

    sigma = norm(s)/(d.lfac**2*(d.rho*c**2+4.*d.p)*norm(v))

    return sigma

def v(d):
    return np.sqrt(d.u1**2+d.u2**2+d.u3**2)/d.lfac

def ptot(d):

    bx=d.b1
    by=d.b2
    bz=d.b3
    lfac=d.lfac
    ux=d.u1/lfac
    uy=d.u2/lfac
    uz=d.u3/lfac
    
    bdash = np.empty((d.ncells,3))
    bdash[:,0] = bx/lfac + lfac/(lfac+1.)*ux*(bx*ux+by*uy+bz*uz)/c/c
    bdash[:,1] = by/lfac + lfac/(lfac+1.)*uy*(bx*ux+by*uy+bz*uz)/c/c
    bdash[:,2] = bz/lfac + lfac/(lfac+1.)*uz*(bx*ux+by*uy+bz*uz)/c/c

    pb = (np.square(bdash[:,0])+np.square(bdash[:,1])+np.square(bdash[:,2]))/(2.*fopi)

    return d.p+pb

def get_bdash(d):

    bx=d.b1
    by=d.b2
    bz=d.b3
    lfac=d.lfac
    ux=d.u1/lfac
    uy=d.u2/lfac
    uz=d.u3/lfac
    
    d.bdash1 = bx/lfac + lfac/(lfac+1.)*ux*(bx*ux+by*uy+bz*uz)/c/c
    d.bdash2 = by/lfac + lfac/(lfac+1.)*uy*(bx*ux+by*uy+bz*uz)/c/c
    d.bdash3 = bz/lfac + lfac/(lfac+1.)*uz*(bx*ux+by*uy+bz*uz)/c/c

    print 'bdash appended to datastructure'



def aniso(d,dir=1):

    print 'just for 3D data, expects direction of slice as second parameter, default is one'
    if dir==1:
        bp2 = d.b2**2+d.b3**2
    if dir==2:
        bp2 = d.b1**2+d.b3**2
    btot = (d.b1**2+d.b2**2+d.b3**2)
    m=btot<1e-11
    return np.ma.filled(np.ma.masked_array(bp2/btot,m),0.)


def beta(d):
    '''
    returns the plasma beta assuming input data is in cgs units
    '''
    m=(d.b1**2+d.b2**2+d.b3**2)<=0.
    beta=np.ma.masked_array(8.*np.pi*d.p/(d.b1**2+d.b2**2+d.b3**2),m)
    return np.ma.filled(beta,beta.max())


def thetam(d,dir=1):
    '''
    returns the half opening angle of the fast Mach cone.
    Expects direction of slice as second parameter, default is one, set dir=-1 for 2.5D simulation
    '''
    bco2=(ptot(d)-d.p)*8.*np.pi
    if dir==1:
        vp2 = (d.u2**2+d.u3**2)/d.lfac**2
    if dir==2:
        vp2 = (d.u1**2+d.u3**2)/d.lfac**2
    if dir==-1:
        vp2 = (d.u1**2+d.u2**2)/d.lfac**2
    tmp=np.sqrt(bco2/(4.*np.pi*d.rho*vp2))/d.lfac
    m = tmp > 1.
    return np.ma.filled(np.ma.masked_array(np.arcsin(tmp),m),np.pi/2.)

def thetav(d,dir=1):
    '''
    returns the half opening angle of the flow.
    Expects direction of slice as second parameter, default is one, set dir=-1 for 2.5D simulation
    '''
    if dir==1:
        ur=d.u2
        uz=d.u3
    if dir==2:
        ur=d.u1
        uz=d.u3
    if dir==-1:
        ur=d.u1
        uz=d.u2
        
    return np.abs(np.arctan(ur/uz))

def em(d):
    '''
    returns the magnetic energy density
    '''
    
    return (d.b1**2+d.b2**2+d.b3**2)/(4.*np.pi) - 1./(8.*np.pi)*(
        (d.b1**2+d.b2**2+d.b3**2)/d.lfac**2 + 
        (d.b1*d.u1+d.b2*d.u2+d.b3*d.u3)**2/c**2/d.lfac**2)

def et(d):
    '''
    returns the thermal energy density
    '''
    
    return 4.*d.lfac**2*d.p-d.p

def ek(d):
    '''
    returns the kinetic energy density
    '''
    
    return (d.lfac-1) * d.lfac*d.rho*c**2

def er(d):
    '''
    returns the restmass energy density
    '''
    
    return d.lfac*d.rho*c**2

def cs(d):
    '''
    returns the sound-speed
    '''

    return np.sqrt(4./3.*d.p/((d.rho*c**2+4.*d.p)))*c


def streaminterp(u,v,x1,x2,var,x1_0,x2_0,dl=0.5):
    '''
    interpolates array of variables var onto the vector field u,v
    assumes evenly gridded data
    returns var_f, x1_f and x2_f
    '''

    xrange=[x1.min(),x1.max()]
    yrange=[x2.min(),x2.max()]

    # Normalize:
    mag=np.sqrt(u**2+v**2)
    u=u/mag
    v=v/mag
    u=np.nan_to_num(u)
    v=np.nan_to_num(v)

    fu = interpolate.RectBivariateSpline(x1, x2, u.transpose(),kx=1,ky=1,s=0)
    fv = interpolate.RectBivariateSpline(x1, x2, v.transpose(),kx=1,ky=1,s=0)
    fvar=[]
    for variable in var:
        fvar.append(interpolate.RectBivariateSpline(x1, x2, variable.transpose(),kx=1,ky=1,s=0))

# Go backwards:
    x1_f=[x1_0]
    x2_f=[x2_0]
    var_f=[list([]) for _ in var]
    for myvar,fvariable in zip(var_f,fvar):
        myvar.append(fvariable(x1_f[-1],x2_f[-1])[0][0])

    while(x1_f[-1] >= xrange[0] and x1_f[-1] <= xrange[1] and
    x2_f[-1] >= yrange[0] and x2_f[-1] <= yrange[1]):
        dt=-dl
        x1_f.append(x1_f[-1]+fu(x1_f[-1],x2_f[-1])[0][0]*dt)
        x2_f.append(x2_f[-1]+fv(x1_f[-1],x2_f[-1])[0][0]*dt)
        for variable,fvariable  in zip(var_f,fvar):
            variable.append(fvariable(x1_f[-1],x2_f[-1])[0][0])

    x1_f.pop()
    x2_f.pop()
    for variable in var_f:
        variable.pop()
        variable.reverse()
#    var_f.pop()
    x1_f.reverse()
    x2_f.reverse()
#    var_f.reverse()

# Go forwards:
    while(x1_f[-1] >= xrange[0] and x1_f[-1] <= xrange[1] and
    x2_f[-1] >= yrange[0] and x2_f[-1] <= yrange[1]):
        dt=dl
        x1_f.append(x1_f[-1]+fu(x1_f[-1],x2_f[-1])[0][0]*dt)
        x2_f.append(x2_f[-1]+fv(x1_f[-1],x2_f[-1])[0][0]*dt)
        for variable,fvariable  in zip(var_f,fvar):
            variable.append(fvariable(x1_f[-1],x2_f[-1])[0][0])


    x1_f.pop()
    x2_f.pop()
    for variable in var_f:
        variable.pop()
        variable=np.array(variable)
#    var_f.pop()

    return [np.array(x1_f),np.array(x2_f),var_f]

