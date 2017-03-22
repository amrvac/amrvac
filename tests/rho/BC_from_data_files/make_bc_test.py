import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d, interp1d

for t in range(0,20): # 20 time snapshots

    Nx=20
    Ny=20
    dx=1./Nx
    dy=1./Ny
    x=np.zeros(Nx)
    y=np.zeros(Ny)
    rho1=np.zeros((Nx,Ny))
    oneline=' {: 0.6E} \n'
    output=open('bc/bc_t'+"%02d" % (t)+'_AMR1.dat','w+')
    for i in range(Nx):
       x[i]=dx/2.+dx*i
       for j in range(Ny):
        y[j]=dy/2.+dy*j
        rho1[i,j]=random.uniform(0.,2.)
        output.write(oneline.format(rho1[i,j]))
    output.close()

    xx,yy=np.meshgrid(np.append(x,x[Nx-1]+x[1]-x[0]),np.append(y,y[Nx-1]+y[1]-y[0]),indexing='ij')
    plt.pcolor(xx,yy,rho1)
    plt.savefig('bc/bc_t'+"%02d" % (t)+'_AMR1.png')
    plt.close()

    # yk=np.linspace(y[0],y[Ny-1],num=np.shape(pif)[0],endpoint=True) # mesh corresponding to the rotated signal, too many pixels
    mxnest=2
    rho2=np.zeros((Nx*2,Ny*2))
    # rho3=np.zeros((Nx*4,Ny*4))
    # for i in range(2,mxnest+1):
    Nx=2*Nx
    Ny=2*Ny
    dx=dx/2.
    dy=dy/2.
    f=interp2d(x,y,rho1,kind='linear') # former signal on former mesh
    x=np.linspace(x[0]/2,x[Nx/2-1]+(x[Nx/2-1]-x[Nx/2-2])/4.,num=Nx,endpoint=True)
    y=np.linspace(y[0]/2,y[Ny/2-1]+(y[Ny/2-1]-y[Ny/2-2])/4.,num=Ny,endpoint=True)
    rho2=f(x,y) # interpolate the rotated signal on the mesh y of the initial size
    output=open('bc/bc_t'+"%02d" % (t)+'_AMR2.dat','w+')
    for i in range(Nx):
       x[i]=dx/2.+dx*i
       for j in range(Ny):
           y[j]=dy/2.+dy*j
           output.write(oneline.format(rho2[i,j]))
    output.close()

    print rho1[0,1], rho1[1,1]
    print rho1[0,0], rho1[1,0]
    print
    print rho2[0,1], rho2[1,1]
    print rho2[0,0], rho2[1,0]
    xx,yy=np.meshgrid(np.append(x,x[Nx-1]+x[1]-x[0]),np.append(y,y[Nx-1]+y[1]-y[0]),indexing='ij')
    plt.pcolor(xx,yy,rho2)
    plt.savefig('bc/bc_t'+"%02d" % (t)+'_AMR2.png')
    plt.close()
