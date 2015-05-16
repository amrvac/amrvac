import read, helper
import numpy as np
import matplotlib.pyplot as plt

''' Converting the vti output into an easily readable 1D ascii-file '''

fontsize=10
dpi=300

filenames   = ['output/data1D','output/data3D','output/data3DFCT']
labels      = [r'1D', r'3D GLM', r'3D, FCT']
markers     = ['ro','gs','bD']
filenameout = 'cuts_'


#slice 3D data at these indices:
ix=16
iy=16

simulations=[]
for filename in filenames:
    simulations.append(read.loadvti(1,file=filename))


for (filename,sim) in zip(filenames,simulations):
    outarr=np.empty((sim.nx,len(sim.getVarnames())+1))
    varnames='x,  '
    outarr[:,0] = sim.x[:]
    i=1
    for varname in sim.getVarnames():
        varnames=varnames+varname+',  '
        if sim.ndim==3:
# We have to extract one line in x-direction
            exec "outarr[:,i]=sim.%s[:,ix,iy]" % (varname)
        else:
            exec "outarr[:,i]=sim.%s[:]" % (varname)
        i=i+1
    np.savetxt(filename+'.csv',outarr,delimiter=',  ', header=varnames)

    


# Also plot the data, sorry for the spaghetti code:

# rho:
fig1=plt.figure(figsize=(3.5,3.5))
ax1  = fig1.add_subplot(1,1,1)

for (marker,label,sim) in zip(markers,labels,simulations):
    if sim.ndim==3:
# We have to extract one line in x-direction
        y=sim.rho[:,ix,iy]
    else:
        y=sim.rho[:]
    ax1.plot(sim.x,y,marker,label=label,markersize=3,fillstyle='none')

leg = ax1.legend(loc='best', fancybox=True,prop={'size':fontsize}, handlelength=1, numpoints=1)
leg.set_alpha(0.8)

ax1.set_xlabel(r'$x$',fontsize=fontsize)
ax1.set_ylabel(r'$\rho$',fontsize=fontsize)

fig1.savefig(filenameout+'rho.pdf', transparent=False,aa=True,dpi=dpi,bbox_inches='tight')
plt.close()


# p:
fig1=plt.figure(figsize=(3.5,3.5))
ax1  = fig1.add_subplot(1,1,1)

for (marker,label,sim) in zip(markers,labels,simulations):
    if sim.ndim==3:
# We have to extract one line in x-direction
        y=sim.p[:,ix,iy]
    else:
        y=sim.p[:]
    ax1.plot(sim.x,y,marker,label=label,markersize=3,fillstyle='none')

leg = ax1.legend(loc='best', fancybox=True,prop={'size':fontsize}, handlelength=1, numpoints=1)
leg.set_alpha(0.8)

ax1.set_xlabel(r'$x$',fontsize=fontsize)
ax1.set_ylabel(r'$p$',fontsize=fontsize)

fig1.savefig(filenameout+'p.pdf', transparent=False,aa=True,dpi=dpi,bbox_inches='tight')
plt.close()


# vx:
fig1=plt.figure(figsize=(3.5,3.5))
ax1  = fig1.add_subplot(1,1,1)

for (marker,label,sim) in zip(markers,labels,simulations):
    if sim.ndim==3:
# We have to extract one line in x-direction
        y=sim.u1[:,ix,iy]/sim.lfac[:,ix,iy]
    else:
        y=sim.u1[:]/sim.lfac[:]
    ax1.plot(sim.x,y,marker,label=label,markersize=3,fillstyle='none')

leg = ax1.legend(loc='best', fancybox=True,prop={'size':fontsize}, handlelength=1, numpoints=1)
leg.set_alpha(0.8)

ax1.set_xlabel(r'$x$',fontsize=fontsize)
ax1.set_ylabel(r'$v_x$',fontsize=fontsize)

fig1.savefig(filenameout+'vx.pdf', transparent=False,aa=True,dpi=dpi,bbox_inches='tight')
plt.close()


# vy:
fig1=plt.figure(figsize=(3.5,3.5))
ax1  = fig1.add_subplot(1,1,1)

for (marker,label,sim) in zip(markers,labels,simulations):
    if sim.ndim==3:
# We have to extract one line in x-direction
        y=sim.u2[:,ix,iy]/sim.lfac[:,ix,iy]
    else:
        y=sim.u2[:]/sim.lfac[:]
    ax1.plot(sim.x,y,marker,label=label,markersize=3,fillstyle='none')

leg = ax1.legend(loc='best', fancybox=True,prop={'size':fontsize}, handlelength=1, numpoints=1)
leg.set_alpha(0.8)

ax1.set_xlabel(r'$x$',fontsize=fontsize)
ax1.set_ylabel(r'$v_y$',fontsize=fontsize)

fig1.savefig(filenameout+'vy.pdf', transparent=False,aa=True,dpi=dpi,bbox_inches='tight')
plt.close()


# Gamma:
fig1=plt.figure(figsize=(3.5,3.5))
ax1  = fig1.add_subplot(1,1,1)

for (marker,label,sim) in zip(markers,labels,simulations):
    if sim.ndim==3:
# We have to extract one line in x-direction
        y=sim.lfac[:,ix,iy]
    else:
        y=sim.lfac[:]
    ax1.plot(sim.x,y,marker,label=label,markersize=3,fillstyle='none')

leg = ax1.legend(loc='best', fancybox=True,prop={'size':fontsize}, handlelength=1, numpoints=1)
leg.set_alpha(0.8)

ax1.set_xlabel(r'$x$',fontsize=fontsize)
ax1.set_ylabel(r'$\Gamma$',fontsize=fontsize)

fig1.savefig(filenameout+'lfac.pdf', transparent=False,aa=True,dpi=dpi,bbox_inches='tight')
plt.close()


# by:
fig1=plt.figure(figsize=(3.5,3.5))
ax1  = fig1.add_subplot(1,1,1)

for (marker,label,sim) in zip(markers,labels,simulations):
    if sim.ndim==3:
# We have to extract one line in x-direction
        y=sim.b2[:,ix,iy]
    else:
        y=sim.b2[:]
    ax1.plot(sim.x,y,marker,label=label,markersize=3,fillstyle='none')

leg = ax1.legend(loc='best', fancybox=True,prop={'size':fontsize}, handlelength=1, numpoints=1)
leg.set_alpha(0.8)

ax1.set_xlabel(r'$x$',fontsize=fontsize)
ax1.set_ylabel(r'$B_y$',fontsize=fontsize)

fig1.savefig(filenameout+'by.pdf', transparent=False,aa=True,dpi=dpi,bbox_inches='tight')
plt.close()
