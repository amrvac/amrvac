import sys
import math
import numpy as np
import matplotlib.pyplot as plt



fname='Blines_0000.txt'
f=open(fname,'r')
f.readline()
sss=f.readline()
dL=float(sss)
f.readline()
sss=f.readline()
nL=int(sss)
nLP=np.zeros(nL,dtype=int)
lengthL=np.zeros(nL)
f.readline()
nPmax=0
for iL in range(0,nL):
  sss=f.readline()
  nLP[iL]=int(sss.split()[1])
  lengthL[iL]=sss.split()[2]
  nPmax=max(nPmax,nLP[iL])

x=np.zeros((nL,nPmax))
y=np.zeros((nL,nPmax))
rho=np.zeros((nL,nPmax))
f.readline()
for iL in range(0,nL):
  for j in range(0,nLP[iL]):
    sss=f.readline()
    x[iL,j]=sss.split()[2]
    y[iL,j]=sss.split()[3]
    rho[iL,j]=sss.split()[4]
f.close()





fig=plt.figure(figsize=(10,10))

sub=plt.axes([0.1,0.1,0.8,0.8])
for iL in range(0,nL):
  colorl='r'
  if (iL>nL/2):
    colorl='b:'
  plt.plot(x[iL,0:nLP[iL]],y[iL,0:nLP[iL]],colorl,linewidth=1)
plt.xlim([-3,3])
plt.ylim([0,6])

plt.show()
