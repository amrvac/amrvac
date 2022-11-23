import numpy as np
import matplotlib.pyplot as plt


fname='DEM_0000.txt'
f=open(fname,'r')
f.readline()
sss=f.readline()
nTe=int(sss)
Telog=np.zeros(nTe)
DEM=np.zeros(nTe)

f.readline()
for iTe in range(0,nTe):
  sss=f.readline()
  Telog[iTe]=sss.split()[1]
  DEM[iTe]=sss.split()[2]
f.close()




fig=plt.figure()

sub=plt.subplot(111)
plt.step(Telog,DEM,where='post',color='k')
plt.yscale('log')
plt.xlabel(r'$\rm log_{10}\ T\ (K)$')
plt.ylabel('EM [cm$^{-5}$]')

plt.show()


