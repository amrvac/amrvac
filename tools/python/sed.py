import matplotlib.pyplot as plt
import numpy as np
from numpy.fft import fft2 as fft2
from numpy.fft import fftshift
import pickle
import read
from parmap import parmap as parmap
from scipy.interpolate import griddata

class handleSed():

    """Obtains power spectra from the simulations"""

    def __init__(self,offsets=None,filenameout='data',type='pvtu',xrange=None,yrange=None):
        

        self.offsetRange=offsets
        self.filenameout=filenameout
        self.type=type
        self.nx=1024
        self.ny=1024
        self.xrange=xrange
        self.yrange=yrange
        self.offset=0
        self.nkBins=None
        self.np=12

        if self.xrange != None:
            self.Lx = self.xrange[1] - self.xrange[0]
        else:
            self.Lx=None

        if self.yrange != None:
            self.Ly = self.yrange[1] - self.yrange[0] 
        else:
            self.Lx=None

    def readData(self):


        if self.type != 'vti':
            if self.xrange[0] <=0. :
                self.data=read.load(self.offset,file=self.filenameout,type=self.type,mirrorPlane=0)
            else:
                self.data=read.load(self.offset,file=self.filenameout,type=self.type)
        else:
            self.data=read.loadvti(self.offset,file=self.filenameout)
            self.nx=self.data.nx
            self.ny=self.data.nx


        if self.xrange == None:
            self.xrange = np.array(self.data.getBounds()[0:2])
        if self.yrange == None:
            self.yrange = np.array(self.data.getBounds()[2:4])

        self.Lx = self.xrange[1] - self.xrange[0]
        self.Ly = self.yrange[1] - self.yrange[0] 


    def getFFT(self,var,varname):

        # directly normalize to what the continuous FT would give:
        exec "self.%sFFT =  self.Lx/np.double(self.nx) * self.Ly/np.double(self.ny) * fftshift(fft2(var))" % (varname)



    def shellAverage(self,var):

        # get center indices:
        kxCenter = self.nx/2
        kyCenter = self.ny/2
        # kmax to be covered:
        kmax1D = min(kxCenter/self.Lx,kyCenter/self.Ly)
        # get the abs(k) value for the data:

        #kx,ky:        
        kx,ky = np.mgrid[0:self.nx,0:self.ny]
        self.kx = (kx - kxCenter)/self.Lx
        self.ky = (ky - kyCenter)/self.Ly
        #k:
        self.k = np.sqrt(self.kx**2+self.ky**2)
        self.getk1D()

        varFlat=var.flatten()
        kClip = np.clip(self.k,0,kmax1D).flatten()
        if self.nkBins == None:
            self.nkBins = min(kxCenter,kyCenter)
        ikBin = ((self.nkBins) * kClip/kmax1D).astype(int)
        dE = np.zeros(self.nkBins+1, dtype=complex)
        
        # just count the elements in bin, naively:
        ne = np.zeros(self.nkBins+1)
        for i in range(len(varFlat)):
            dE[ikBin[i]] = dE[ikBin[i]] + varFlat[i]
            ne[ikBin[i]] = ne[ikBin[i]] + 1.

        return dE[0:self.nkBins]/ne[0:self.nkBins] *2.*np.pi * self.k1D # \Delta E / \Delta k

    def ensembleAverage(self):

        avSpectrum = np.zeros(len(self.spectra[0]), dtype=complex)
        for spectrum in self.spectra:
            avSpectrum = avSpectrum + spectrum
        self.avSpectrum = avSpectrum / len(self.spectra)

        avDataFFT = np.zeros((self.nx,self.ny), dtype=complex)
        for data in self.dataFFT:
            avDataFFT = avDataFFT + data
        self.avDataFFT = avDataFFT / len(self.dataFFT)


    def regrid(self,var):

        CC=self.data.getCenterPoints()

        tmp0=np.complex(0,self.nx)
        self.x=np.linspace(self.xrange[0],self.xrange[1],self.nx)
        tmp1=np.complex(0,self.ny)
        self.y=np.linspace(self.yrange[0],self.yrange[1],self.ny)
        self.grid_x, self.grid_y = np.mgrid[self.xrange[0]:self.xrange[1]:tmp0, self.yrange[0]:self.yrange[1]:tmp1]
        gridvar = griddata(CC, var, (self.grid_x, self.grid_y), method='linear',fill_value=float(np.NaN))

        return gridvar

    def mirror(self,var,sign=1,side=1):
        varMirror = sign * np.append(var[::-1,::1],var,axis=0)
        return varMirror

    def addvar(self,varname):
        if self.type != 'vti':
            exec "self.%s = self.regrid(self.data.%s)" % (varname,varname)
        else:
            exec "self.%s = self.data.%s" % (varname,varname)

    def getSpectra(self,Quantity):

        self.offsets=np.arange(self.offsetRange[1]-self.offsetRange[0]+1)+self.offsetRange[0]

        if self.Lx == None or self.Ly == None: # We need to read the data once on the master to get the ranges.
            self.offset=self.offsets[0]
            self.readData()

        self.dataFFT = parmap(Quantity,self.offsets,np=self.np)
        print '========== got dataFFT =========='
        self.spectra = parmap(self.shellAverage,self.dataFFT,np=self.np)
        self.getk1D()

    def getk1D(self):
        # get center indices:
        kxCenter = np.double(self.nx)/2. / self.Lx
        kyCenter = np.double(self.ny)/2. / self.Ly
        if self.nkBins == None:
            self.nkBins = min(self.nx/2,self.ny/2)
        # kmax to be covered:
        kmax1D = min(kxCenter,kyCenter)
        dk = np.double(kmax1D)/np.double(self.nkBins)
        self.k1D=np.linspace(0,kmax1D,self.nkBins)+dk/2.
        self.Lmax1D = min(self.nx/2.,self.ny/2.) / kmax1D

    def Em(self,offset):
        
        self.offset=offset
        self.readData()
        if self.xrange[0] <=0. :
            self.data.reflectVar(self.data.b3)
        self.addvar('b3')
#        self.addvar('lfac')
#exclude wind zone:
#        iwind = self.lfac > 9
#        self.b3[iwind] = 0.
        self.getFFT(self.b3,'b3')
        self.EmFFT = 1./(8.*np.pi)*(self.b3FFT*self.b3FFT.conjugate())
        self.EmSpectrum = self.shellAverage(self.EmFFT)

        return self.EmFFT

    def Ek(self,offset):
        
        self.offset=offset
        self.readData()
        self.addvar('u1')
        self.addvar('u2')
#        self.addvar('u3')
        if self.xrange[0] <=0. :
            self.data.reflectVar(self.data.u1)
#            self.data.reflectVar(self.data.u3)
        self.getFFT(self.u1,'u1')
        self.getFFT(self.u2,'u2')
#        self.getFFT(self.u3,'u3')

        self.EkFFT = 0.5*(self.u1FFT*self.u1FFT.conjugate()
                           + self.u2FFT*self.u2FFT.conjugate())
#                           + self.u3FFT*self.u3FFT.conjugate()))
        self.EkSpectrum = self.shellAverage(self.EkFFT)

        return self.EkFFT

    def drho(self,offset):
        
        self.offset=offset
        self.readData()
        self.addvar('rho')
        self.rho=self.rho*1e18

        if self.xrange[0] <=0. :
            self.data.reflectVar(self.data.rho)
        self.getFFT(self.rho-self.rho.mean(),'rho')

        self.rhoFFT = np.sqrt(self.rhoFFT*self.rhoFFT.conjugate())
        self.rhoSpectrum = self.shellAverage(self.rhoFFT)

        return self.rhoFFT



    def EkFull(self,offset):
        
        self.offset=offset
        self.readData()
        self.addvar('u1')
        self.addvar('u2')
        self.addvar('rho')
#        self.rho=np.abs(self.rho)
        if self.xrange[0] <=0. :
            self.data.reflectVar(self.data.u1)
        self.getFFT(np.sqrt(self.rho)*self.u1,'u1')
        self.getFFT(np.sqrt(self.rho)*self.u2,'u2')

        self.EkFFT = 0.5*(self.u1FFT*self.u1FFT.conjugate()
                          + self.u2FFT*self.u2FFT.conjugate())
        self.EkSpectrum = self.shellAverage(self.EkFFT)

        return self.EkFFT

def getMagneticSpectra(offsets=[18,28],filenameout='data',xrange=[0e17,9e17],yrange=[-6e17,6e17],nx=512,ny=1024,np=6):
    se=handleSed(offsets=offsets,filenameout=filenameout,xrange=xrange,yrange=yrange)
    se.np=np
    se.nx=nx
    se.ny=ny
    se.getSpectra(se.Em)
    se.ensembleAverage()
    return se

def getVelocitySpectra(offsets=[18,28],filenameout='data',xrange=[3e17,9e17],yrange=[-6e17,6e17],nx=512,ny=1024,np=6):
    se=handleSed(offsets=offsets,filenameout=filenameout,xrange=xrange,yrange=yrange)
    se.np=np  
    se.nx=nx
    se.ny=ny
    se.getSpectra(se.EkFull)
    se.ensembleAverage()
    return se


def compareSpectra(offsets=[20,30],get=getMagneticSpectra,dumpfile='spectra.dat'):

#    folders=['45d_s1_2D','45d_s1_2Dhr','45d_s1_2Dvhr','45d_s1_2Duhr','45d_s1_2Dehr']
    folders=['45d_s1_2D','45d_s1_2Dhr','45d_s1_2Dvhr','45d_s1_2Duhr']
    nxs = np.array((100,200,400,800,1600))
    nys = (nxs * 12./9.).astype(int)
    npes = np.array((12,12,12,6,3))

    data=[]
    i=0
    for folder in folders:
        data.append(get(offsets=offsets,filenameout=''.join([folder,'/data']),nx=nxs[i],ny=nys[i],np=npes[i]))
        i=i+1

    dump(data,dumpfile=dumpfile)

    return data


def restore(dumpfile='spectra.dat'):

    file=open(dumpfile,'rb')
    data = pickle.load(file)
    file.close()
    return data

def dump(data,dumpfile='spectra.dat'):

    file=open(dumpfile,'wb') 
    pickle.dump(data,file, pickle.HIGHEST_PROTOCOL)
    file.close()
