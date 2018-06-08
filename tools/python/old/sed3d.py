import matplotlib.pyplot as plt
import numpy as np
from numpy.fft import fftn as fft
from numpy.fft import fftshift
import pickle
import read
from parmap import parmap as parmap
from scipy.interpolate import griddata
from tvtk import array_handler as ah
import vtk as v


class handleSed():

    """Obtains power spectra from the 3D simulations"""

    def __init__(self,offsets=None,filenameout='data',type='pvtu',xrange=None,yrange=None,zrange=None):

        self.offsetRange=offsets
        self.filenameout=filenameout
        self.type=type
        self.nx=256
        self.ny=256
        self.nz=256
        self.xrange=xrange
        self.yrange=yrange
        self.zrange=zrange
        self.offset=0
        self.nkBins=None
        self.np=1

        if self.xrange != None:
            self.Lx = self.xrange[1] - self.xrange[0]
        else:
            self.Lx=None

        if self.yrange != None:
            self.Ly = self.yrange[1] - self.yrange[0] 
        else:
            self.Ly=None

        if self.zrange != None:
            self.Lz = self.zrange[1] - self.zrange[0] 
        else:
            self.Lz=None


    def readData(self):
      
        if self.type == 'vtu':
            filename=''.join([self.filenameout,repr(self.offset).zfill(4),'.vtu'])
            datareader = v.vtkXMLUnstructuredGridReader()
        elif self.type == 'pvtu':
            filename=''.join([self.filenameout,repr(self.offset).zfill(4),'.pvtu'])
            datareader = v.vtkXMLPUnstructuredGridReader()

        print '=== Reading ',filename,' ==='

        datareader.SetFileName(filename)
        datareader.Update()
        self.data = datareader.GetOutput()
        vtk_points=self.data.GetPoints().GetData()
        self.points=ah.vtk2array(vtk_points)


    def addvar(self,varname):

        print '=== Obtaining variable ',varname,' ==='

        exec "var = read.extract(self.data,'%s',attribute_mode='point')" % (varname)
        exec "self.%s = self.regrid(var)" % varname

    def regrid(self,var):

        CC=self.points
        self.makeGrid()
        gridvar = griddata(CC, var, (self.grid_x, self.grid_y, self.grid_z), method='nearest',fill_value=float(np.NaN))

        return gridvar

    def makeGrid(self):

        tmp0=np.complex(0,self.nx)
        tmp1=np.complex(0,self.ny)
        tmp2=np.complex(0,self.nz)
        self.x=np.linspace(self.xrange[0],self.xrange[1],self.nx)
        self.y=np.linspace(self.yrange[0],self.yrange[1],self.ny)
        self.z=np.linspace(self.zrange[0],self.zrange[1],self.nz)
        self.grid_x, self.grid_y, self.grid_z = np.mgrid[self.xrange[0]:self.xrange[1]:tmp0, self.yrange[0]:self.yrange[1]:tmp1, self.zrange[0]:self.zrange[1]:tmp2]



    def getFFT(self,var,varname):
        print '=== FFT of variable ',varname,' ==='
        # directly normalize to what the continuous FT would give:
        exec "self.%sFFT =  self.Lx/np.double(self.nx) * self.Ly/np.double(self.ny) * self.Lz/np.double(self.nz) * fftshift(fft(var))" % (varname)

    def Em(self,offset):
        
        self.offset=offset
        self.readData()
        self.addvar('b1')
        self.getFFT(self.b1,'b1')
        self.addvar('b2')
        self.getFFT(self.b2,'b2')
        self.addvar('b3')
        self.getFFT(self.b3,'b3')
        self.EmFFT = 1./(8.*np.pi)*(self.b1FFT*self.b1FFT.conjugate()
                                    + self.b2FFT*self.b2FFT.conjugate()
                                    + self.b3FFT*self.b3FFT.conjugate()
                                    )
        self.EmSpectrum = self.shellAverage(self.EmFFT)

        return self.EmFFT

    def Ek(self,offset):
        
        self.offset=offset
        self.readData()
        self.addvar('u1')
        self.getFFT(self.u1,'u1')
        self.addvar('u2')
        self.getFFT(self.u2,'u2')
        self.addvar('u3')
        self.getFFT(self.u3,'u3')

        self.EkFFT = (0.5*(self.u1FFT*self.u1FFT.conjugate()
                           + self.u2FFT*self.u2FFT.conjugate()
                           + self.u3FFT*self.u3FFT.conjugate()))
        self.EkSpectrum = self.shellAverage(self.EkFFT)

        return self.EkFFT

    def generateSpectrum(self,k0,alpha):
        

        # get center indices:
        kxCenter = self.nx/2
        kyCenter = self.ny/2
        kzCenter = self.nz/2
        # kmax to be covered:
        kmax1D = min(kxCenter/self.Lx,kyCenter/self.Ly,kzCenter/self.Lz)
        # get the abs(k) value for the data:

        #kx,ky,kz:        
        kx,ky,kz = np.mgrid[0:self.nx,0:self.ny,0:self.nz]
        self.kx = (kx - kxCenter)/self.Lx
        self.ky = (ky - kyCenter)/self.Ly
        self.kz = (kz - kzCenter)/self.Lz
        #k:
        self.k = np.sqrt(self.kx**2+self.ky**2+self.kz**2)

        phase1 = np.random.rand(self.nx,self.ny,self.nz)
        phase2 = np.random.rand(self.nx,self.ny,self.nz)
        phase3 = np.random.rand(self.nx,self.ny,self.nz)

        self.b1FFT = (self.k/k0)**(-alpha/2.-1.) * np.exp(2.*np.pi*1j*phase1)
        self.b2FFT = (self.k/k0)**(-alpha/2.-1.) * np.exp(2.*np.pi*1j*phase2)
        self.b3FFT = (self.k/k0)**(-alpha/2.-1.) * np.exp(2.*np.pi*1j*phase3)

        self.b1FFT[self.k<k0] = 0.
        self.b2FFT[self.k<k0] = 0.
        self.b3FFT[self.k<k0] = 0.

        self.b1 = np.fft.ifftn(np.fft.ifftshift(self.b1FFT))
        self.b2 = np.fft.ifftn(np.fft.ifftshift(self.b2FFT))
        self.b3 = np.fft.ifftn(np.fft.ifftshift(self.b3FFT))

    def shellAverage(self,var):
        print '=== Shell averaging spectra ==='

        # get center indices:
        kxCenter = self.nx/2
        kyCenter = self.ny/2
        kzCenter = self.nz/2
        # kmax to be covered:
        kmax1D = min(kxCenter/self.Lx,kyCenter/self.Ly,kzCenter/self.Lz)
        # get the abs(k) value for the data:

        #kx,ky,kz:        
        kx,ky,kz = np.mgrid[0:self.nx,0:self.ny,0:self.nz]
        self.kx = (kx - kxCenter)/self.Lx
        self.ky = (ky - kyCenter)/self.Ly
        self.kz = (kz - kzCenter)/self.Lz
        #k:
        self.k = np.sqrt(self.kx**2+self.ky**2+self.kz**2)
        self.getk1D()

        varFlat=var.flatten()
        kClip = np.clip(self.k,0,kmax1D).flatten()
        if self.nkBins == None:
            self.nkBins = min(kxCenter,kyCenter,kzCenter)
        ikBin = ((self.nkBins) * kClip/kmax1D).astype(int)
        dE = np.zeros(self.nkBins+1, dtype=complex)
        
        # just count the elements in bin, naively:
        ne = np.zeros(self.nkBins+1)
        for i in range(len(varFlat)):
            dE[ikBin[i]] = dE[ikBin[i]] + varFlat[i]
            ne[ikBin[i]] = ne[ikBin[i]] + 1.

        return dE[0:self.nkBins]/ne[0:self.nkBins] *4.*np.pi * self.k1D**2 # \Delta E / \Delta k


    def getk1D(self):
        # get center indices:
        kxCenter = np.double(self.nx)/2. / self.Lx
        kyCenter = np.double(self.ny)/2. / self.Ly
        kzCenter = np.double(self.nz)/2. / self.Lz
        if self.nkBins == None:
            self.nkBins = min(self.nx/2,self.ny/2,self.nz/2)
        # kmax to be covered:
        kmax1D = min(kxCenter,kyCenter,kzCenter)
        dk = np.double(kmax1D)/np.double(self.nkBins)
        self.k1D=np.linspace(0,kmax1D,self.nkBins)+dk/2.
        self.Lmax1D = min(self.nx/2.,self.ny/2.,self.nz/2.) / kmax1D

    def getSpectra(self,Quantity):

        self.offsets=np.arange(self.offsetRange[1]-self.offsetRange[0]+1)+self.offsetRange[0]

        if self.Lx == None or self.Ly == None or self.Lz == None: # We need to read the data once on the master to get the ranges.
            self.offset=self.offsets[0]
            self.readData()

        self.dataFFT = parmap(Quantity,self.offsets,np=self.np)
        print '========== got dataFFT =========='
        self.spectra = parmap(self.shellAverage,self.dataFFT,np=self.np)
        self.getk1D()

    def ensembleAverage(self):

        avSpectrum = np.zeros(len(self.spectra[0]), dtype=complex)
        for spectrum in self.spectra:
            avSpectrum = avSpectrum + spectrum
        self.avSpectrum = avSpectrum / len(self.spectra)

        avDataFFT = np.zeros((self.nx,self.ny,self.nz), dtype=complex)
        for data in self.dataFFT:
            avDataFFT = avDataFFT + data
        self.avDataFFT = avDataFFT / len(self.dataFFT)


def getMagneticSpectra(offsets=[40,50],filenameout='dataRG',xrange=[-1.5e18,1.5e18],yrange=[-1.5e18,1.5e18],zrange=[-1.5e18,1.5e18]):
    se=handleSed(offsets=offsets,filenameout=filenameout,xrange=xrange,yrange=yrange,zrange=zrange)
    se.np=2
    se.nx=300
    se.ny=300
    se.nz=300
    se.getSpectra(se.Em)
    se.ensembleAverage()
    return se

def getVelocitySpectra(offsets=[40,50],filenameout='dataRG',xrange=[-1.5e18,1.5e18],yrange=[-1.5e18,1.5e18],zrange=[-1.5e18,1.5e18]):
    se=handleSed(offsets=offsets,filenameout=filenameout,xrange=xrange,yrange=yrange,zrange=zrange)
    se.np=2
    se.nx=300
    se.ny=300
    se.nz=300
    se.getSpectra(se.Ek)
    se.ensembleAverage()
    return se


def compareSpectra(offsets=[40,50],filenameout='dataRG',get=getMagneticSpectra,dumpfile='spectra.dat'):

    folders=['.']

    data=[]
    for folder in folders:
        data.append(get(offsets=offsets,filenameout=''.join([folder,'/',filenameout])))

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
