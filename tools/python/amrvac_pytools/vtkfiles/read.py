"""Modules used for reading AMRVAC data"""
#=============================================================================
import numpy as np
import vtk as v
from numpy import loadtxt as pylabload
import sys, time
import struct

from amrvac_pytools.vtkfiles import numpy_support as ah


if sys.platform == "win32":
# On Windows, the best timer is time.clock()
    default_timer = time.clock
else:
# On most other platforms the best timer is time.time()
    default_timer = time.time


#=============================================================================
def extract(data,varname,attribute_mode='cell'):
    """Extracts variable "varname" from vtkfiles datastructure 'data'"""
    if attribute_mode == 'cell':
        vtk_values = data.GetCellData().GetArray(varname)
    elif attribute_mode == 'point':
        if data.GetPointData().GetNumberOfArrays() > 0:
            vtk_values = data.GetPointData().GetArray(varname)
        else:
            # Convert to pointdata first
            c2p = v.vtkCellDataToPointData()
            c2p.SetInput(data)
            pointdata=c2p.GetOutput()
            pointdata.Update()
            vtk_values = pointdata.GetPointData().GetScalars(varname)
    elif attribute_mode == 'topoint':
        c2p = v.vtkCellDataToPointData()
        c2p.SetInput(data)
        pointdata=c2p.GetOutput()
        pointdata.Update()
        vtk_values = pointdata.GetPointData().GetScalars(varname)
    else:
        print("attribute_mode is either 'cell' or 'point'")
    
    return ah.vtk2array(vtk_values)

#=============================================================================
class load_vtkfile():
    """
    Loader class for vtu and pvtu files.
    """
    def __init__(self,offset,get=1,file='data',type='vtu',mirrorPlane=None,
                 rotateX=0,rotateY=0,rotateZ=0,
                 scaleX=1,scaleY=1,scaleZ=1,
                 silent=0):
        self.offset=offset
        self.filenameout = file
        self.type = type
        self.isLoaded = False
        self.mirrorPlane=mirrorPlane
        self.silent = silent
        
        self.rotateX = rotateX
        self.rotateY = rotateY
        self.rotateZ = rotateZ

        self.scaleX = scaleX
        self.scaleY = scaleY
        self.scaleZ = scaleZ

        if type == 'vtu':
            self.filename=''.join([self.filenameout,repr(offset).zfill(4),'.vtu'])
            self.datareader = v.vtkXMLUnstructuredGridReader()
        elif type == 'pvtu':
            self.filename=''.join([self.filenameout,repr(offset).zfill(4),'.pvtu'])
            self.datareader = v.vtkXMLPUnstructuredGridReader()
        elif type == 'vti':
            self.filename=''.join([self.filenameout,repr(offset).zfill(4),'.vti'])
            self.datareader = v.vtkXMLImageDataReader()
        else:
            print('Unknown filetype')
            
        if (self.silent == 0): print('========================================')
        if (self.silent == 0): print('loading file %s' % (self.filename))

        if get is not None:
            self.getAll()


    def getTime(self):
        try:
            self.time=ah.vtk2array(self.data.GetFieldData().GetArray(0))[0]
        except AttributeError:
            self.time = np.nan
        return self.time


    def getVar(self,varname):
        try:
            tmp = None
            exec("tmp=self.%s" % (varname))
            return tmp
        except AttributeError:
            print("Unknown variable", varname)
    
    def getData(self):
        self.datareader.SetFileName(self.filename)
        self.datareader.Update()
        self.data = self.datareader.GetOutput()
        self.ncells = self.data.GetNumberOfCells()
        self.isLoaded = True

        if (self.rotateX != 0 or self.rotateY != 0 or self.rotateZ != 0
            or self.scaleX != 1 or self.scaleY != 1 or self.scaleZ != 1):
            transform = v.vtkTransform()
            transform.RotateX(self.rotateX)
            transform.RotateY(self.rotateY)
            transform.RotateZ(self.rotateZ)
            transform.Scale(self.scaleX, self.scaleY, self.scaleZ)
            transfilter = v.vtkTransformFilter()
            transfilter.SetTransform(transform)
            transfilter.SetInputData(self.data)
            self.data = transfilter.GetOutput()
            transfilter.Update()


    def getPointData(self):
        self.getData()
        if self.data.GetPointData().GetNumberOfArrays() == 0:
            c2p = v.vtkCellDataToPointData()
            c2p.SetInputData(self.data)
            self.pointdata=c2p.GetOutput()
            c2p.Update()
        else:
            self.pointdata=self.data


    def getVars(self):
        nvars= self.data.GetCellData().GetNumberOfArrays()
        for i in range(nvars):
            varname = self.data.GetCellData().GetArrayName(i)
            if (self.silent == 0): print("Assigning variable:", varname)
            vtk_values = self.data.GetCellData().GetArray(varname)
            exec("self.%s = ah.vtk2array(vtk_values)[0:self.ncells]" % (varname))


    def getVarnames(self):
        nvars= self.data.GetCellData().GetNumberOfArrays()
        varnames=[]
        for i in range(nvars):
            varnames.append(self.data.GetCellData().GetArrayName(i))
        return varnames


    def getBounds(self):
        return self.data.GetBounds()


    def getVert(self,icell):
        if self.data.GetCell(icell).GetCellType() == 8 :
            pts=ah.vtk2array(self.data.GetCell(icell).GetPoints().GetData())
            return np.array((pts[0][0:2],pts[1][0:2],pts[3][0:2],pts[2][0:2]))
        if self.data.GetCell(icell).GetCellType() == 3 :
            pts=ah.vtk2array(self.data.GetCell(icell).GetPoints().GetData())
            return np.array((pts[0][0],pts[1][0]))
        else: 
            if (self.silent == 0): print("Can handle only type 3 or type 8")
        

    def getPointList(self):
        tstart = default_timer()
        try:
            self.data
        except AttributeError:
            self.getData()
        try:
            [self.xlist,self.ylist]
        except AttributeError:
            if self.data.GetCell(0).GetCellType() != 8 :
                if (self.silent == 0): print("Can handle pixel types only")
                pass
            self.xlist = []
            self.ylist = []
            for icell in range(self.ncells):
                pts=ah.vtk2array(self.data.GetCell(icell).GetPoints().GetData())
                self.xlist.extend((pts[0][0],pts[1][0],pts[3][0],pts[2][0],None))
                self.ylist.extend((pts[0][1],pts[1][1],pts[3][1],pts[2][1],None))
        tend = default_timer()
        if (self.silent == 0): print('Getting formatted pointlist time=%f sec' % (tend-tstart))
        return [self.xlist,self.ylist]        


    def getCenterPoints(self):
        tstart = default_timer()
        firstget = False
        try:
            self.data
        except AttributeError:
            self.getData()
            firstget = True
        try:
            self.centerpoints
        except AttributeError:
            if self.getNdim() == 2 or self.getNdim() == 3:
                self.centerpoints=np.empty((self.ncells,2))
                for icell in range(self.ncells):
                    vert=self.getVert(icell)
                    self.centerpoints[icell,0]=vert[:,0].mean()
                    self.centerpoints[icell,1]=vert[:,1].mean()
            if self.getNdim() == 1 :
                self.centerpoints=np.empty((self.ncells))
                for icell in range(self.ncells):
                    vert=self.getVert(icell)
                    self.centerpoints[icell]=vert.mean()
            tend = default_timer()
        if firstget:
            if (self.silent == 0): print('Getting cell center coordiantes time=%f sec' % (tend-tstart))
        return self.centerpoints


    def getSurface(self):
        def calcSurface(icell):
            vert=self.getVert(icell)
            return np.sqrt(((vert[0][0]-vert[1][0])**2+(vert[0][1]-vert[1][1])**2)*
                                              ((vert[1][0]-vert[2][0])**2+(vert[1][1]-vert[2][1])**2))

        tstart = default_timer()
        try:
            self.data
        except AttributeError:
            self.getData()
        try:
            self.surface
        except AttributeError:
            # Assuming cells are rectangles here.
            surface=list(map(calcSurface,list(range(self.ncells))))
            self.surface=np.array(surface)
        tend = default_timer()
        if (self.silent == 0): print('Getting cell surface (assuming rectangles) time=%f sec' % (tend-tstart))
        return self.surface

    def getPieces(self,nBWidth,nBHeight,nIgnore):
        tstart = default_timer()
        try:
            self.data
        except AttributeError:
            self.getData()
        try:
            [self.xBlockList,self.yBlockList]
        except AttributeError:

            self.xBlockList=[]
            self.yBlockList=[]

            pts=self.getPoints()
            nPoints = np.shape(pts)[0]

            # Calculate number of blocks
            nBSize=(nBWidth+1)*(nBHeight+1)
            nBlocks=int(nPoints/nBSize)

            # For all blocks, draw the perimeter
            for iBlock in range(nBlocks):
                start=iBlock*nBSize
                end  =(iBlock+1)*nBSize-1
                # Side 1
                for i in range(start,start+nBWidth+1):
                    self.xBlockList.append(pts[i][0])
                    self.yBlockList.append(pts[i][1])
                # Side 2
                for i in range(start+nBWidth,end+1,nBWidth+1):
                    self.xBlockList.append(pts[i][0])
                    self.yBlockList.append(pts[i][1])
                # Side 3
                for i in range(end,end-nBWidth-1,-1):
                    self.xBlockList.append(pts[i][0])
                    self.yBlockList.append(pts[i][1])
                # Side 4
                for i in range(end-nBWidth,start-1,-nBWidth-1):
                    self.xBlockList.append(pts[i][0])
                    self.yBlockList.append(pts[i][1])
                self.xBlockList.append(None)
                self.yBlockList.append(None)

        tend = default_timer()
        if (self.silent == 0): print('Getting formatted blocklist time=%f sec' % (tend-tstart))
        return [self.xBlockList,self.yBlockList]

    def showValues(self,icell):
        if (self.silent == 0): print('=======================================================')
        if (self.silent == 0): print('icell= %d; x=%e; y=%e' % (icell,self.getCenterPoints()[icell,0],self.getCenterPoints()[icell,1]))
        for varname in self.getVarnames():
            exec("if (self.silent == 0): print('{} = ', self.{}[icell])".format(varname,varname))


    def getIcellByPoint(self,x,y):
         radii2 = (self.getCenterPoints()[:,0]-x)**2 + (self.getCenterPoints()[:,1]-y)**2
         icell=radii2.argmin()
         return icell


    def getPoints(self):
        try:
            self.data
        except AttributeError:
            self.getData()
        try:
            self.points
        except AttributeError:
            vtk_points=self.data.GetPoints().GetData()
            self.points=ah.vtk2array(vtk_points)
        return self.points


    def mirror(self):
        """
        Called when mirrorPlane != None
        The reflection plane is labeled as follows: From the vtkfiles documentation:
        ReflectionPlane {
        USE_X_MIN = 0, USE_Y_MIN = 1, USE_Z_MIN = 2, USE_X_MAX = 3,
        USE_Y_MAX = 4, USE_Z_MAX = 5, USE_X = 6, USE_Y = 7,
        USE_Z = 8
        }
        """

        vr=v.vtkReflectionFilter()
        vr.SetInputData(self.data)
        vr.SetPlane(self.mirrorPlane)
        self.data=vr.GetOutput()
        vr.Update()
        #self.data.Update()
        self.ncells = self.data.GetNumberOfCells()


    def reflectVar(self,var):
        if self.mirrorPlane == 0:
            CC=self.getCenterPoints()
            im = CC[:,0] < 0
            var[im] = -var[im]
        else:
            if (self.silent == 0): print('reflection of this plane not yet implemented, sorry')

        return var


    def getNdim(self):
        smalldouble = 1e-10
        self.ndim = 3
        if abs(self.data.GetBounds()[1] - self.data.GetBounds()[0]) <= smalldouble:
            self.ndim=self.ndim - 1
        if abs(self.data.GetBounds()[3] - self.data.GetBounds()[2]) <= smalldouble:
            self.ndim=self.ndim - 1
        if abs(self.data.GetBounds()[5] - self.data.GetBounds()[4]) <= smalldouble:
            self.ndim=self.ndim - 1
        return self.ndim


    def getAll(self):
        t0 = default_timer()
        self.getData()
        tdata = default_timer()
        if (self.silent == 0): print('Reading data time= %f sec' % (tdata-t0))

        if self.mirrorPlane is not None:
            if (self.silent == 0): print('========== Mirror about plane ',self.mirrorPlane,' ... ============')
            self.mirror()


        if (self.silent == 0): print('========== Initializing ... ============')

        self.getVars()
        tvars = default_timer()
        if (self.silent == 0): print('Getting vars time= %f sec' % (tvars-tdata))

        self.getPoints()
        tpoints = default_timer()
        if (self.silent == 0): print('Getting points time= %f sec' % (tpoints-tvars))

        self.getTime()

        tend = default_timer()
        if (self.silent == 0): print('========== Finished loading %d cells in %f sec, have a nice day! ===========' % (self.ncells, (tend-t0) ))

#=============================================================================
class loadvti(load_vtkfile):
    """Loader class for vti data"""
    def __init__(self,offset,get=1,file='data',mirrorPlane=None,silent=0,
                 rotateX=0,rotateY=0,rotateZ=0,
                 scaleX=1,scaleY=1,scaleZ=1):
        super(loadvti, self).__init__(offset=offset, get=get, file=file, type='vti', mirrorPlane=mirrorPlane,
                                      silent=silent, rotateX=rotateX, rotateY=rotateY, rotateZ=rotateZ, scaleX=scaleX,
                                      scaleY=scaleY, scaleZ=scaleZ)

    def getNdim(self):
        self.ndim = 3
        if self.data.GetExtent()[1] - self.data.GetExtent()[0] == 0.:
            self.ndim=self.ndim - 1
        if self.data.GetExtent()[3] - self.data.GetExtent()[2] == 0.:
            self.ndim=self.ndim - 1
        if self.data.GetExtent()[5] - self.data.GetExtent()[4] == 0.:
            self.ndim=self.ndim - 1
        return self.ndim

    def getX(self):
        self.dx=self.data.GetSpacing()[0]
        self.nx=self.data.GetExtent()[1]-self.data.GetExtent()[0]
        self.x=self.data.GetOrigin()[0] + np.arange(self.data.GetExtent()[0],self.data.GetExtent()[1])*self.dx+self.dx/2.

    def getY(self):
        self.dy=self.data.GetSpacing()[1]
        self.ny=self.data.GetExtent()[3]-self.data.GetExtent()[2]
        self.y=self.data.GetOrigin()[1] + np.arange(self.data.GetExtent()[2],self.data.GetExtent()[3])*self.dy+self.dy/2.

    def getZ(self):
        self.dz=self.data.GetSpacing()[2]
        self.nz=self.data.GetExtent()[5]-self.data.GetExtent()[4]
        self.z=self.data.GetOrigin()[2] + np.arange(self.data.GetExtent()[4],self.data.GetExtent()[5])*self.dz+self.dz/2.


    def getVars(self):
        nvars= self.data.GetCellData().GetNumberOfArrays()
        for i in range(nvars):
            varname = self.data.GetCellData().GetArrayName(i)
            if (self.silent == 0): print("Assigning variable:", varname)
            #            vtk_values = self.pointdata.GetPointData().GetArray(varname)
            vtk_values = self.data.GetCellData().GetArray(varname)
            if self.getNdim() == 1:
                exec("self.%s = ah.vtk2array(vtk_values)[0:self.ncells].reshape((self.nx),order='F')" % (varname))
            if self.getNdim() == 2:
                exec("self.%s = ah.vtk2array(vtk_values)[0:self.ncells].reshape((self.nx,self.ny),order='F')" % (varname))
            if self.getNdim() == 3:
                exec("self.%s = ah.vtk2array(vtk_values)[0:self.ncells].reshape((self.nx,self.ny,self.nz),order='F')" % (varname))


    def getAll(self):
        t0 = default_timer()
        self.getData()
        tdata = default_timer()
        if (self.silent == 0): print('Reading data time= %f sec' % (tdata-t0))

        self.getNdim()
        self.getX()
        self.getY()
        self.getZ()

        if self.mirrorPlane is not None:
            if (self.silent == 0): print('========== Mirror about plane ',self.mirrorPlane,' ... ============')
            self.mirror()


        if (self.silent == 0): print('========== Initializing ... ============')


        self.getVars()
        tvars = default_timer()
        if (self.silent == 0): print('Getting vars time= %f sec' % (tvars-tdata))

        self.getTime()

        tend = default_timer()
        if (self.silent == 0): print('========== Finished loading %d cells in %f sec, have a nice day! ===========' % (self.ncells, (tend-t0) ))

#=============================================================================
class loadcsv():
    """Load 1D comma separated list from a cut"""
    def __init__(self,offset,get=1,file='data',dir=1,coord=0.,silent=0):
        self.silent = silent
        self.offset=offset
        self.filenameout = file
        self.coord = coord
        self.dir   = dir
        self.isLoaded = False
        self.makefilename()

        if (self.silent == 0): print('========================================')
        if (self.silent == 0): print('loading file %s' % (self.filename))

        if get is not None:
            self.getAll()


    def makefilename(self):
        number='%+.2e' % self.coord
        if self.coord >= 0. and abs(self.coord) >1.:
            self.filename = self.filenameout+'_d'+str(self.dir)+'_x+'+'0.'+number[1]+number[3]+'D+'+str(int(np.log10(abs(self.coord)))+1).zfill(2)+'_n'+str(self.offset).zfill(4)+'.csv'
        if self.coord < 0. and abs(self.coord) >1.:
            self.filename = self.filenameout+'_d'+str(self.dir)+'_x-'+'0.'+number[1]+number[3]+'D+'+str(int(np.log10(abs(self.coord)))+1).zfill(2)+'_n'+str(self.offset).zfill(4)+'.csv'
        if self.coord > 0. and abs(self.coord) <=1.:
            self.filename = self.filenameout+'_d'+str(self.dir)+'_x+'+'0.'+number[1]+number[3]+'D+'+str(int(np.log10(abs(self.coord)))).zfill(2)+'_n'+str(self.offset).zfill(4)+'.csv'
        if self.coord < 0. and abs(self.coord) <=1.:
            self.filename = self.filenameout+'_d'+str(self.dir)+'_x-'+'0.'+number[1]+number[3]+'D+'+str(int(np.log10(abs(self.coord)))).zfill(2)+'_n'+str(self.offset).zfill(4)+'.csv'
        if self.coord == 0.:
            self.filename = self.filenameout+'_d'+str(self.dir)+'_x+'+'0.'+number[1]+number[3]+'D+'+str(0).zfill(2)+'_n'+str(self.offset).zfill(4)+'.csv'


    def getHeader(self):
        self.file.seek(0)
        self.header=self.file.readline().split(',')
        self.header[-1]=self.header[-1].rstrip('\n')

        return self.header


    def getAll(self):
        t0 = default_timer()
        self.file = open(self.filename, 'r')
        if (self.silent == 0): print('========== Initializing ... ============')
        self.getHeader()
        self.getData()
        self.file.close()
        tdata = default_timer()
        if (self.silent == 0): print('Reading data time= %f sec' % (tdata-t0))


    def getData(self):
        i=0
        for var in self.header:
            exec("self.%s = pylabload(self.filename,usecols=[i],skiprows=1, delimiter=',')" % (var))
            i=i+1
        exec("self.ncells = len(self.%s)" % (self.getHeader()[0]))

#=============================================================================
class particles():
    """Load binary particles data"""

    def __init__(self,offset,file='data',components=3):
        self.offset=offset
        self.filenameout = file
        self.isLoaded = False
        self.components = components

        self.data=[]
        self.mynparticles = 0

        self.makefilename()
        self.openfile()
        self.read()
        self.closefile()


    def makefilename(self):
        self.filename = self.filenameout+'_particles'+str(self.offset).zfill(4)+'.datfiles'


    def openfile(self):
        self.file = open(self.filename,'rb')
        self.isLoaded = True


    def closefile(self):
        self.file.close()
        self.isLoaded = False


    def readheader(self):
        self.file.seek(0)
        (self.nparticles,) = struct.unpack('i',self.file.read(4))
        (self.itparticles,) = struct.unpack('i',self.file.read(4))
        (self.npayload,) = struct.unpack('i',self.file.read(4))


    def read_next_particle(self):
        x       = np.empty(self.components)
        u       = np.empty(self.components)
        payload = np.empty(self.npayload)

        (index,)  = struct.unpack('i',self.file.read(4))

        (ifollow,) = struct.unpack('i',self.file.read(4))
        if (ifollow == -1):
            follow = True
        else:
            follow = False

        (q,)      = struct.unpack('d',self.file.read(8))
        (m,)      = struct.unpack('d',self.file.read(8))
        (t,)      = struct.unpack('d',self.file.read(8))
        (dt,)     = struct.unpack('d',self.file.read(8))

        for icomp in range(self.components):
            (x[icomp],) = struct.unpack('d',self.file.read(8))

        for icomp in range(self.components):
            (u[icomp],) = struct.unpack('d',self.file.read(8))

        for ipayload in range(self.npayload):
            (payload[ipayload],) = struct.unpack('d',self.file.read(8))
        
        self.data.append({'index':index, 'q':q, 
                          'follow':follow, 'm':m, 't':t, 
                          'dt':dt, 'x':x, 'u':u, 'payload':payload})
        
        self.mynparticles= self.mynparticles + 1


    def read(self):
        self.readheader()
        while (self.mynparticles < self.nparticles):
            self.read_next_particle()


    def particle(self,index):
        for particle in self.data:
            if particle['index'] == index :
                return particle
        return False


    def get_sorted(self,i):
        try:
            self.sorted
        except AttributeError:
            self.index=np.array(self.nparticles)
            index=[]
            for particle in self.data:
                index.append(particle['index'])
            index=np.array(index)
            self.sorted=index.argsort()
        return self.sorted[i-1]

#=============================================================================
class ensemble():
    """Load particle ensemble .csv files"""

    def __init__(self,offset,file='data0000',npayload=1,components=3,delimiter=','):
        self.offset=offset
        self.filenameout = file
        self.isLoaded = False
        self.components = components
        self.npayload = npayload
        self.delimiter = delimiter

        self.data=[]

        self.makefilename()
        self.read()


    def makefilename(self):
        self.filename = self.filenameout+'_ensemble'+str(self.offset).zfill(6)+'.csv'


    def read(self):
        x = pylabload(self.filename, comments='#', delimiter=self.delimiter)

        t_         = 0
        dt_        = 1
        x1_        = 2
        u1_        = x1_+self.components
        payload_   = u1_+self.components
        ipe_       = payload_+self.npayload
        iteration_ = ipe_+1
        index_     = iteration_+1

        if x.shape[1] == index_+1 : 

            for particle in x:
                self.data.append({'t':particle[t_],  
                                  'dt':particle[dt_], 'x':np.array(particle[x1_:x1_+self.components]), 
                                  'u':np.array(particle[u1_:u1_+self.components]), 
                                  'payload':np.array(particle[payload_:payload_+self.npayload]),
                                  'ipe':particle[ipe_].astype(np.int), 'iteration':particle[iteration_].astype(np.int), 'index':particle[index_].astype(np.int)})

        else:
            print('File inconsistent, assuming iteration counter overflow')
            for particle in x:
                self.data.append({'t':particle[t_],  
                                  'dt':particle[dt_], 'x':np.array(particle[x1_:x1_+self.components]), 
                                  'u':np.array(particle[u1_:u1_+self.components]), 
                                  'payload':np.array(particle[payload_:payload_+self.npayload]),
                                  'garbage':particle[ipe_].astype(np.int), 'index':particle[index_-1].astype(np.int)})
            

        self.isLoaded = True


    def particle(self,index):
        for particle in self.data:
            if particle['index'] == index :
                return particle
        return False
#=============================================================================
class postrad():
    """Load binary postrad data"""

    def __init__(self,offset,file='data',headersize=256,nvars=8):
        self.offset=offset
        self.filenameout = file
        self.headersize = headersize
        self.nvars =8

        self.indices={'r':0,'theta':1,'phi':2,'rho':0,'vr':1,'vtheta':2,'vphi':3,'p':4,'br':5,'btheta':6,'bphi':7}

        self.data=-1
        self.grid=-1
        self.header={}

        self.__makefilenames()
        self.__openfiles()
        self.__read()
        self.__closefiles()
        
    def __makefilenames(self):
        self.filename_vars = self.filenameout+str(self.offset).zfill(4)+'.blk'
        self.filename_grid = self.filenameout+'_grid.blk'

    def __openfiles(self):
        self.file_vars = open(self.filename_vars,'rb')
        self.file_grid = open(self.filename_grid,'rb')

    def __closefiles(self):
        self.file_vars.close()
        self.file_grid.close()

    def __read(self):
        self.__read_header()
        self.__set_dtype()
        self.__read_grid()
        self.__read_data()
        
    def __read_header(self):
        file_vars = self.file_vars
        file_vars.seek(0)
        nr = struct.unpack('i',file_vars.read(4)) [0]
        ntheta = struct.unpack('i',file_vars.read(4)) [0]
        nphi = struct.unpack('i',file_vars.read(4)) [0]
        t = struct.unpack('d',file_vars.read(8)) [0]
        a = struct.unpack('d',file_vars.read(8)) [0]
        r_in = struct.unpack('d',file_vars.read(8)) [0]
        theta_in = struct.unpack('d',file_vars.read(8)) [0]
        phi_in = struct.unpack('d',file_vars.read(8)) [0]
        r_out = struct.unpack('d',file_vars.read(8)) [0]
        theta_out = struct.unpack('d',file_vars.read(8)) [0]
        phi_out = struct.unpack('d',file_vars.read(8)) [0]
        gammahat = struct.unpack('d',file_vars.read(8)) [0]
        #    hslope  = struct.unpack('d',file_vars.read(8)) [0]
        hslope  = 0.25 # was incorrectly set in codeComparison files
        imetric = struct.unpack('i',file_vars.read(4)) [0]
        dummy = struct.unpack('i',file_vars.read(4)) [0]
        dim = struct.unpack('i',file_vars.read(4)) [0]
        isp = struct.unpack('i',self.file_vars.read(4)) [0]
        if (isp == -1):
            single_precision = True
        else:
            single_precision = False
        
        self.header = {'nr':nr, 'ntheta':ntheta, 'nphi':nphi,'t':t,'a':a,'r_in':r_in,'theta_in':theta_in,'phi_in':phi_in,'r_out':r_out,'theta_out':theta_out,'phi_out':phi_out,'gammahat':gammahat,'hslope':hslope,'imetric':imetric,'dim':dim,'single_precision':single_precision}
    
    def __set_dtype(self):
        if (self.header['single_precision']):
            self.__dtype_str = 'f'
        else:
            self.__dtype_str = 'd'

    def __read_data(self):
        file_vars = self.file_vars
        header = self.header
        file_vars.seek(self.headersize)
            
        len_bytes = header['nphi']*header['ntheta']*header['nr']*self.nvars
        data = np.array(struct.unpack(self.__dtype_str*len_bytes,file_vars.read())).reshape(header['nphi'],header['ntheta'],header['nr'],self.nvars)
        
        self.data = data

    def __read_grid(self):
        file_grid = self.file_grid
        header = self.header
        file_grid.seek(0)
        len_bytes = header['nphi']*header['ntheta']*header['nr']*3
        data = np.array(struct.unpack(self.__dtype_str*len_bytes,file_grid.read())).reshape(header['nphi'],header['ntheta'],header['nr'],3)
        
        self.grid = data
