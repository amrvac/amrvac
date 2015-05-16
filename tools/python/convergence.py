import vtk as v
from tvtk import array_handler as ah
import matplotlib.pyplot as plt
import numpy as np
import fileinput
import sys
import os


folderstub='rk54_fd_mp5'
filenamestub = 'data'
offset = 1
resolutions = np.array([64,128,256,512])

def makefoldername(folderstub,resolution):
    return ''.join([folderstub,'_',str(resolution)])

def makefilename(folderstub,resolution,filenamestub,offset):
    return ''.join([makefoldername(folderstub,resolution),'/',filenamestub,str(offset).zfill(4),'.vti'])


def load(file):
    datareader = v.vtkXMLImageDataReader()
    datareader.SetFileName(file)
    datareader.Update()
    data = datareader.GetOutput()

    return data


def getX(data):
    dx=data.GetSpacing()[0]
    nx=data.GetExtent()[1]
    x=np.arange(0,1,dx)+dx/2.

    return x

def getRho(data):
    rho=(ah.vtk2array(data.GetCellData().GetArray('rho'))).astype(np.float64)
    return rho

def getVar(data,var='rho'):
    var=(ah.vtk2array(data.GetCellData().GetArray(var))).astype(np.float64)
    return var


def plot(data):
    x=getX(data)
    rho=getRho(data)

    plt.plot(x,rho,'k+-')


def getError(resolutions,norm,var=None,offset=1,folderstub=folderstub,filenamestub=filenamestub):

    error=[]
    for res in resolutions:
        solutionDataFilename = makefilename(folderstub,res,filenamestub,offset)
        solutionData = load(solutionDataFilename)
        referenceDataFilename = makefilename(folderstub,res,filenamestub,0)
        referenceData = load(referenceDataFilename)

        if var != None:
            sol = getVar(solutionData,var=var)
            sol_0 = getVar(referenceData,var=var)
        else:
            sol = sol_0 = []
#            for ivar in range(referenceData.GetCellData().GetNumberOfArrays()):
            for ivar in range(9):
                varname = referenceData.GetCellData().GetArrayName(ivar)
                print 'Resolution: ',res, ', appending to solution vector: ', varname
                sol = np.append(sol,getVar(solutionData,var=varname))
                sol_0 = np.append(sol_0,getVar(referenceData,var=varname))
        error.append(norm(sol,sol_0))

    error = np.array(error)

    return error


def plotError(resolutions,error):
    
    plt.loglog(resolutions,error)


def l1(solution,reference):
    '''Calculates the L1 norm between two solutions'''

    if np.shape(solution) != np.shape(reference):
        print 'solution datapoints dont match reference shape'
        return

    lref= len(reference.flatten())

    l1 = 1./lref*np.sum(np.abs(solution-reference))

    return l1


def linfty(solution,reference):
    '''Calculates the Linfinity norm between two solutions'''

    if np.shape(solution) != np.shape(reference):
        print 'solution datapoints dont match reference shape'
        return

    lref= len(reference.flatten())

    linfty = np.max(np.abs(solution-reference))

    return linfty


def run(resolutions,folderstub=folderstub,filenamestub=filenamestub):
    for res in resolutions:
        parfile = ''.join([makefoldername(folderstub,res),'/amrvac.par'])
        jobfile = ''.join(['job_',makefoldername(folderstub,res),'.sh'])
#        command = ''.join(['nohup mpirun -np 2 ./amrvac -i ',parfile,' >> ',filenamestub,'_',str(res),'.out'])
        os.system(''.join(['mkdir ',makefoldername(folderstub,res)]))

        os.system(''.join(['cp amrvac.par.template ', parfile]))
        replaceAll(parfile,'XRES',str(res))
        replaceAll(parfile,'YRES',str(int(res/2.)))
        replaceAll(parfile,'FILE',filenamestub)
        replaceAll(parfile,'FOLDER',makefoldername(folderstub,res))

        os.system(''.join(['cp job.sh.template ', jobfile]))
        replaceAll(jobfile,'NCPU',str(min([int(res/8.*res/8.*res/8./4.),120])))
        replaceAll(jobfile,'YRES',str(res))
        replaceAll(jobfile,'FILE',filenamestub)
        replaceAll(jobfile,'FOLDER',makefoldername(folderstub,res))


#        os.system(command)

def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)
