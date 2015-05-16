import vtk as v
from tvtk import array_handler as ah
import matplotlib.pyplot as plt
import numpy as np
import fileinput
import sys
import os

filenamestub = 'rk54_cfl08_mp5'
offset = 10
resolutions = np.array([32,64,128,256,1024,2048])

def makefilename(filenamestub,offset,resolution):
    return ''.join([filenamestub,str(resolution),'_',str(offset).zfill(4),'.vti'])

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
    rho=ah.vtk2array(data.GetCellData().GetArray('rho'))

    return rho

def plot(data):
    x=getX(data)
    rho=getRho(data)

    plt.plot(x,rho,'k+-')


def getError(resolutions,norm,offset=10,filenamestub=filenamestub):

    error=[]
    for res in resolutions:
        solutionDataFilename = makefilename(''.join([filenamestub,'_']),offset,res)
        solutionData = load(solutionDataFilename)
        rho = getRho(solutionData)

        referenceDataFilename = makefilename(''.join([filenamestub,'_']),0,res)
        referenceData = load(referenceDataFilename)
        rho_0 = getRho(referenceData)

        error.append(norm(rho,rho_0))

    error = np.array(error)

    return error


def plotError(resolutions,error):
    
    plt.loglog(resolutions,error)

def l1(solution,reference):

    if len(solution) < len(reference):
        print 'solution datapoints dont match reference shape'
        return

    lref= len(reference)

    l1 = 1./lref*np.sum(np.abs(solution[0:lref]-reference[0:lref]))

    return l1

def linfty(solution,reference):

    if len(solution) < len(reference):
        print 'solution datapoints dont match reference shape'
        return

    lref= len(reference)

    linfty = np.max(np.abs(solution[0:lref]-reference[0:lref]))

    return linfty


def run(resolutions,filenamestub=filenamestub):
    for res in resolutions:
        parfile = ''.join(['amrvac_',str(res),'.par'])
        command = ''.join(['nohup ./amrvac -i ',parfile,' >> ',filenamestub,'_',str(res),'.out','&'])
        os.system(''.join(['cp amrvac.par.template ', parfile]))
        replaceAll(parfile,'RES',str(res))
        replaceAll(parfile,'FILE',filenamestub)
        os.system(command)

def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)
