# -------------------------------------------------
import numpy as np
from scipy.interpolate import griddata
import vtk as v
#from tvtk import array_handler as ah
import numpy_support as ah

def vtu_values(offset,varname,filenameout='data',attribute_mode='cell'):

    filename=''.join([filenameout,repr(offset).zfill(4),'.vtu'])

    datareader = v.vtkXMLUnstructuredGridReader()
    datareader.SetFileName(filename)
    datareader.Update()
    data = datareader.GetOutput()
    
    if attribute_mode == 'cell':
        vtk_values = data.GetCellData().GetArray(varname)
    elif attribute_mode == 'point':
        vtk_values = data.GetPointData().GetArray(varname)
    else:
        print "attribute_mode is either 'cell' or 'point'"
        
# this is convenient to convert vtkarrays to numpy arrays
    value = ah.vtk2array(vtk_values)
    return value
# ------------------------------------------------------
def vtu_points(offset,filenameout='data',):

    filename=''.join([filenameout,repr(offset).zfill(4),'.vtu'])

    datareader = v.vtkXMLUnstructuredGridReader()
    datareader.SetFileName(filename)
    datareader.Update()
    data = datareader.GetOutput()
    
    vtk_points=data.GetPoints().GetData()
    points=ah.vtk2array(vtk_points)
    return points
# ------------------------------------------------------
def vtu_regrid(offset,varname,nregrid,filenameout='data',attribute_mode='point'):        

    if attribute_mode != 'point':
        print 'regrid needs point data'
        return None
    
    var = vtu_values(offset,varname,filenameout=filenameout,attribute_mode=attribute_mode)
    points = vtu_points(offset,filenameout=filenameout)

    x=points[:,0]
    y=points[:,1]
    z=points[:,2]

    # reduce dimension of points:
    ndim = points.ndim

    # this is needed because griddata freaks out otherwise:
    smalldouble=1.e-8

    if ndim == 1:
        xi = np.linspace(x.min(),x.max(),nregrid[0])
        data = griddata(x, var, xi, method='linear')
        data_dict={varname: data, 'x': xi}

    if ndim == 2:
        tmp=(x.max()-x.min())*smalldouble
        xi = np.linspace(x.min()-tmp,x.max()+tmp,nregrid[0])
        tmp=(y.max()-y.min())*smalldouble
        yi = np.linspace(y.min()-tmp,y.max()+tmp,nregrid[1])
        data = griddata((x,y), var, (xi[None,:],yi[:,None]), method='linear')
        data_dict={varname: data, 'x': xi, 'y': yi}

    # Not tested (must be very slow): 
    if ndim == 3:
        tmp=(x.max()-x.min())*smalldouble
        xi = np.linspace(x.min()-tmp,x.max()+tmp,nregrid[0])
        tmp=(y.max()-y.min())*smalldouble
        yi = np.linspace(y.min()-tmp,y.max()+tmp,nregrid[1])
        tmp=(z.max()-z.min())*smalldouble
        zi = np.linspace(z.min()-tmp,z.max()+tmp,nregrid[2])
        data = griddata((x,y,z), var, (xi[None,:,:],yi[:,:,None],zi[:,:,None]), method='linear')        
        data_dict={varname: data, 'x': xi, 'y': yi, 'z': zi}
    
    return data_dict
# ------------------------------------------------------
