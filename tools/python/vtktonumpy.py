from numpy import *
from vtk import *
import types
from vtk.util import vtkConstants

# Useful constants for VTK arrays.
VTK_ID_TYPE_SIZE = vtkIdTypeArray().GetDataTypeSize()
if VTK_ID_TYPE_SIZE == 4:
    ID_TYPE_CODE = int32
elif VTK_ID_TYPE_SIZE == 8:
    ID_TYPE_CODE = int64

VTK_LONG_TYPE_SIZE = vtkLongArray().GetDataTypeSize()
if VTK_LONG_TYPE_SIZE == 4:
    LONG_TYPE_CODE = int32S
    ULONG_TYPE_CODE = uint32
elif VTK_LONG_TYPE_SIZE == 8:
    LONG_TYPE_CODE = int64
    ULONG_TYPE_CODE = uint64

BASE_REFERENCE_COUNT = vtkObject().GetReferenceCount()

#def create_points(array):
def array2vtkPoints(array):
    """Create vtkPoints from double array"""
    vtk_points = vtkPoints()
    double_array = vtkDoubleArray()
    double_array.SetVoidArray(array, len(array), 1)
    double_array.SetNumberOfComponents(3)
    vtk_points.SetData(double_array)
    return vtk_points

def array2vtk(num_array, vtk_array=None):
    """Converts a real numpy Array (or a Python list) to a VTK array
    object.

    This function only works for real arrays.  Complex arrays are NOT
    handled.  It also works for multi-component arrays.  However, only
    1, and 2 dimensional arrays are supported.  This function is very
    efficient, so large arrays should not be a problem.

    Even in cases when no copy of the numpy array data is performed,
    a reference to the array is cached.  The passed array can
    therefore be deleted safely in all circumstances.

    Parameters
    ----------

    - num_array : numpy array or Python list/tuple

      The input array must be 1 or 2D.  A copy of the numeric array
      data passed is made in the following circumstances:

       1. A Python list/tuple was passed.
       2. A non-contiguous numpy array was passed.
       3. A `vtkBitArray` instance was passed as the second argument.
       4. The types of the `vtk_array` and the `num_array` are not
          equivalent to each other.  For example if one is an integer
          array and the other a float.


    - vtk_array : `vtkDataArray` (default: `None`)

      If an optional `vtkDataArray` instance, is passed as an argument
      then a new array is not created and returned.  The passed array
      is itself returned.

    """

    z = asarray(num_array)

    shape = z.shape
    assert len(shape) < 3, \
           "Only arrays of dimensionality 2 or lower are allowed!"
    assert not issubdtype(z.dtype, complex), \
           "Complex numpy arrays cannot be converted to vtk arrays."\
           "Use real() or imag() to get a component of the array before"\
           " passing it to "

    # First create an array of the right type by using the typecode.
    # Bit arrays need special casing.
    bit_array = False
    if vtk_array is None:
        vtk_typecode = get_vtk_array_type(z.dtype)
        result_array = create_vtk_array(vtk_typecode)
    elif vtk_array.GetDataType() == vtkConstants.VTK_BIT:
        vtk_typecode = vtkConstants.VTK_CHAR
        result_array = create_vtk_array(vtkConstants.VTK_CHAR)
        bit_array = True
    else:
        vtk_typecode = vtk_array.GetDataType()
        result_array = vtk_array

    # Find the shape and set number of components.
    if len(shape) == 1:
        result_array.SetNumberOfComponents(1)
    else:
        result_array.SetNumberOfComponents(shape[1])

    result_array.SetNumberOfTuples(shape[0])

    # Ravel the array appropriately.
    arr_dtype = get_numeric_array_type(vtk_typecode)
    if issubdtype(z.dtype, arr_dtype):
        z_flat = ravel(z)
    else:
        z_flat = ravel(z).astype(arr_dtype)

    # Point the VTK array to the numpy data.  The last argument (1)
    # tells the array not to deallocate.
    result_array.SetVoidArray(getbuffer(z_flat), len(z_flat), 1)

    if bit_array:
        # Handle bit arrays -- they have to be copied.  Note that bit
        # arrays are used ONLY when the user has passed one as an
        # argument to this function.
        vtk_array.SetNumberOfTuples(result_array.GetNumberOfTuples())
        vtk_array.SetNumberOfComponents(result_array.GetNumberOfComponents())
        for i in range(result_array.GetNumberOfComponents()):
            vtk_array.CopyComponent(i, result_array, i)
        result_array = vtk_array

    return result_array

def create_vtk_array(vtk_arr_type):
    """Internal function used to create a VTK data array from another
    VTK array given the VTK array type.
    """
    tmp = vtkDataArray.CreateDataArray(vtk_arr_type)
    # CreateDataArray sets the refcount to 3 and this causes a severe
    # memory leak.
    tmp.SetReferenceCount(BASE_REFERENCE_COUNT)
    return tmp

def get_vtk_to_numeric_typemap():
    """Returns the VTK array type to numpy array type mapping."""
    _vtk_arr = {vtkConstants.VTK_BIT:bool,
                vtkConstants.VTK_CHAR:int8,
                vtkConstants.VTK_UNSIGNED_CHAR:uint8,
                vtkConstants.VTK_SHORT:int16,
                vtkConstants.VTK_UNSIGNED_SHORT:uint16,
                vtkConstants.VTK_INT:int32,
                vtkConstants.VTK_UNSIGNED_INT:uint32,
                vtkConstants.VTK_LONG:LONG_TYPE_CODE,
                vtkConstants.VTK_UNSIGNED_LONG:ULONG_TYPE_CODE,
                vtkConstants.VTK_ID_TYPE:ID_TYPE_CODE,
                vtkConstants.VTK_FLOAT:float32,
                vtkConstants.VTK_DOUBLE:float64}
    return _vtk_arr

def get_numeric_array_type(vtk_array_type):
    """Returns a numpy array typecode given a VTK array type."""
    return get_vtk_to_numeric_typemap()[vtk_array_type]


def create_cells(array):
    """Create a vtkCellArray from long array"""
    vtk_cells = vtkCellArray()
    vtk_id_array = vtkIdTypeArray()
    vtk_id_array.SetVoidArray(array, len(array), 1)
    vtk_cells.SetCells(len(array)/4, vtk_id_array)
    return vtk_cells

def vtk2array(vtk_array):
    at = vtk_array.GetDataType()
    if at == 10:
        #vtkFloatArray
        pt=float32
    if at == 11:
        #vtkDoubleArray
        pt=float64
    elif at == 12:
        #vtkIdTypeArray
        pt=int
    ndim=vtk_array.GetNumberOfComponents()
    if ndim > 1:
        r = empty((vtk_array.GetSize()/ndim,ndim),dtype=pt)
    else: 
        r = empty(vtk_array.GetSize(),dtype=pt)
    vtk_array.ExportToVoidPointer(r)
    return r

######################################################################
# Array conversion functions.
######################################################################
def get_vtk_array_type(numeric_array_type):
    """Returns a VTK typecode given a numpy array."""
    # This is a Mapping from numpy array types to VTK array types.
    _arr_vtk = {dtype(character):vtkConstants.VTK_UNSIGNED_CHAR,
                dtype(uint8):vtkConstants.VTK_UNSIGNED_CHAR,
                dtype(uint16):vtkConstants.VTK_UNSIGNED_SHORT,
                dtype(int8):vtkConstants.VTK_CHAR,
                dtype(int16):vtkConstants.VTK_SHORT,
                dtype(int32):vtkConstants.VTK_INT,
                dtype(uint32):vtkConstants.VTK_UNSIGNED_INT,
                dtype(float32):vtkConstants.VTK_FLOAT,
                dtype(float64):vtkConstants.VTK_DOUBLE,
                dtype(complex64):vtkConstants.VTK_FLOAT,
                dtype(complex128):vtkConstants.VTK_DOUBLE,
                }
    _extra = {dtype(ID_TYPE_CODE):vtkConstants.VTK_ID_TYPE,
              dtype(ULONG_TYPE_CODE):vtkConstants.VTK_UNSIGNED_LONG,
              dtype(LONG_TYPE_CODE):vtkConstants.VTK_LONG,
             }
    for t in _extra:
        if t not in _arr_vtk:
            _arr_vtk[t] = _extra[t]

    try:
        return _arr_vtk[numeric_array_type]
    except KeyError:
        for key in _arr_vtk:
            if issubdtype(numeric_array_type, key):
                return _arr_vtk[key]
    raise TypeError, "Couldn't translate array's type to VTK"
