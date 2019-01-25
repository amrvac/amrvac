"""
Created 4 Dec 2018

Module to read in MPI-AMRVAC data (.dat files) for 1D, 2D and 3D. Currently implemented methods are:
== BASIC METHODS ==
- get_header        : Reads the header from an MPI-AMRVAC 2.0 snapshot.
- get_block_data    : Returns an array with block data from a snapshot.

==FOR UNIFORM GRID DATA==
- get_uniform_data  : Reads in the block data for a uniformely refined grid and returns
                      the raw data as a dictionary.

==FOR NONUNIFORM GRID DATA==
- interpolate_block_1d : Interpolates a single 1D block to the maximum refinement level over the grid.
- interpolate_block_2d : Interpolates a single 2D block to the maximum refinement level over the grid.
- interpolate_block_3d : Interpolates a single 3D block to the maximum refinement level over the grid.
- print_regrid_amount  : Prints the total amount of blocks in the grid + the amount needing regridding.
- get_amr_data         : Interpolates the blocks to the maximum refinement level, returns the raw data.
                         Method optimized to run on a single CPU.
- save_regridded_data  : Method to save the regridded raw data as Numpy .npy files which can be loaded in later.
- increment_progress_mp: Small method to track progress when using multiple cores.
- multiprocessing_init : Initialization method for the multiprocessing module. Defines global variables
                         that can be accessed across multiple cores, which is not the case for settings.py variables.
- get_amr_data_multiprocessing(dat): Similar to get_amr_data, except this method is optimized for usage
                                     with the multiprocessing module. Performs regridding of multiple blocks
                                     simultaneously using parallelization.
- methods with _unpack(): The parallelization routine uses multiprocessing.starmap(), which is only supported from
                          Python 3.3 onwards. The _unpack() routines are used for compatibility with Python 2,
                          these unpack the multiple arguments which are otherwise passed on to starmap().
                          Now map() is used instead.

@author: Jannis Teunissen (original)
@author: Niels Claes (modifications)
         Last edit: 07 January 2019

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
# from __future__ import unicode_literals

import struct
import numpy as np
import os, sys
import scipy.interpolate as interp
import print_tools
import settings
import multiprocessing

# Size of basic types (in bytes)
size_logical = 4
size_int = 4
size_double = 8
name_len = 16

# For un-aligned data, use '=' (for aligned data set to '')
align = '='


def get_header(dat): 
    """
    Reads header from MPI-AMRVAC 2.0 snapshot.
    :param dat: .dat file opened in binary mode.
    :return: Dictionary containing header data from snapshot.
    """

    dat.seek(0)
    h = {}

    fmt = align + 'i'
    [h['version']] = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))

    if h['version'] < 3:
        print("Unsupported .dat file version: ", h['version'])
        sys.exit()

    # Read scalar data at beginning of file
    fmt = align + 9 * 'i' + 'd'
    hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    [h['offset_tree'], h['offset_blocks'], h['nw'],
     h['ndir'], h['ndim'], h['levmax'], h['nleafs'], h['nparents'],
     h['it'], h['time']] = hdr

    # Read min/max coordinates
    fmt = align + h['ndim'] * 'd'
    h['xmin'] = np.array(struct.unpack(fmt, dat.read(struct.calcsize(fmt))))
    h['xmax'] = np.array(struct.unpack(fmt, dat.read(struct.calcsize(fmt))))

    # Read domain and block size (in number of cells)
    fmt = align + h['ndim'] * 'i'
    h['domain_nx'] = np.array(
        struct.unpack(fmt, dat.read(struct.calcsize(fmt))))
    h['block_nx'] = np.array(
        struct.unpack(fmt, dat.read(struct.calcsize(fmt))))

    # Read w_names
    w_names = []
    for i in range(h['nw']):
        fmt = align + name_len * 'c'
        hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
        w_names.append(b''.join(hdr).strip().decode())
    h['w_names'] = w_names

    # Read physics type
    fmt = align + name_len * 'c'
    hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    h['physics_type'] = b''.join(hdr).strip().decode()

    # Read number of physics-defined parameters
    fmt = align + 'i'
    [n_pars] = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))

    # First physics-parameter values are given, then their names
    fmt = align + n_pars * 'd'
    vals = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))

    fmt = align + n_pars * name_len * 'c'
    names = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    # Split and join the name strings (from one character array)
    names = [b''.join(names[i:i+name_len]).strip().decode()
             for i in range(0, len(names), name_len)]

    # Store the values corresponding to the names
    for val, name in zip(vals, names):
        h[name] = val
    return h


def get_block_data(dat):
    """
    Reads block data from an MPI-AMRVAC 2.0 snapshot.
    :param dat: .dat file opened in binary mode.
    :return: Dictionary containing block data.
    """

    dat.seek(0)
    h = get_header(dat)
    nw = h['nw']
    block_nx = np.array(h['block_nx'])
    domain_nx = np.array(h['domain_nx'])
    xmax = np.array(h['xmax'])
    xmin = np.array(h['xmin'])
    nleafs = h['nleafs']
    nparents = h['nparents']
    
    # Read tree info. Skip 'leaf' array
    dat.seek(h['offset_tree'] + (nleafs+nparents) * size_logical)

    # Read block levels
    fmt = align + nleafs * 'i'
    block_lvls = np.array(struct.unpack(fmt, dat.read(struct.calcsize(fmt))))

    # Read block indices
    fmt = align + nleafs * h['ndim'] * 'i'
    block_ixs = np.reshape(
        struct.unpack(fmt, dat.read(struct.calcsize(fmt))),
        [nleafs, h['ndim']])

    # Start reading data blocks
    dat.seek(h['offset_blocks'])

    blocks = []

    for i in range(nleafs):
        lvl = block_lvls[i]
        ix = block_ixs[i]

        # Read number of ghost cells
        fmt = align + h['ndim'] * 'i'
        gc_lo = np.array(struct.unpack(fmt, dat.read(struct.calcsize(fmt))))
        gc_hi = np.array(struct.unpack(fmt, dat.read(struct.calcsize(fmt))))

        # Read actual data
        block_shape = np.append(gc_lo + block_nx + gc_hi, nw)
        fmt = align + np.prod(block_shape) * 'd'
        d = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
        w = np.reshape(d, block_shape, order='F') # Fortran ordering

        b = {}
        b['lvl'] = lvl
        b['ix'] = ix
        b['w'] = w
        blocks.append(b)
        
    return blocks


def get_uniform_data(dat):
    """
    Reads block data from an MPI-AMRVAC 2.0 snapshot
    and returns the data as a dictionary. Assumes a uniformely refined grid.
    :param dat: .dat file opened in binary mode.
    :return: Dictionary containing grid data
    :raise IOError: If data is not uniformly refined.
    """
    h = get_header(dat)
    blocks = get_block_data(dat)

    # Check if grid is uniformly refined
    refined_nx = 2**(h['levmax']-1) * h['domain_nx']
    nleafs_uniform = np.prod(refined_nx/h['block_nx'])

    if h['nleafs'] == nleafs_uniform:
        print("Data is uniformely refined.")
        domain_shape = np.append(refined_nx, h['nw'])
        d = np.zeros(domain_shape, order='F')

        for b in blocks:
            i0 = (b['ix'] - 1) * h['block_nx']
            i1 = i0 + h['block_nx']
            if h['ndim'] == 1:
                d[i0[0]:i1[0], :] = b['w']
            elif h['ndim'] == 2:
                d[i0[0]:i1[0], i0[1]:i1[1], :] = b['w']
            elif h['ndim'] == 3:
                d[i0[0]:i1[0], i0[1]:i1[1], i0[2]:i1[2], :] = b['w']
        return d
    raise IOError('Data in .dat file is not uniformly refined.')


def interpolate_block_1d(b, hdr):
    """
    Interpolates a 1-dimensional block to a new 1D-array, depending
    on the current block level and the maximum refinement level of the snapshot.
    This information is contained in the header and block dictionary from get_block_data().
    :param b: Current block
    :param hdr: Header of the current snapshot.
    :return: Block interpolated to the new required 1D block dimension according to the difference
             in mesh refinement.
             Type is np.ndarray of shape (nx,nw) where nx is the number of gridpoints
             at maximum refinement level, and nw is the number of conservative variables
             contained in the header.
    """
    block_lvl = b['lvl']
    max_lvl   = hdr['levmax']
    # Check for multiprocessing method
    if block_lvl == max_lvl:
        increment_progress_mp()
        return b

    grid_diff = 2**(max_lvl - block_lvl)
    lvl1_block_width = hdr['block_nx']
    
    curr_x   = lvl1_block_width[0]
    regrid_x = lvl1_block_width[0] * grid_diff
    
    b_interpolated = np.zeros([regrid_x, hdr['nw']])
    
    for cons_var in range(0, hdr['nw']):
        block_spline = interp.interp1d(np.arange(curr_x), b[:, cons_var])
        block_result = block_spline(np.linspace(0, b[:, cons_var].size-1, regrid_x))
        
        b_interpolated[:, cons_var] = block_result

    if mp_activated.value:
        increment_progress_mp()
    
    return b_interpolated


def interpolate_block_2d(b, hdr):
    """
    Interpolates a 2-dimensional block to a new 2D-matrix, depending
    on the current block level and the maximum refinement level of the snapshot.
    This information is contained in the header and block dictionary from get_block_data().
    Regridding process:
        1) Get current block level and maximum grid level
        2) Calculate the grid difference, which is the multiplication factor for regridding.
           For a block 1 level behind max_lvl this means 1 additional refinement level,
           eg. regridding block_nx from 16x16 to 32x32. Grid difference in this case is 2.
           For a block 2 levels behind max_lvl, the block is regridded from 16x16 to 64x64.
           Now the grid difference is 4, etc.
        3) Regridding is performed for each conservative variable, using scipy.interpolate.griddata
        4) The new interpolated block is returned, the appropriate indices to fill the regridded matrix
           are calculated in get_amr_data() before calling this method. 
    :param b: Current block
    :param hdr: The header from the current snapshot.
    :return: Block interpolated to the new required 2D block dimension according to the difference in
             mesh refinement.
             Type is np.ndarray of shape (nx, ny, nw) where nx and ny are the number of gridpoints
             in each direction at maximum refinement level, and nw is the number of conservative variables
             contained in the header.                       
    """
    block_lvl = b['lvl']
    max_lvl   = hdr['levmax']
    # Check for multiprocessing method
    if block_lvl == max_lvl:
        increment_progress_mp()
        return b

    grid_diff = 2**(max_lvl - block_lvl)
    lvl1_block_width = hdr['block_nx']
    
    curr_x   = lvl1_block_width[0]
    curr_y   = lvl1_block_width[1]
    regrid_x = lvl1_block_width[0] * grid_diff
    regrid_y = lvl1_block_width[1] * grid_diff
    
    nb_elements = curr_x * curr_y
    b_interpolated = np.zeros([regrid_x, regrid_y, hdr['nw']])
    
    for cons_var in range(0, hdr['nw']):
        vals = np.reshape(b['w'][:, :, cons_var], nb_elements)
        pts  = np.array(  [[i, j] for i in np.linspace(0, 1, curr_x) for j in np.linspace(0, 1, curr_y)]  )
        
        grid_x, grid_y = np.mgrid[0:1:regrid_x*1j, 0:1:regrid_y*1j]
        #regrid mesh to required values, use linear interpolation
        grid_interpolated = interp.griddata(pts, vals, (grid_x, grid_y), method="linear")
        
        b_interpolated[:, :, cons_var] = grid_interpolated

    if mp_activated.value:
        increment_progress_mp()

    return b_interpolated


def interpolate_block_3d(b, hdr):
    """
    Interpolates a 3-dimensional block to a new 3D-matrix, depending
    on the current block level and the maximum refinement level of the snapshot.
    This information is contained in the header and block dictionary from get_block_data().
    :param b: Current block
    :param hdr: The header from the current snapshot.
    :return: Block interpolated to the new required 3D block dimension according to the difference in
             mesh refinement.
             Type is np.ndarray of shape (nx, ny, nz, nw) where nx, ny and nz are the number of gridpoints
             in each direction at maximum refinement level, and nw is the number of conservative variables
             contained in the header.
    """
    block_lvl = b['lvl']
    max_lvl   = hdr['levmax']
    # Check for multiprocessing method
    if block_lvl == max_lvl:
        increment_progress_mp()
        return b

    grid_diff = 2**(max_lvl - block_lvl)
    lvl1_block_width = hdr['block_nx']
    
    curr_x   = lvl1_block_width[0]
    curr_y   = lvl1_block_width[1]
    curr_z   = lvl1_block_width[2]
    regrid_x = lvl1_block_width[0] * grid_diff
    regrid_y = lvl1_block_width[1] * grid_diff
    regrid_z = lvl1_block_width[1] * grid_diff
    
    nb_elements = curr_x * curr_y * curr_z
    b_interpolated = np.zeros([regrid_x, regrid_y, regrid_z, hdr['nw']])
    
    for cons_var in range(0, hdr['nw']):
        vals = np.reshape(b['w'][:, :, :, cons_var], nb_elements)
        pts  = np.array(  [[i, j, k] for i in np.linspace(0, 1, curr_x)
                                     for j in np.linspace(0, 1, curr_y)
                                     for k in np.linspace(0, 1, curr_z)]   )
        
        grid_x, grid_y, grid_z = np.mgrid[0:1:regrid_x*1j, 0:1:regrid_y*1j, 0:1:regrid_z*1j]
        grid_interpolated = interp.griddata(pts, vals, (grid_x, grid_y, grid_z), method='linear')
        
        b_interpolated[:, :, :, cons_var] = grid_interpolated

    if mp_activated.value:
        increment_progress_mp()
    
    return b_interpolated


def print_regrid_amount(blocks, max_lvl):
    """
    Simple routine to count number of blocks in file and checks
    how many need regridding.
    :param blocks: List of blocks, output from get_block_data
    :param max_lvl: Maximum refinement level in the grid
    """
    block_regrid_needed = 0
    print("    Number of blocks in file:", len(blocks))
    for b in blocks:
        if not b['lvl'] == max_lvl:
            block_regrid_needed += 1
    print("    Number of blocks needing regridding: %s" % block_regrid_needed)

    
def get_amr_data(dat):
    """
    Returns a uniform grid in the case the mesh is not uniformely refined, hence
    when the method call to get_uniform_data() throws an IOError.
    This method calculates the maximum refinement level present in the grid, and regrids
    the entire mesh to this level. Blocks at a higher level than the maximum are refined using
    linear interpolation.
    :param dat: .dat file, opened in binary mode.
    :return: Dictionary containing grid data.
    :raise IOError: If number of dimensions in the header is not equal to 1, 2 or 3 for some reason.
    """
    h = get_header(dat)
    blocks = get_block_data(dat)
    
    refined_nx = 2**(h['levmax'] - 1) * h['domain_nx']
    #Perform regridding to finest level
    domain_shape = np.append(refined_nx, h['nw'])
    d = np.zeros(domain_shape, order='F')
    
    max_lvl    = h['levmax']
    #Get amount of blocks that need regridding
    print_regrid_amount(blocks, max_lvl)

    # No multiprocessing, so set mp_activated to False
    mp_bool = multiprocessing.Value("i", False)
    multiprocessing_init(0, mp_bool, len(blocks))
    
    counter = 0
    print_tools.progress(counter, len(blocks), status='-- iterating over blocks...')
    for b in blocks:
        
        block_lvl = b['lvl']
        block_idx = b['ix']
        
        grid_diff = 2**(max_lvl - block_lvl)
        
        max_idx = block_idx * grid_diff
        min_idx = max_idx - grid_diff

        if h['ndim'] == 1:
            idx0 = min_idx * h['block_nx']
            if block_lvl == max_lvl:
                idx1 = idx0 + h['block_nx']
                d[idx0[0]:idx1[0], :] = b['w']
            else:
                idx1 = idx0 + (h['block_nx'] * grid_diff)
                d[idx0[0]:idx1[0], :] = interpolate_block_1d(b, h)
                
        elif h['ndim'] == 2:
            idx0 = min_idx * h['block_nx']
            if block_lvl == max_lvl:
                #block is on finest level, return block
                idx1 = idx0 + h['block_nx']
                d[idx0[0]:idx1[0], idx0[1]:idx1[1], :] = b['w']
            else:
                #block is not on finest level, so interpolate
                idx1 = idx0 + (h['block_nx'] * grid_diff)
                d[idx0[0]:idx1[0], idx0[1]:idx1[1], :] = interpolate_block_2d(b, h)
                
        elif h['ndim'] == 3:
            idx0 = min_idx * h['block_nx']
            if block_lvl == max_lvl:
                idx1 = idx0 + h['block_nx']
                d[idx0[0]:idx1[0], idx0[1]:idx1[1], idx0[2]:idx1[2], :] = b['w']
            else:
                idx1 = idx0 + (h['block_nx'] * grid_diff)
                d[idx0[0]:idx1[0], idx0[1]:idx1[1], idx0[2]:idx1[2], :] = interpolate_block_3d(b, h)
        else:
            raise IOError("Unknown number of dimensions %s" % h['ndim'])
        
        counter += 1
        if counter % 10 == 0 or counter == len(blocks):
            print_tools.progress(counter, len(blocks), status='-- iterating over blocks...')
    print("\n")
            
    save_regridded_data(d)
    
    return d    

def save_regridded_data(regrid_data):
    """
    Saves the regridded data as a Numpy file.
    :param regrid_data: The regridded data, output from get_amr_data().
    """
    if settings.saveFiles:
        if not os.path.isdir("dat_files"):
            os.mkdir("dat_files")
        fn = 'dat_files/' + print_tools.trim_filename(settings.filename) + "_regridded_dat"
        np.save(fn, regrid_data)
        print("Regridded data saved to %s.npy" % fn)
    return


def increment_progress_mp():
    """
    Small method to increment progress when using multiple processors
    """
    progress.value += 1
    if progress.value % 10 == 0 or progress.value == total.value:
        print_tools.progress(progress.value, total.value, status='-- iterating over blocks...')
    return


def interpolate_block_1d_unpack(args):
    """
    1D block interpolation unpacking function for multiprocessing PY2 compatibility.
    :param args: (Iterable) Arguments for function call
    :return: Function call with unpacked arguments
    """
    return interpolate_block_1d(*args)


def interpolate_block_2d_unpack(args):
    """
    2D block interpolation unpacking function for multiprocessing PY2 compatibility.
    :param args: (Iterable) Arguments for function call
    :return: Function call with unpacked arguments
    """
    return interpolate_block_2d(*args)


def interpolate_block_3d_unpack(args):
    """
    3D block interpolation unpacking function for multiprocessing PY2 compatibility.
    :param args: (Iterable) Arguments for function call
    :return: Function call with unpacked arguments
    """
    return interpolate_block_3d(*args)


def multiprocessing_init(t, mp_bool, total_blocks):
    """
    Initialization method for multiprocessing progress tracking. Also create a
    boolean value as variables defined in other modules are read as their
    default values during multiprocessing, even when redefined.
    :param t: Multiprocessing.Value object at initialization of worker Pool
    :param mp_bool: (Boolean) True if multiprocessing is called, False otherwise
    :param total_blocks: (int) Total number of blocks, used to calculate progress
    """
    global progress, mp_activated, total
    progress = t
    mp_activated = mp_bool
    total = total_blocks


def get_amr_data_multiprocessing(dat):
    """
    Method to regrid entire mesh using the multiprocessing module.
    Same principle as the get_amr_data() method, except now each block that
    needs regridding is passed on to one of the multiple processors in use.
    :param dat: .dat file, opened in binary mode.
    :return: Dictionary containing grid data.
    """
    # Perform version check
    PY2 = sys.version_info[0] == 2

    h = get_header(dat)
    blocks = get_block_data(dat)

    refined_nx = 2 ** (h['levmax'] - 1) * h['domain_nx']
    domain_shape = np.append(refined_nx, h['nw'])
    d = np.zeros(domain_shape, order='F')

    max_lvl = h['levmax']
    # Get amount of blocks that need regridding
    print_regrid_amount(blocks, max_lvl)

    # Create multiprocessing iterable
    block_iterable = [(b, h) for b in blocks]
    # Create variable for multiprocess progress tracking
    init_progress = multiprocessing.Value("i", 0)
    mp_bool = multiprocessing.Value("i", True)
    total_blocks = multiprocessing.Value("i", len(blocks))

    #Initialize multiprocessing pool
    pool = multiprocessing.Pool(initializer=multiprocessing_init,
                                initargs=[init_progress, mp_bool, total_blocks],
                                processes=settings.nb_of_procs)

    print_tools.progress(0, 100, status='-- iterating over blocks...')
    #The aray blocks_regridded contains the data for each regridded block
    #:note: pool.(star)map obeys the array order during parallelization, i.e.
    #       blocks_regridded[i] equals the calculation for blocks[i]
    if h['ndim'] == 1:
        if PY2:
            blocks_regridded = np.array(pool.map(interpolate_block_1d_unpack, block_iterable))
        else:
            blocks_regridded = np.array(pool.starmap(interpolate_block_1d, block_iterable))
    elif h['ndim'] == 2:
        if PY2:
            blocks_regridded = np.array(pool.map(interpolate_block_2d_unpack, block_iterable))
        else:
            blocks_regridded = np.array(pool.starmap(interpolate_block_2d, block_iterable))
    else:
        if PY2:
            blocks_regridded = np.array(pool.map(interpolate_block_3d_unpack, block_iterable))
        else:
            blocks_regridded = np.array(pool.starmap(interpolate_block_3d, block_iterable))
    pool.close()
    pool.join()
    print_tools.progress(100, 100, status='-- completed.')
    print("\n")

    # Fill arrays with regridded data
    for i in range(len(blocks)):
        b = blocks[i]
        block_lvl = b['lvl']
        block_idx = b['ix']

        grid_diff = 2 ** (max_lvl - block_lvl)

        max_idx = block_idx * grid_diff
        min_idx = max_idx - grid_diff

        idx0 = min_idx * h['block_nx']

        if h['ndim'] == 1:
            if block_lvl == max_lvl:
                idx1 = idx0 + h['block_nx']
                d[idx0[0]:idx1[0], :] = b['w']
            else:
                idx1 = idx0 + (h['block_nx'] * grid_diff)
                d[idx0[0]:idx1[0], :] = blocks_regridded[i]

        elif h['ndim'] == 2:
            if block_lvl == max_lvl:
                idx1 = idx0 + h['block_nx']
                d[idx0[0]:idx1[0], idx0[1]:idx1[1], :] = b['w']
            else:
                idx1 = idx0 + (h['block_nx'] * grid_diff)
                d[idx0[0]:idx1[0], idx0[1]:idx1[1], :] = blocks_regridded[i]
        elif h['ndim'] == 3:
            if block_lvl == max_lvl:
                idx1 = idx0 + h['block_nx']
                d[idx0[0]:idx1[0], idx0[1]:idx1[1], idx0[2]:idx1[2], :] = b['w']
            else:
                idx1 = idx0 + (h['block_nx'] * grid_diff)
                d[idx0[0]:idx1[0], idx0[1]:idx1[1], idx0[2]:idx1[2], :] = blocks_regridded[i]
        else:
            raise IOError("Unknown number of dimensions %s" % h['ndim'])

    save_regridded_data(d)

    return d

