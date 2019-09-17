import sys
import multiprocessing
import scipy.interpolate as interp
import numpy as np

from amrvac_tools.datfiles.reading import datfile_utilities


def regrid_amr_data(istream, hdr, nbprocs):
    """
    Retrieves the data for a non-uniform data set by performing regridding.
    :param istream   open datfile buffer in 'rb' mode
    :param hdr       the .datfiles file header.
    :param nbprocs   the number of processors to use when regridding.
    :return: The raw data as a NumPy array.
    """
    # version check
    PY2 = sys.version_info[0] == 2
    if PY2:
        print("Regridding only has Python 3 support due to methods from the multiprocessing module.")
        sys.exit(1)

    if nbprocs is None:
        nbprocs = multiprocessing.cpu_count() - 2
    print("[INFO] Regridding using {} processors.".format(nbprocs))

    blocks = datfile_utilities.get_blocks(istream)
    refined_nx = 2 ** (hdr['levmax'] - 1) * hdr['domain_nx']
    domain_shape = np.append(refined_nx, hdr['nw'])
    d = np.zeros(domain_shape, order='F')

    max_lvl = hdr['levmax']
    block_iterable = [(b, hdr) for b in blocks]

    # Progress tracking
    init_progress = multiprocessing.Value("i", 0)
    nb_blocks  = multiprocessing.Value("i", len(blocks))

    print_progress(0, 100)

    # Initialize pool
    pool = multiprocessing.Pool(initializer=mp_init,
                                initargs=[init_progress, nb_blocks],
                                processes=nbprocs)

    # Execute multiprocessing pool
    blocks_regridded = np.array(pool.starmap(interpolate_block, block_iterable))
    pool.close()
    pool.join()
    print_progress(100, 100)
    print("")

    # fill arrays with regridded data
    for i in range(len(blocks)):
        b = blocks[i]
        block_lvl = b['lvl']
        block_idx = b['ix']

        grid_diff = 2 ** (max_lvl - block_lvl)

        max_idx = block_idx * grid_diff
        min_idx = max_idx - grid_diff

        idx0 = min_idx * hdr['block_nx']

        if hdr['ndim'] == 1:
            if block_lvl == max_lvl:
                idx1 = idx0 + hdr['block_nx']
                d[idx0[0]:idx1[0], :] = b['w']
            else:
                idx1 = idx0 + (hdr['block_nx'] * grid_diff)
                d[idx0[0]:idx1[0], :] = blocks_regridded[i]

        elif hdr['ndim'] == 2:
            if block_lvl == max_lvl:
                idx1 = idx0 + hdr['block_nx']
                d[idx0[0]:idx1[0], idx0[1]:idx1[1], :] = b['w']
            else:
                idx1 = idx0 + (hdr['block_nx'] * grid_diff)
                d[idx0[0]:idx1[0], idx0[1]:idx1[1], :] = blocks_regridded[i]
        else:
            if block_lvl == max_lvl:
                idx1 = idx0 + hdr['block_nx']
                d[idx0[0]:idx1[0], idx0[1]:idx1[1], idx0[2]:idx1[2], :] = b['w']
            else:
                idx1 = idx0 + (hdr['block_nx'] * grid_diff)
                d[idx0[0]:idx1[0], idx0[1]:idx1[1], idx0[2]:idx1[2], :] = blocks_regridded[i]

    return d


def print_progress(count, total):
    """
    Small method to print the current progress.
    :param count: Number of blocks done.
    :param total: Number of blocks in total.
    """
    percentage = round(100.0 * count / float(total), 1)
    print("Regridding...    {}%".format(percentage), end="\r")

def add_progress():
    """
    Adds progress to the multiprocessing variables
    """
    progress.value += 1
    if progress.value % 10 == 0:
        print_progress(progress.value, total_blocks.value)
    return

def mp_init(t, nb_blocks):
    """
    Initialiser method passed to the multiprocessing pool.
    :param t: progress
    :param nb_blocks: number of blocks
    """
    global progress, total_blocks
    progress = t
    total_blocks = nb_blocks

def interpolate_block(b, hdr):
    """
    Interpolates a given block to the maximum refinement level using flat interpolation.
    :param b: The block to refine.
    :param hdr: The .datfiles file header.
    :return: NumPy array containing the refined block data.
    """
    block_lvl = b['lvl']
    max_lvl   = hdr['levmax']
    if block_lvl == max_lvl:
        add_progress()
        return b
    ndim = hdr['ndim']
    curr_width = hdr['block_nx']

    grid_diff = 2**(max_lvl - block_lvl)
    regrid_width = curr_width * grid_diff
    nb_elements = np.prod(hdr['block_nx'])

    b_interpolated = np.zeros([*regrid_width, hdr['nw']])

    for var in range(0, hdr['nw']):
        if ndim == 1:
            block_spline = interp.interp1d(np.arange(curr_width), b['w'][:, var])
            block_result = block_spline(np.linspace(0, b['w'][:, var].size-1, regrid_width[0]))
            b_interpolated[:, var] = block_result
        elif ndim == 2:
            vals = np.reshape(b['w'][:, :, var], nb_elements)
            pts  = np.array(  [[i, j] for i in np.linspace(0, 1, curr_width[0])
                                      for j in np.linspace(0, 1, curr_width[1])]  )
            grid_x, grid_y = np.mgrid[0:1:regrid_width[0]*1j,
                                      0:1:regrid_width[1]*1j]
            grid_interpolated = interp.griddata(pts, vals, (grid_x, grid_y), method="linear")
            b_interpolated[:, :, var] = grid_interpolated
        else:
            vals = np.reshape(b['w'][:, :, :, var], nb_elements)
            pts  = np.array(  [[i, j, k] for i in np.linspace(0, 1, curr_width[0])
                                         for j in np.linspace(0, 1, curr_width[1])
                                         for k in np.linspace(0, 1, curr_width[2])]  )
            grid_x, grid_y, grid_z = np.mgrid[0:1:regrid_width[0]*1j,
                                              0:1:regrid_width[1]*1j,
                                              0:1:regrid_width[2]*1j]
            grid_interpolated = interp.griddata(pts, vals, (grid_x, grid_y, grid_z), method="linear")
            b_interpolated[:, :, :, var] = grid_interpolated

    add_progress()
    return b_interpolated

def regrid_2dmatrix(matrix, new_shape):
    if matrix.shape == tuple(new_shape):
        return matrix
    nb_elements = np.prod(matrix.shape)
    vals = np.reshape(matrix, nb_elements)
    pts = np.array([[i, j] for i in np.linspace(0, 1, matrix.shape[0])
                           for j in np.linspace(0, 1, matrix.shape[1])])
    grid_x, grid_y = np.mgrid[0:1:new_shape[0]*1j,
                              0:1:new_shape[1]*1j]
    block_regridded = interp.griddata(pts, vals, (grid_x, grid_y), method='linear')
    return block_regridded
