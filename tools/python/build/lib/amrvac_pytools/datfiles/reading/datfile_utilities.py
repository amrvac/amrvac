"""
Module containing reading and processing methods for an MPI-AMRVAC .datfiles file.

@author: Jannis Theunissen (original)
         Clément Robert (extensions, modifications)
"""
import struct
import numpy as np

# Size of basic types (in bytes)
SIZE_LOGICAL = 4
SIZE_INT = 4
SIZE_DOUBLE = 8
NAME_LEN = 16

# For un-aligned data, use '=' (for aligned data set to '')
ALIGN = '='


def get_header(istream):
    """Read header from an MPI-AMRVAC 2.1 snapshot. This is compatible with versions down to 2.0.
    :param: istream     open datfile buffer with 'rb' mode
    :return: h          header information contained in a dictionary
    """
    istream.seek(0)
    h = {}

    fmt = ALIGN + 'i'
    [h['datfile_version']] = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))

    if h['datfile_version'] < 3:
        raise IOError("Unsupported AMRVAC .datfiles file version: %d", h['datfile_version'])

    # Read scalar data at beginning of file
    fmt = ALIGN + 9 * 'i' + 'd'
    hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    [h['offset_tree'], h['offset_blocks'], h['nw'],
     h['ndir'], h['ndim'], h['levmax'], h['nleafs'], h['nparents'],
     h['it'], h['time']] = hdr

    # Read min/max coordinates
    fmt = ALIGN + h['ndim'] * 'd'
    h['xmin'] = np.array(
        struct.unpack(fmt, istream.read(struct.calcsize(fmt))))
    h['xmax'] = np.array(
        struct.unpack(fmt, istream.read(struct.calcsize(fmt))))

    # Read domain and block size (in number of cells)
    fmt = ALIGN + h['ndim'] * 'i'
    h['domain_nx'] = np.array(
        struct.unpack(fmt, istream.read(struct.calcsize(fmt))))
    h['block_nx'] = np.array(
        struct.unpack(fmt, istream.read(struct.calcsize(fmt))))

    if h['datfile_version'] >= 5:
        # Read periodicity
        fmt = ALIGN + h['ndim'] * 'i' # Fortran logical is 4 byte int
        h['periodic'] = np.array(
            struct.unpack(fmt, istream.read(struct.calcsize(fmt))), dtype=bool)

        # Read geometry name
        fmt = ALIGN + NAME_LEN * 'c'
        hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
        h['geometry'] = b''.join(hdr).strip().decode()

        # Read staggered flag
        fmt = ALIGN + 'i' # Fortran logical is 4 byte int
        h['staggered'] = bool(
            struct.unpack(fmt, istream.read(struct.calcsize(fmt)))[0])

    # Read w_names
    w_names = []
    for i in range(h['nw']):
        fmt = ALIGN + NAME_LEN * 'c'
        hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
        w_names.append(b''.join(hdr).strip().decode())
    h['w_names'] = w_names

    # Read physics type
    fmt = ALIGN + NAME_LEN * 'c'
    hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    h['physics_type'] = b''.join(hdr).strip().decode()

    # Read number of physics-defined parameters
    fmt = ALIGN + 'i'
    [n_pars] = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))

    # First physics-parameter values are given, then their names
    fmt = ALIGN + n_pars * 'd'
    vals = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))

    fmt = ALIGN + n_pars * NAME_LEN * 'c'
    names = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    # Split and join the name strings (from one character array)
    names = [b''.join(names[i:i+NAME_LEN]).strip().decode()
             for i in range(0, len(names), NAME_LEN)]

    # Store the values corresponding to the names
    for val, name in zip(vals, names):
        h[name] = val
    return h


def get_tree_info(istream):
    """
    Read levels, morton-curve indices, and byte offsets for each block as stored in the datfile
    :param: istream         open datfile buffer with 'rb' mode
    :return: block_lvls     numpy array with block levels
             block_ixs      numpy array with block morton-curve indices
             block_offsets  numpy array with block offset in the datfile
    """
    istream.seek(0)
    header = get_header(istream)
    nleafs = header['nleafs']
    nparents = header['nparents']

    # Read tree info. Skip 'leaf' array
    istream.seek(header['offset_tree'] + (nleafs+nparents) * SIZE_LOGICAL)

    # Read block levels
    fmt = ALIGN + nleafs * 'i'
    block_lvls = np.array(struct.unpack(fmt, istream.read(struct.calcsize(fmt))))

    # Read block indices
    fmt = ALIGN + nleafs * header['ndim'] * 'i'
    block_ixs = np.reshape(struct.unpack(fmt, istream.read(struct.calcsize(fmt))),
                           [nleafs, header['ndim']])

    # Read block offsets (skip ghost cells !)
    bcfmt = ALIGN + header['ndim'] * 'i'
    bcsize = struct.calcsize(bcfmt) * 2

    fmt = ALIGN + nleafs * 'q'
    block_offsets = np.array(struct.unpack(fmt, istream.read(struct.calcsize(fmt)))) + bcsize
    return block_lvls, block_ixs, block_offsets


def get_single_block_data(istream, byte_offset, block_shape):
    """"
    Retrieve a specific block from the datfile
    :param: istream       open datfile buffer in 'rb' mode
    :param: byte_offset   offset of the given block in the datfile
    :param: block_shape   the shape of the block (list containing dimensions + number of variables)
    :return: block_data   numpy array containing the block data, with dimensions equal to block_shape
    """
    istream.seek(byte_offset)
    # Read actual data
    fmt = ALIGN + np.prod(block_shape) * 'd'
    d = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    # Fortran ordering
    block_data = np.reshape(d, block_shape, order='F')
    return block_data


def get_blocks(dataset):
    """
    Reads all block data from an MPI-AMRVAC 2.0 snapshot.
    :param dataset   instance of 'amrvac_reader.load_file' class
    :return list containing block data as dictionaries with level, morton index and data.
    """
    blocks = []
    for ileaf, offset in enumerate(dataset.block_offsets):
        block = get_single_block_data(dataset.file, offset, dataset.block_shape)
        lvl = dataset.block_lvls[ileaf]
        ix = dataset.block_ixs[ileaf]
        b = {'lvl': lvl, 'ix': ix, 'w': block}
        blocks.append(b)

    return blocks



def get_uniform_data(dataset):
    """
    Retrieves the data for a uniform data set.
    :param dataset: instance of 'amrvac_reader.load_file' class
    :return The raw data as a NumPy array.
    """
    blocks = get_blocks(dataset)
    hdr = dataset.header

    refined_nx = 2 ** (hdr['levmax'] - 1) * hdr['domain_nx']
    domain_shape = np.append(refined_nx, hdr['nw'])
    d = np.zeros(domain_shape, order='F')

    for b in blocks:
        i0 = (b['ix'] - 1) * hdr['block_nx']
        i1 = i0 + hdr['block_nx']
        if hdr['ndim'] == 1:
            d[i0[0]:i1[0], :] = b['w']
        elif hdr['ndim'] == 2:
            d[i0[0]:i1[0], i0[1]:i1[1], :] = b['w']
        elif hdr['ndim'] == 3:
            d[i0[0]:i1[0], i0[1]:i1[1], i0[2]:i1[2], :] = b['w']
    return d
