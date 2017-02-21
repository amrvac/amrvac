#!/usr/bin/env python

# Author: Jannis Teunissen
# This file should work in both Python 2 and 3

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
# from __future__ import unicode_literals

import argparse
import struct
import sys
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np

size_logical = 4
size_int = 4
size_double = 8

def get_args():
    # Get and parse the command line arguments
    parser = argparse.ArgumentParser(
        description='''Plot .dat files of 1D MPI-AMRVAC simulations''',
        epilog='''Example:
        ./plot_1d.py ...''')
    parser.add_argument('dat_file', type=argparse.FileType('rb'),
                        help='MPI-AMRVAC .dat files to read')
    parser.add_argument('-iw', type=int, default=1,
                        help='Index of variable to plot')
    return parser.parse_args()

def read_header(dat):
    h = {}

    fmt = 'iiiiiiiiiid'
    hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    [h['version'], h['offset_tree'], h['offset_blocks'], h['nw'],
     h['ndir'], h['ndim'], h['levmax'], h['nleafs'], h['nparents'],
     h['it'], h['time']] = hdr

    fmt = 2 * 'd' + 2 * 'i'
    hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    [h['xmin'], h['xmax'], h['domain_nx'], h['block_nx']] = hdr

    # Read w_names
    w_names = []
    for i in range(h['nw']):
        fmt = 16 * 'c'
        hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
        w_names.append(''.join(hdr).strip())
    h['w_names'] = w_names

    fmt = 16 * 'c'
    hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    h['physics_type'] = ''.join(hdr).strip()

    return h

def get_data_1d(dat, args):
    h = read_header(dat)
    nw = h['nw']
    nx = h['block_nx']
    nleafs = h['nleafs']
    nparents = h['nparents']

    # Read tree info. Skip 'leaf' array
    dat.seek(h['offset_tree'] + (nleafs+nparents) * size_logical)

    # Read block levels
    fmt = nleafs * 'i'
    block_lvls = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))

    # Read block indices
    fmt = nleafs * 'i'
    block_ixs = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))

    # Determine coarse grid spacing
    dx0 = (h['xmax'] - h['xmin']) / h['domain_nx']

    ### Start reading data blocks
    dat.seek(h['offset_blocks'])

    out_data = np.zeros((nx * nleafs, nw+1))

    for i in range(nleafs):
        lvl = block_lvls[i]
        ix = block_ixs[i]

        x0 = (ix-1) * nx * dx0 * 0.5**(lvl-1)
        x1 = ix * nx * dx0 * 0.5**(lvl-1)
        x = np.linspace(x0, x1, nx)

        # Read number of ghost cells
        fmt = 2 * 'i'
        [gc_lo, gc_hi] = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))

        # Read actual data
        block_size = gc_lo + nw * nx + gc_hi
        fmt = block_size * nw * 'd'
        d = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
        w = np.reshape(d, [block_size, nw], 'F') # Fortran ordering

        out_data[i*nx:(i+1)*nx, 0] = x
        out_data[i*nx:(i+1)*nx, 1:] = w[gc_lo:gc_lo+nx, :]

    return h, out_data

args = get_args()
fig, ax = plt.subplots()
h, data_1d = get_data_1d(args.dat_file, args)
line, = ax.plot(data_1d[:,0], data_1d[:,args.iw])
ax.set_xlabel('x')
ax.set_ylabel(h['w_names'][args.iw-1])
plt.show()
