#!/usr/bin/env python

# Author: Jannis Teunissen
# This file should work in both Python 2 and 3

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import struct
import sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

def get_args():
    # Get and parse the command line arguments
    parser = argparse.ArgumentParser(
        description='''Plot .dat files of 1D MPI-AMRVAC simulations''',
        epilog='''Example:
        ./plot_1d.py ...''')
    parser.add_argument('dat_files', nargs='+',
                        type=argparse.FileType('rb'),
                        help='MPI-AMRVAC .dat files to read')
    parser.add_argument('-iw', type=int, default=1,
                        help='Index of variable to plot')
    parser.add_argument('-outfile', type=argparse.FileType('w'),
                        default=sys.stdout, help='MPI-AMRVAC .dat file')
    return parser.parse_args()

def main():
    ani = animation.FuncAnimation(fig, animate, args.dat_files, interval=50)
    plt.show()

def animate(f):
    data_1d = get_data(f, args)
    line.set_xdata(data_1d[:,0])
    line.set_ydata(data_1d[:,args.iw])

def get_data(dat, args):
    size_logical = 4
    size_int = 4
    size_double = 8

    fmt = 'iiiiiiiiiid'
    hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    [version, offset_tree, offset_blocks, nw, ndir, ndim,
     levmax, nleafs, nparents, it, time] = hdr

    fmt = 2 * 'd' + 2 * 'i'
    hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    [xmin, xmax, domain_nx, block_nx] = hdr

    # Read tree info. Skip 'leaf' array
    dat.seek(offset_tree + (nleafs+nparents) * size_logical)

    # Read block levels
    fmt = nleafs * 'i'
    block_lvls = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))

    # Read block indices
    fmt = nleafs * 'i'
    block_ixs = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))

    # Determine block coordinates
    dx0 = (xmax - xmin) / domain_nx

    ### Start reading data blocks
    dat.seek(offset_blocks)

    for i in range(nleafs):
        lvl = block_lvls[i]
        ix = block_ixs[i]

        x0 = (ix-1) * block_nx * dx0 * 0.5**(lvl-1)
        x1 = ix * block_nx * dx0 * 0.5**(lvl-1)

        # Read ghost cell info
        fmt = 2 * 'i'
        n_ghost = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))

        # Read actual data
        fmt = nw * block_nx * 'd'
        data = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
        print(x0, x1, data)

    sys.exit(0)

def calc_coord(c,n_cell,lb,dxc_block):
    # calculate the cell coordinates from the cell number
    local_ig = 1 + c % n_cell
    x = lb + (local_ig - 0.5) * dxc_block
    return x


args = get_args()
fig, ax = plt.subplots()
data_1d = get_data(args.dat_files[0], args)
line, = ax.plot(data_1d[:,0], data_1d[:,args.iw])
ax.set_xlabel('x')

if __name__ == "__main__":
    main()
