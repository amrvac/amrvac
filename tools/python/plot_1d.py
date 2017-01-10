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

    dat.seek(0,2)       # goto EOF
    offset = -32        # 32 = 6*size_int + size_double = 6*4 + 8
    dat.seek(offset,1)  # go back 36 bytes

    ### Read 6 footer ints + 1 double
    int_in = dat.read(4)
    nleafs = struct.unpack('i',int_in)[0]    # number of active tree leafs nleafs (= #blocks)
    int_in = dat.read(4)
    levmaxini = struct.unpack('i',int_in)[0] # maximal refinement level present levmax
    int_in = dat.read(4)
    ndimini = struct.unpack('i',int_in)[0]   # dimensionality NDIM
    int_in = dat.read(4)
    ndirini = struct.unpack('i',int_in)[0]   # number of vector components NDIR
    int_in = dat.read(4)
    nwini = struct.unpack('i',int_in)[0]     # number of variables nw
    int_in = dat.read(4)
    it = struct.unpack('i',int_in)[0]        # integer time counter it
    flt_in = dat.read(8)
    t = struct.unpack('d',flt_in)[0]         # global time t

    # TODO: Check ndimini

    # read block size in each dimension and all eqpars
    offset = offset-int(2*ndimini*4 + 2*ndimini*8) # int size = 4, double = 8
    dat.seek(offset,1)

    # Get block size
    int_in = dat.read(4)
    n_cell = struct.unpack('i',int_in)[0]
    size_block = n_cell*nwini*8   # block size in bytes
    # Get size of level one mesh
    int_in = dat.read(4)
    n_mesh_lev1= struct.unpack('i',int_in)[0]
    # Get x minimum location
    flt_in = dat.read(8)
    xmin = struct.unpack('d',flt_in)[0]
    # Get x maximum location
    flt_in = dat.read(8)
    xmax = struct.unpack('d',flt_in)[0]
    

    # Read forest
    dat.seek(nleafs*size_block,0)   # go to end of block stuctures, start of the forest

    n_blocks = n_mesh_lev1 // n_cell     # number of blocks in each direction

    block_info = []
    igrid = 0
    refine = 0
    forest = []
    level = 1

    for i in range(1, n_blocks+1):
        (igrid,refine) = read_node(dat,forest,igrid,refine,ndimini,level,block_info, i)

    outfile = args.outfile

    # Calculate physical sizes of blocks and cells on level one

    # physical block size in each direction (on the first level)
    dx = (xmax - xmin) / n_blocks

    # physical cell size in each direction (on the first level)
    dxc = (xmax - xmin) / n_mesh_lev1

    data_1d = np.zeros((n_cell * nleafs, nwini+1))
    i_cell = 0
    ### Start reading data blocks

    dat.seek(0,0)  # go to start of file
    for i in range(nleafs):
        # coordinates of bottom left corner of block
        lb = xmin + (block_info[i][1] - 1) * dx / (2**(block_info[i][0]-1))

        # local cell size
        dxc_block = dxc/(2**(block_info[i][0]-1))
        values = [ [] for i in range(nwini) ]
        for nw in range(nwini):
            for c in range(n_cell):
                flt_in = dat.read(8)
                dbl = struct.unpack('d',flt_in)[0]
                values[nw].append(dbl)

        for c in range(n_cell):
            x = calc_coord(c,n_cell,lb,dxc_block)
            data_1d[i_cell][0] = x
            data_1d[i_cell][1:] = values[nw][c]
            i_cell += 1

    return data_1d

def read_node(dat,forest,igrid,refine,ndimini,level,block_info,ix):
    # Recursive function that checks if the block is a leaf or not.
    # If leaf, write block info. If not, run recursively for all
    # possible child nodes

    b = dat.read(4)
    leaf = struct.unpack('i',b)[0]

    if (leaf):
        forest.append(1)
        igrid += 1
        block_info.append([level,ix])
    else:
        forest.append(0)
        refine += 1
        for i in range(2):
            child_index = new_index(ndimini, i, ix)
            (igrid,refine) = read_node(dat,forest,igrid,refine,ndimini,level+1,block_info,child_index)
    return igrid,refine


def new_index(dims, i, ig):
    # Calculate new grid indices for child nodes if node is refined.
    # See Keppens (2012) section 3.3
    return 2 * (ig-1) + 1 + i


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
