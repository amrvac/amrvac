#!/usr/bin/env python

# Author: Jannis Teunissen
# This file should work in both Python 2 and 3

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
# from __future__ import unicode_literals

import argparse
import struct
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import numpy as np

# Size of basic types (in bytes)
size_logical = 4
size_int = 4
size_double = 8

# For un-aligned data, use '=' (for aligned data set to '')
align = '='

def get_args():
    # Get and parse the command line arguments
    parser = argparse.ArgumentParser(
        description='''Plot .dat files of 1D MPI-AMRVAC simulations''',
        epilog='''Example:
        ./plot_1d.py ...''')
    parser.add_argument('dat_files', nargs='+',
                        type=argparse.FileType('rb'),
                        help='MPI-AMRVAC .dat files to read')
    return parser.parse_args()

def read_header(dat):
    h = {}

    fmt = align + 'iiiiiiiiiid'
    hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    [h['version'], h['offset_tree'], h['offset_blocks'], h['nw'],
     h['ndir'], h['ndim'], h['levmax'], h['nleafs'], h['nparents'],
     h['it'], h['time']] = hdr

    fmt = align + 2 * 'd' + 2 * 'i'
    hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    [h['xmin'], h['xmax'], h['domain_nx'], h['block_nx']] = hdr

    # Read w_names
    w_names = []
    for i in range(h['nw']):
        fmt = align + 16 * 'c'
        hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
        w_names.append(''.join(hdr).strip())
    h['w_names'] = w_names

    fmt = align + 16 * 'c'
    hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    h['physics_type'] = ''.join(hdr).strip()

    fmt = align + 'i'
    n_pars, = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))

    for i in range(n_pars):
        fmt = align + 'd'
        val = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
        fmt = align + 16 * 'c'
        name = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
        name = ''.join(name).strip()
        h[name] = val

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
    fmt = align + nleafs * 'i'
    block_lvls = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))

    # Read block indices
    fmt = align + nleafs * 'i'
    block_ixs = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))

    # Determine coarse grid spacing
    dx0 = (h['xmax'] - h['xmin']) / h['domain_nx']

    ### Start reading data blocks
    dat.seek(h['offset_blocks'])

    out_data = np.zeros((nx * nleafs, nw+1))

    for i in range(nleafs):
        lvl = block_lvls[i]
        ix = block_ixs[i]

        x0 = (ix-1) * nx * dx0 * 0.5**(lvl-1) + dx0 * 0.5**lvl + h['xmin']
        x1 = ix * nx * dx0 * 0.5**(lvl-1) - dx0 * 0.5**lvl + h['xmin']
        x = np.linspace(x0, x1, nx)

        # Read number of ghost cells
        fmt = align + 2 * 'i'
        [gc_lo, gc_hi] = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))

        # Read actual data
        block_size = gc_lo + nx + gc_hi
        fmt = align + block_size * nw * 'd'
        d = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
        w = np.reshape(d, [block_size, nw], 'F') # Fortran ordering
        out_data[i*nx:(i+1)*nx, 0] = x
        out_data[i*nx:(i+1)*nx, 1:] = w[gc_lo:gc_lo+nx, :]

    return h, out_data

def next_var(event):
    global i_var
    i_var = (i_var + 1) % h['nw']
    update(i_file)

def prev_var(event):
    global i_var
    i_var = (i_var - 1) % h['nw']
    update(i_file)

def update(val):
    i_file = int(val+0.5)           # Get file index
    line.set_data(all_data[i_file][:,0], all_data[i_file][:,i_var])
    ax.relim()
    ax.autoscale_view()
    ax.set_ylabel(h['w_names'][i_var-1])
    ax.set_title('{} t = {:.3e}'.format(args.dat_files[i_file].name,
                                        all_times[i_file]))
    fig.canvas.draw_idle()

args = get_args()
all_data = []
all_times = []

# Read in all data
for f in args.dat_files:
    h, data_1d = get_data_1d(f, args)
    all_data.append(data_1d)
    all_times.append(h['time'])

fig, ax = plt.subplots()

# Add a time slider
if len(all_times) > 1:
    axfreq = plt.axes([0.1, 0.0, 0.5, 0.03])
    sfreq = Slider(axfreq, 'index', 0., len(all_times)-1, 0.)
    sfreq.on_changed(update)

# Add buttons to go to next/previous variable
axprev = plt.axes([0.7, 0.0, 0.07, 0.05])
axnext = plt.axes([0.8, 0.0, 0.07, 0.05])
bnext = Button(axnext, '+iw')
bnext.on_clicked(next_var)
bprev = Button(axprev, '-iw')
bprev.on_clicked(prev_var)

# Initial state
i_var = 1
i_file = 0

line, = ax.plot(all_data[0][:,0], all_data[0][:,i_var])
ax.set_title('{} t = {}'.format(args.dat_files[0].name, all_times[0]))
ax.set_xlabel('x')
ax.set_ylabel(h['w_names'][i_var-1])

plt.show()

