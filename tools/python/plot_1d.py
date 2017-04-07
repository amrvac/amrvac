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
        [h[name]] = val

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
    global i_var, i_file
    i_var += 1
    if i_var == h['nw'] + 1:
        i_var = 1
    update(i_file)

def prev_var(event):
    global i_var, i_file
    i_var -= 1
    if i_var == 0:
        i_var = h['nw']
    update(i_file)

def next_file(event):
    global i_file
    time_slider.set_val((i_file+1) % len(all_times))

def prev_file(event):
    global i_file
    time_slider.set_val((i_file+len(all_times)-1) % len(all_times))

def update(val):
    global i_file, i_var, is_conservative

    i_file = int(val+0.5)           # Get file index

    if is_conservative:
        line.set_data(all_cons[i_file][:,0], all_cons[i_file][:,i_var])
        ax.set_ylabel(cons_names[i_var])
    else:
        line.set_data(all_prim[i_file][:,0], all_prim[i_file][:,i_var])
        ax.set_ylabel(prim_names[i_var])

    ax.relim()
    ax.autoscale_view()
    ax.set_title('{} t = {:.3e}'.format(args.dat_files[i_file].name,
                                        all_times[i_file]))
    fig.canvas.draw_idle()

def switch_cons_prim(event):
    global i_file, is_conservative
    is_conservative = not is_conservative
    update(i_file)

def get_primitive_data(cc, wnames):
    pp = np.copy(cc)
    prim_names = wnames[:]
    ndir = h['ndir']
    i_rho = wnames.index('rho')

    # Convert momentum
    try:
        i_m1 = wnames.index('m1')
        for i in range(ndir):
            pp[:, i_m1+i] = pp[:, i_m1+i] / pp[:, i_rho]
            prim_names[i_m1+i] = 'v' + str(i+1)
    except ValueError:
        pass

    # Convert energy
    try:
        i_e = wnames.index('e')
        prim_names[i_e] = 'p'
        kin_en = 0.5 * pp[:, i_rho] * np.sum(pp[:, i_m1:i_m1+ndir]**2, axis=1)

        if h['physics_type'] == 'hd':
            pp[:, i_e] = (h['gamma'] - 1.0) * (pp[:, i_e] - kin_en)
        elif h['physics_type'] == 'mhd':
            i_b1 = wnames.index('b1')
            mag_en = 0.5 * np.sum(pp[:, i_b1:i_b1+ndir]**2, axis=1)
            pp[:, i_e] = (h['gamma'] - 1.0) * (pp[:, i_e] - kin_en - mag_en)
        else:
            print("Unknown physics type, cannot convert to pressure")
    except ValueError:
        pass

    return [pp, prim_names]

args = get_args()
all_cons = []
all_prim = []
all_times = []

# Read in all data
for f in args.dat_files:
    h, data_1d = get_data_1d(f, args)
    all_cons.append(data_1d)
    cons_names = ['x'] + h['w_names']
    prim_1d, prim_names = get_primitive_data(data_1d, cons_names)
    all_prim.append(prim_1d)
    all_times.append(h['time'])

fig, ax = plt.subplots()

# Add a time slider
if len(all_times) > 1:
    axfreq = plt.axes([0.1, 0.0, 0.5, 0.03])
    time_slider = Slider(axfreq, 'ix', 0., len(all_times)-1, 0.)
    time_slider.on_changed(update)

# Add buttons to go to next/previous variable and snapshot, and one to switch
# between primitive
axprev_iw = plt.axes([0.8, 0.0, 0.075, 0.05])
axnext_iw = plt.axes([0.875, 0.0, 0.075, 0.05])
axprim = plt.axes([0.95, 0.0, 0.05, 0.05])
axprev_file = plt.axes([0.0, 0.0, 0.05, 0.05])
axnext_file = plt.axes([0.7, 0.0, 0.05, 0.05])

bnext_iw = Button(axnext_iw, '+iw')
bnext_iw.on_clicked(next_var)
bprev_iw = Button(axprev_iw, '-iw')
bprev_iw.on_clicked(prev_var)
bnext_file = Button(axnext_file, '+')
bnext_file.on_clicked(next_file)
bprev_file = Button(axprev_file, '-')
bprev_file.on_clicked(prev_file)
bprim = Button(axprim, 'P/C')
bprim.on_clicked(switch_cons_prim)

# Initial state
i_var = 1
i_file = 0
is_conservative = True

line, = ax.plot(all_cons[0][:,0], all_cons[0][:,i_var])
ax.set_xlabel('x')
update(i_file)
plt.show()

