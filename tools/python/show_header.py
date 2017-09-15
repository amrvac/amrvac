#!/usr/bin/env python

# Author: Jannis Teunissen
# This file should work in both Python 2 and 3

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
# from __future__ import unicode_literals

import argparse
import struct
import array

# For un-aligned data, use '=' (for aligned data set to '')
align = '='

def get_args():
    # Get and parse the command line arguments
    parser = argparse.ArgumentParser(
        description='''Show header of MPI-AMRVAC .dat file''')
    parser.add_argument('dat_file',
                        type=argparse.FileType('rb'),
                        help='MPI-AMRVAC .dat file to read')
    return parser.parse_args()

def read_header(dat):
    h = {}

    fmt = align + 'iiiiiiiiiid'
    hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    [h['version'], h['offset_tree'], h['offset_blocks'], h['nw'],
     h['ndir'], h['ndim'], h['levmax'], h['nleafs'], h['nparents'],
     h['it'], h['time']] = hdr

    ndim = h['ndim']
    fmt = align + 2 * ndim * 'd' + 2 * ndim * 'i'
    hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    h['xmin'] = hdr[0:ndim]
    h['xmax'] = hdr[ndim:2*ndim]
    h['domain_nx'] = hdr[2*ndim:3*ndim]
    h['block_nx'] = hdr[3*ndim:4*ndim]

    # Read w_names
    w_names = []
    for i in range(h['nw']):
        fmt = align + 16 * 'c'
        hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
        name = b''.join(hdr).decode("utf-8").strip()
        w_names.append(name)
    h['w_names'] = w_names

    fmt = align + '16s'
    hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    name = b''.join(hdr).decode("utf-8").strip()
    h['physics_type'] = name

    fmt = align + 'i'
    n_pars, = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))

    fmt = align + n_pars * 'd'
    vals = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    fmt = align + n_pars * 16 * 'c'
    names = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    # Split and join the name strings (from one character array)
    names = [b''.join(names[i:i+16]).decode("utf-8").strip()
             for i in range(0, len(names), 16)]

    for val, name in zip(vals, names):
        h[name] = val

    return h

if __name__ == '__main__':
    args = get_args()
    h = read_header(args.dat_file)

    print("{:20s} | {}".format("parameter", "value(s)"))
    print("----------------------------------------")

    for k,v in h.items():
    	print("{:20s} = {}".format(k,v))
