#!/usr/bin/env python

# Author: Jannis Teunissen
# This file should work in both Python 2 and 3

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import sys
import numpy as np

def get_args():
    # Get and parse the command line arguments
    parser = argparse.ArgumentParser(
        description='''Tool to compare two MPI-AMRVAC log files.
        An error message is shown if the physical results differ.''',
        epilog='''Example:
        ./compare_logfiles.py amrvac_1.log amrvac_2.log''')

    parser.add_argument('log', nargs=2, type=str,
                        help='The log files')
    parser.add_argument('-comments', type=str, default='#',
                        help='String indicating start of comments')
    parser.add_argument('-atol', type=float, default=1e-8,
                        help='Absolute tolerance for differences')
    parser.add_argument('-rtol', type=float, default=1e-5,
                        help='Relative tolerance for differences')

    return parser.parse_args()

if __name__ == '__main__':
    args = get_args()

    d = []

    # Don't show error messages, just return 0 (success) or 1 (fail)
    try:
        for fname in args.log:
            d.append(np.genfromtxt(fname, comments=args.comments))
        if not np.allclose(d[0], d[1], rtol=args.rtol, atol=args.atol):
            sys.exit(1)
    except:
        sys.exit(1)
