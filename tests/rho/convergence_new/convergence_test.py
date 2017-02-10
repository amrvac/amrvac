#!/usr/bin/env python

# Author: Jannis Teunissen
# This file should work in both Python 2 and 3

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import vtk
from tvtk import array_handler as ah
from multiprocessing import Pool, Manager, cpu_count
from subprocess import check_output
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def get_args():
    p = argparse.ArgumentParser()

    p.add_argument('-only_plot', action='store_true',
                   help='Only plot, do not run amrvac')
    p.add_argument('-np', type=int,
                   default=max(1, cpu_count()//2),
                   help='Number of parallel amrvac proccesses')
    p.add_argument('-fig_size', type=float, nargs=2,
                   default=[10,6], help='Figure size in inches')
    p.add_argument('-fig_name', type=str, default='result.pdf',
                   help='Store the result in this figure')
    p.add_argument('-show_powers', type=float, nargs='+',
                   default=[2.0],
                   help='Show these powers in the graph')
    p.add_argument('-base_par', type=str, default='conv.par',
                   help='Base .par file for the runs')
    p.add_argument('-grids', type=int, nargs='+',
                   default=[16, 32, 64, 128, 256, 512, 1024, 2048],
                   help='Grid sizes to use (.par files have to exist)')
    p.add_argument('-scheme_dir', type=str,
                   default='../../schemes',
                   help='Directory for the schemes')
    p.add_argument('-norm', type=str, default='inf',
                   help='Which norm to use (1, 2, inf)')
    p.add_argument('-schemes', type=str, nargs='+',
                   default=['1step_tvd', '2step_tvdlf_mm',
                            '2step_tvdmu_al', '3step_hll_cada',
                            '4step_hll_mc', '4step_hllc_ko',
                            'rk4_tvdlf_cada', 'ssprk43_fd_mp5',
                            'ssprk54_fd_mp5'],
                   help='Which schemes to use')
    return p.parse_args()

def load(file):
    datareader = vtk.vtkXMLImageDataReader()
    datareader.SetFileName(file)
    datareader.Update()
    return datareader.GetOutput()


def getX(data):
    dx=data.GetSpacing()[0]
    nx=data.GetExtent()[1]
    x=np.arange(0,1,dx)+dx/2.
    return x


def getRho(data):
    rho = ah.vtk2array(data.GetCellData().GetArray('rho'))
    return rho


def plot_rho(data):
    x=getX(data)
    rho=getRho(data)
    plt.plot(x,rho,'k+-')


def getDifference(file_a, file_b, norm_type):
    data_a = load(file_a)
    rho_a = getRho(data_a)

    data_b = load(file_b)
    rho_b = getRho(data_b)

    if norm_type == 'inf':
        norm = np.linalg.norm(rho_a-rho_b, np.inf)
    elif norm_type == '1':
        norm = np.linalg.norm(rho_a-rho_b, 1) / rho_a.size

    return norm

# Wrapper script to run amrvac
def amrvac_wrapper(par_files):
    cmd = ['./amrvac', '-i ' + par_files]
    try:
        res = check_output(cmd)
    except:
        print('Error when running:')
        print(' '.join(cmd))
        sys.exit(1)
    return res

if __name__ == '__main__':
    args = get_args()

    grids = np.array(args.grids)
    scheme_pars = [args.scheme_dir + '/' + sch + '.par' for
                   sch in args.schemes]

    if not args.only_plot:
        pool = Pool(processes=args.np)
        cmd_list = []
        for sch, par_sch  in zip(args.schemes, scheme_pars):
            for res in grids:
                par_files = '{} {:d}.par {}'.format(
                    args.base_par, res, par_sch)
                cmd_list.append(par_files)

        run = pool.map(amrvac_wrapper, cmd_list)

    plt.figure(figsize=args.fig_size)

    for sch, par_sch in zip(args.schemes, scheme_pars):
        diff = []
        for res in grids:
            fname0 = 'conv_{:d}_{}{:04d}.vti'.format(res, sch, 0)
            fname1 = 'conv_{:d}_{}{:04d}.vti'.format(res, sch, 1)
            diff.append(getDifference(fname0, fname1, args.norm))
        line, = plt.loglog(grids, diff, linestyle='-',
                           marker='.', label=sch)

    for p in args.show_powers:
        plt.loglog(grids, 0.1 * (float(grids[0])/grids)**p,
                   linestyle='--', label='order-{:.1f}'.format(p))

    plt.ylabel('Error {}-norm'.format(args.norm))
    plt.xlabel('Grid points')
    plt.legend(loc='lower left')
    plt.xlim([grids[0]*0.9, grids[-1]*1.1])
    plt.savefig(args.fig_name)
    print('Saved figure {}'.format(args.fig_name))
