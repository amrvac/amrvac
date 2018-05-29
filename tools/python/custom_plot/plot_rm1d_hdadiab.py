#!/usr/bin/env python

# Author: Jannis Teunissen
# This file should work in both Python 2 and 3

# Modified by Wenzhi Ruan
# The function to get the analytic solution of Riemann 1D problem is added

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
# from __future__ import unicode_literals

import sys
import math
import argparse
import struct
import matplotlib.pyplot as plt
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
        description='''Plot .dat files of 1D MPI-AMRVAC simulations [for the results of adiabatic hd module only]. This code will read the initial condition at %0000.dat automatically and work out the analytic solutions of Riemann problem. Please set itsave(1,2)=0 in .par files before you run the simulation to output the initial condition into %0000.dat.''',
        epilog='''Example:
        ./plot_1d.py ...''')
    parser.add_argument('dat_file', type=argparse.FileType('rb'),
                        help='MPI-AMRVAC .dat files to read')
    parser.add_argument('-iw', type=int, default=1,
                        help='Index of variable to plot')
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


#get analytic solutions of Riemann 1D problem
def get_analytic_solution(t, x_t, args_t):
    #get data of step0
    filename_t0 = get_filename_t0(args_t) 
    file_t0 = open(filename_t0, 'r')
    h_t0, data_t0 = get_data_1d(file_t0, args_t)
    file_t0.close()
    t0 = h_t0['time']
    global gamma
    x_t0 = data_t0[:,0]
    rho_t0 = data_t0[:,1]
    v_t0 = data_t0[:,2]/rho_t0[:]

    if (t0==t):
        return x_t0, rho_t0, v_t0
    else:
        #three region of adiabatic Riemann 1D problem: left|center|right
        #get the values of rho, v at the left & right region
        num_cell_t0 = len(x_t0)
        rho_left = rho_t0[0]
        rho_right = rho_t0[num_cell_t0-1]
        v_left = v_t0[0]
        v_right = v_t0[num_cell_t0-1]
        cs_left = math.sqrt(gamma*rho_left**(gamma-1))
        cs_right = math.sqrt(gamma*rho_right**(gamma-1))
        rho_center = get_rho_center(rho_left, v_left, rho_right, v_right)
        cs_center = math.sqrt(gamma*rho_center**(gamma-1))

        #location of discontinuity
        x_dis = get_discountinuity(x_t0, rho_t0)
        if (x_dis == -1):
            x_dis = get_discountinuity(x_t0, v_t0)
        if (x_dis == -1):
            x_dis = 0
            print('no initial discountinuity!')

        #region x/t < v_center
        if (rho_center >= rho_left):
            #shock
            f_left = get_f(rho_center, rho_left)
            v_center = v_left - f_left
            x_center = v_center*t + x_dis
            Z_left = (rho_center*v_center-rho_left*v_left)/(rho_center-rho_left)
            x_left = x_dis + Z_left*t
        else:
            #rarefaction wave
            f_left = get_f(rho_center, rho_left)
            v_center = v_left - f_left
            x_center = v_center*t + x_dis
            Z_left1 = v_left - cs_left
            Z_left2 = v_center - cs_center
            x_left1 = x_dis + Z_left1*t
            x_left2 = x_dis + Z_left2*t
            x_left = x_left1

        #region x/t > v_center
        if (rho_center >= rho_right):
            #shock
            f_right = get_f(rho_center, rho_right)
            Z_right = (rho_center*v_center-rho_right*v_right)/(rho_center-rho_right)
            x_right = x_dis + Z_right*t
        else:
            #rarefaction wave
            f_left = get_f(rho_center, rho_left)
            Z_right1 = v_right + cs_right
            Z_right2 = v_center + cs_center
            x_right1 = x_dis + Z_right1*t
            x_right2 = x_dis + Z_right2*t
            x_right = x_right1

        #analytic solutions
        x_analy = x_t
        len_analy = len(x_analy)
        rho_analy = np.zeros(len_analy)
        v_analy = np.zeros(len_analy)

        for i in range(0, len_analy):
            if (x_analy[i] < x_left):
                rho_analy[i] = rho_left
                v_analy[i] = v_left

            elif (x_right <= x_analy[i]):
                rho_analy[i] = rho_right
                v_analy[i] = v_right

            elif (x_left <= x_analy[i] < x_center):
                if (rho_center >= rho_left):
                    #shock
                    rho_analy[i] = rho_center
                    v_analy[i] = v_center
                else:
                    #rarefaction wave
                    if (x_left2 <= x_analy[i]): 
                        #center region
                        rho_analy[i] = rho_center
                        v_analy[i] = v_center
                    else: 
                        #wave region
                        temp = (gamma-1)*(v_left-(x_analy[i]-x_dis)/t)/(gamma+1)
                        cs_rare = temp + 2*cs_left/(gamma+1)
                        v_analy[i] = (x_analy[i]-x_dis)/t + cs_rare
                        rho_analy[i] = (cs_rare**2/gamma)**(1/(gamma-1))

            else:    #x_center <= x_analy[i] < x_right
                if (rho_center > rho_right):
                    #shock
                    rho_analy[i] = rho_center
                    v_analy[i] = v_center
                else:    #rarefaction wave
                    if (x_analy[i] <= x_right2):
                        #center region
                        rho_analy[i] = rho_center
                        v_analy[i] = v_center
                    else:
                        #wave region
                        temp1 = (x_analy[i]-x_dis)/t-v_right
                        temp2 = (gamma-1)/(gamma+1) * temp1
                        cs_rare = temp2 + 2*cs_right/(gamma+1)
                        v_analy[i] = (x_analy[i]-x_dis)/t - cs_rare
                        rho_analy[i] = (cs_rare**2/gamma)**(1/(gamma-1))

        return x_analy, rho_analy, v_analy


#get rho in center region
def get_rho_center(rho_l, v_l, rho_r, v_r):
    rho1 = 0
    rho2 = 100*(rho_l+rho_r)
    F = get_f((rho1+rho2)/2, rho_l) + get_f((rho1+rho2)/2, rho_r) - v_l + v_r
    k = 0
    while (abs(F) > 0.0001*(abs(v_l)+abs(v_r)+1) and k < 50):
        if (F > 0):
            rho2 = (rho1+rho2)/2
        else:
            rho1 = (rho1+rho2)/2
        F = get_f((rho1+rho2)/2, rho_l) + get_f((rho1+rho2)/2, rho_r)
        F = F - v_l + v_r
        k = k + 1
        if (k > 30):
            print('please check the values of p1 & p2 in function get_pcenter!')
            sys.exit()

    return (rho1+rho2)/2



def get_f(rho_cent, rho_j):
    global gamma
    global c_adiab
    if (rho_cent >= rho_j):
        temp_j = (rho_cent-rho_j)*(rho_cent**gamma-rho_j**gamma)
        f_j = math.sqrt(c_adiab*temp_j/(rho_cent*rho_j))
    else:
        temp_j = rho_cent**((gamma-1)/2) - rho_j**((gamma-1)/2)
        f_j = 2*math.sqrt(gamma)*temp_j/(gamma-1)

    return f_j

    

#get the name of output file at step0
def get_filename_t0(args_temp):
    string = "{}".format(args_temp.dat_file)
    lenstring = len(string)
    for i in range(0, lenstring):
        if (string[i]=="'"):
            i1 = i
            break
    for i in range(i1+1, lenstring-i1-1):
        if (string[i]=="'"):
            i2 = i
            break
    name_temp = string[i1+1:i2]

    len_name = len(name_temp)
    for i in range(0, len_name):
        if (name_temp[i]=='.'):
            i3 = i
            break
    name = name_temp[0:i3-4] + '0000' + name_temp[i3:len_name]

    return name


#get initial condition
def get_initial(args_t):
    filename_t0 = get_filename_t0(args_t) 
    file_t0 = open(filename_t0, 'r')
    h_t0, data_t0 = get_data_1d(file_t0, args_t)
    file_t0.close()
    global gamma
    x_t0 = data_t0[:,0]
    rho_t0 = data_t0[:,1]
    v_t0 = data_t0[:,2]/rho_t0[:]

    return x_t0, rho_t0, v_t0


#find initial discountinuity
def get_discountinuity(x_in, y_in):
    lenx = len(x_in)
    ymax = max(y_in[0], y_in[lenx-1])
    ymin = min(y_in[0], y_in[lenx-1])
    ydiff = ymax - ymin
    x_d = -1.0

    for i in range(1, lenx-1):
        if (y_in[i-1] < ymin+0.5*ydiff):
            if (ymin+0.5*ydiff <= y_in[i]):
                x_d = x_in[i]
                break
        if (ymin+0.5*ydiff < y_in[i-1]):
            if (y_in[i+1] <= ymin+0.5*ydiff):
                x_d = x_in[i]
                break

    return x_d


#get scheme name
def get_name_scheme(args_temp):
    string = "{}".format(args_temp.dat_file)
    lenstring = len(string)
    for i in range(0, lenstring):
        if (string[i]=="'"):
            i1 = i
            break
    for i in range(i1+1, lenstring-i1-1):
        if (string[i]=="'"):
            i2 = i
            break
    name_temp = string[i1+1:i2]
    name = name_temp[6:len(name_temp)-8]

    return name




args = get_args()
h, data_1d = get_data_1d(args.dat_file, args)

time = h['time']

#convert variables
global c_adiab
c_adiab = 1
global gamma
gamma = h['gamma'][0]
x = data_1d[:,0]
rho = data_1d[:,1]
v = data_1d[:,2]/rho[:]

#get the analytic solution of Riemann 1D problem
x_analytic, rho_analytic, v_analytic = get_analytic_solution(time, x, args)

#get initial condition
x_initial, rho_initial, v_initial = get_initial(args)

scheme = get_name_scheme(args)

fig, ax = plt.subplots()
#plt.suptitle('analytic: __     numerical: **    initial: ..', fontsize=18)
#line, = ax.plot(data_1d[:,0], data_1d[:,args.iw])
if (args.iw==1):    
    line, = ax.plot(x, rho, 'k.', label='numerical')
    line, = ax.plot(x_analytic, rho_analytic, 'r', label='analytic')
    line, = ax.plot(x_initial, rho_initial, 'b:', label='initial')
    plt.ylim([0.4, 1.601])
    ax.set_ylabel('rho')
    w_name = 'rho'
    ymax, ymin = max(rho), min(rho)
elif (args.iw==2):    
    line, = ax.plot(x, v, 'k.', label='numerical')
    line, = ax.plot(x_analytic, v_analytic, 'r', label='analytic')
    line, = ax.plot(x_initial, v_initial, 'b:', label='initial')
    plt.ylim([-0.101, 0.801])
    ax.set_ylabel('v')
    w_name = 'v'
    ymax, ymin = max(v), min(v)
else:
    print('Error! Wrong iw!') 
    sys.exit()   

if (ymax == ymin):
    plt.ylim([ymin-0.5, ymax+0.5])
else:
    plt.ylim([1.2*ymin-0.2*ymax, 1.2*ymax-0.2*ymin])

ax.set_xlabel('x')
ax.set_title('{}    t = {} s'.format(scheme, time))
plt.legend(loc=1)
#ax.set_ylabel(h['w_names'][args.iw-1])
#plt.savefig('rm_1d_{}_{}.png'.format(scheme, w_name))
plt.show()
