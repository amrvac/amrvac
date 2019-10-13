"""
Script for various data manipulations, dictionary creation, primitive variable addition etc.

@author: Niels Claes
"""

import sys
import numpy as np
import copy
from collections import OrderedDict


def get_block_edges(ileaf, dataset):
    lvl = dataset.block_lvls[ileaf]
    morton_idx = dataset.block_ixs[ileaf]
    block_nx = dataset.header["block_nx"]

    # dx at coarsest grid level
    dx0 = dataset.domain_width / dataset.header["domain_nx"]
    # dx at current block level
    dx = dx0 * 0.5 ** (lvl - 1)

    # calculate actual edges of the block
    l_edge = dataset.header["xmin"] + (morton_idx - 1) * block_nx * dx
    r_edge = l_edge + block_nx * dx
    return l_edge, r_edge

def create_data_dict(raw_data, header, adiab_cte=1):
    default_fields = copy.deepcopy(header['w_names'])
    data_dict = OrderedDict()

    for var in default_fields:
        idx = default_fields.index(var)
        data_dict[var] = raw_data[..., idx]

    # check if gamma is present (hd or mhd), otherwise return dictionary
    if not 'gamma' in header.keys():
        return data_dict

    # add primitive variables to the data dictionary
    # pressure is always present for hd/mhd
    p = _pressure(data_dict, header, adiab_cte)
    data_dict.update({'p': p})

    # add velocities if momenta are present
    if 'm1' in default_fields:
        v1 = _velocity1(data_dict)
        data_dict.update({'v1': v1})
    if 'm2' in default_fields:
        v2 = _velocity2(data_dict)
        data_dict.update({'v2': v2})
    if 'm3' in default_fields:
        v3 = _velocity3(data_dict)
        data_dict.update({'v3': v3})

    # add temperature if there is an energy equation
    if 'e' in default_fields:
        T = data_dict['p'] / data_dict['rho']
        data_dict.update({'T': T})
    return data_dict


def _get_ekin(data_dict, header):
    ndim = header['ndim']

    ekin = 0.5 * data_dict['m1']**2 / data_dict['rho']
    if ndim > 1:
        ekin = ekin + 0.5 * data_dict['m2'] / data_dict['rho']
    if ndim > 2:
        ekin = ekin + 0.5 * data_dict['m3'] / data_dict['rho']
    return ekin


def _get_emag(data_dict, header):
    ndim = header['ndim']

    emag = 0.5 * data_dict['b1']**2
    if ndim > 1:
        emag = emag + 0.5 * data_dict['b2']**2
    if ndim > 2:
        emag = emag + 0.5 * data_dict['b3']**2
    return emag


def _pressure(data_dict, header, adiab_cte=1):
    gamma = header['gamma']
    phys_type = header['physics_type']
    default_fields = copy.deepcopy(header['w_names'])

    # pressure for HD/MHD datasets
    if 'e' in default_fields:
        if phys_type == 'hd':
            return (gamma - 1) * (data_dict['e'] - _get_ekin(data_dict, header))
        else:
            return (gamma - 1) * (data_dict['e'] - _get_ekin(data_dict, header) - _get_emag(data_dict, header))
    else:
        return adiab_cte * data_dict['rho'] ** gamma


def _velocity1(data_dict):
    return data_dict['m1'] / data_dict['rho']
def _velocity2(data_dict):
    return data_dict['m2'] / data_dict['rho']
def _velocity3(data_dict):
    return data_dict['m3'] / data_dict['rho']