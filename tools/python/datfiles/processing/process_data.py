"""
Class to process the MPI-AMRVAC raw data into a useful Python dictionary.

@author: Niels Claes
"""

import sys
import numpy as np
import copy
from amrvac_tools.datfiles.processing import convert
from collections import OrderedDict


def add_primitives_to_single_block(block, dataset, add_velocities=False):
    """
    This function is used by the synthetic views to calculate pressure and temperature on-the-go.
    For now, only add pressure and temperature
    :param block: data of a single block, numpy array of shape dataset.block_shape
    :param dataset: amrvac_reader instance, containing the dataset info
    :param add_velocities: if True, adds velocities to the block as well. Default is false
    :return: block with pressure and temperature added
    """
    # calculate pressure
    ndim = dataset.header['ndim']
    phys_type = dataset.header['physics_type']
    block_fields = copy.deepcopy(dataset.header['w_names'])
    if not (phys_type == 'hd' or phys_type == 'mhd'):
        return block, block_fields

    gamma = dataset.header['gamma']
    rho_idx = dataset.header['w_names'].index('rho')
    if ndim == 2:
        if phys_type == 'hd':
            v1, v2, p = convert.hd_to_primitive_2d(*(block[:, :, idx] for idx in range(0, len(block_fields))), gamma)
        else:
            v1, v2, p = convert.mhd_to_conserved_2d(*(block[:, :, idx] for idx in range(0, len(block_fields))), gamma)
        temp = p / block[:, :, rho_idx]
    else:
        if phys_type == 'hd':
            v1, v2, v3, p = convert.hd_to_primitive_3d(*(block[:, :, :, idx] for idx in range(0, len(block_fields))), gamma)
        else:
            v1, v2, v3, p = convert.mhd_to_primitive_3d(*(block[:, :, :, idx] for idx in range(0, len(block_fields))), gamma)
    temp = p / block[..., rho_idx]
    # add pressure and temperature to block_fields for consistent retrieval of index later on
    block = np.concatenate((block, p[..., np.newaxis]), axis=ndim)      # add pressure
    block = np.concatenate((block, temp[..., np.newaxis]), axis=ndim)   # add temperature
    block_fields += ['p', 'T']
    if add_velocities:
        block = np.concatenate((block, v1[..., np.newaxis]), axis=ndim)
        block_fields += ['v1']
        if ndim > 1:
            block = np.concatenate((block, v2[..., np.newaxis]), axis=ndim)
            block_fields += ['v2']
        if ndim > 2:
            block = np.concatenate((block, v3[..., np.newaxis]), axis=ndim)
            block_fields += ['v3']
    # additional check
    assert block.shape[-1] == len(block_fields)
    return block, block_fields


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


class create_amrvac_dict():
    def __init__(self, raw_data, header):
        self._header  = header
        self._standard_fields = copy.deepcopy(header["w_names"])
        self._ndim    = header["ndim"]

        self.data = OrderedDict()

        for var in self._standard_fields:
            idx = self._standard_fields.index(var)
            if self._ndim == 1:
                self.data[var] = raw_data[:, idx]
            elif self._ndim == 2:
                self.data[var] = raw_data[:, :, idx]
            elif self._ndim == 3:
                self.data[var] = raw_data[:, :, :, idx]
            else:
                print("Something wrong with the ndim parameter in the header?")
                print("Current value ndim = {}".format(self._ndim))
                sys.exit(1)


    def add_primitives(self):
        phys_type = self._header["physics_type"]
        gamma = self._header["gamma"]

        if self._ndim == 1:
            if phys_type == "hd":
                v1, p = convert.hd_to_primitive_1d(*(self.data[var] for var in self._standard_fields), gamma=gamma)
            else:
                v1, p = convert.mhd_to_primitive_1d(*(self.data[var] for var in self._standard_fields), gamma=gamma)
            self.data.update({'v1': v1, 'p': p, 'T': p/self.data['rho']})
        elif self._ndim == 2:
            if phys_type == "hd":
                v1, v2, p = convert.hd_to_primitive_2d(*(self.data[var] for var in self._standard_fields), gamma=gamma)
            else:
                v1, v2, p = convert.mhd_to_primitive_2d(*(self.data[var] for var in self._standard_fields), gamma=gamma)
            self.data.update({'v1': v1, 'v2': v2, 'p': p, 'T': p/self.data['rho']})
        else:
            if phys_type == "hd":
                v1, v2, v3, p = convert.hd_to_primitive_3d(*(self.data[var] for var in self._standard_fields), gamma=gamma)
            else:
                v1, v2, v3, p = convert.mhd_to_primitive_3d(*(self.data[var] for var in self._standard_fields), gamma=gamma)
            self.data.update({'v1': v1, 'v2': v2, 'v3': v3, 'p': p, 'T': p/self.data['rho']})