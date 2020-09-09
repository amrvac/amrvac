"""
Main class to load in MPI-AMRVAC native .dat files.

@author: Niels Claes    (niels.claes@kuleuven.be)
@author: Clément Robert (clement.robert@oca.eu)

         Last edit: 13 October 2019
"""
import sys, os
import numpy as np
import copy
import matplotlib.pyplot as plt

from amrvac_pytools.datfiles.reading import datfile_utilities
from amrvac_pytools.datfiles.processing import regridding, process_data
from amrvac_pytools.datfiles.physics import physical_constants
from amrvac_pytools.datfiles.plotting import amrvac_plotter
from amrvac_pytools.datfiles.views import synthetic


class load_datfile():
    # Following methods provide easy access to functionalities in other scripts
    def amrplot(self, var, **kwargs):
        return amrvac_plotter.amrplot(self, var, **kwargs)
    def rgplot(self, data, **kwargs):
        return amrvac_plotter.rgplot(self, data, **kwargs)
    def halpha(self, **kwargs):
        return synthetic.h_alpha(self, **kwargs)
    def faraday(self, **kwargs):
        return synthetic.faraday(self, **kwargs)
    @staticmethod
    def show():
        plt.show()

    def __init__(self, filename):
        try:
            file = open(filename, "rb")
        except:
            raise FileNotFoundError("Unable to open '{}', is this an MPI-AMRVAC .dat file?".format(filename))

        self._filename = filename
        # Trim filename, assumes format .../.../filenameXXXX.datfiles
        for c in range(len(filename) - 1, -1, -1):
            if filename[c] == '/' or filename[c] == '\\':
                self._filename = filename[c+1:]
                break
        print(">> Reading {}".format(filename))

        # setup basic attributes
        self.file = file
        self.header = datfile_utilities.get_header(file)
        self.data_dict = None
        self._uniform = self._data_is_uniform()
        self._regriddir = "regridded_files"

        # load field information, including primitive variables
        self.known_fields = self._get_known_fields()

        # check if dataset is isothermal
        self.isothermal = self._data_is_isothermal()
        self.adiab_constant = 1

        # load blocktree information
        self.block_lvls, self.block_ixs, self.block_offsets = datfile_utilities.get_tree_info(file)
        self.block_shape = np.append(self.header["block_nx"], self.header["nw"])
        self.domain_width = self.header["xmax"] - self.header["xmin"]

        # setup units
        self.units = physical_constants.units(self.header)

    def get_info(self):
        """
        Prints the current snapshot info (obtained from the header) to the console.
        """
        if self.isothermal:
            info_isothermal = ' (isothermal) '
        else:
            info_isothermal = ''
        print("")
        print("[INFO] Current file      : {}".format(self._filename))
        print("[INFO] Datfile version   : {}".format(self.header["datfile_version"]))
        print("[INFO] Current time      : {}".format(self.header["time"]))
        print("[INFO] Physics type      : {}".format(self.header["physics_type"] + info_isothermal))
        print("[INFO] Boundaries        : {} -> {}".format(self.header["xmin"],
                                                           self.header["xmax"]))
        print("[INFO] Max AMR level     : {}".format(self.header["levmax"]))
        print("[INFO] Uniform mesh      : {}".format(self._uniform))
        print("[INFO] Block size        : {}".format(self.header["block_nx"]))
        print("[INFO] Number of blocks  : {}".format(len(self.block_offsets)))
        print("-" * 40)
        print("Known variables: {}".format(self.known_fields))

    def get_bounds(self):
        """
        Get the bounds of the calculation domain.
        :return: A list containing the bounds in the form [[x1, x2], [y1, y2], [z1, z2]].
                 Depending on the dimension the list size is either 1, 2 or 3.
        """
        bounds = []
        for i in range(len(self.header["xmin"])):
            bounds.append([self.header["xmin"][i], self.header["xmax"][i]])
        return bounds

    def get_coordinate_arrays(self):
        """
        Calculates the domain discretisation.
        :return: A list containing the coordinate axes in the form [x, y, z] where
                 eg. 'x' is a NumPy array representing the domain discretization along the x-axis.
        """
        self._check_datadict_exists()

        w_names = self.header["w_names"]
        xmax = self.header["xmax"]
        xmin = self.header["xmin"]
        nx = np.asarray(self.data_dict[w_names[0]].shape)

        coordinate_arrays = []
        for i in range(len(nx)):
            coordinate_arrays.append(np.linspace(xmin[i], xmax[i], nx[i]))

        return np.asarray(coordinate_arrays, dtype=object)

    def load_all_data(self, nbprocs=None, regriddir=None):
        """
        Loads in all the data to RAM, this can take up quite some space for huge datasets.
        Data will be regridded if this is not already done.
        """
        if self.data_dict is not None:
            return self.data_dict
        if self._uniform:
            data = datfile_utilities.get_uniform_data(self)
        else:
            data = self._regrid_data(nbprocs, regriddir)
        # create dictionary containing the data
        self.data_dict = process_data.create_data_dict(data, self.header, self.adiab_constant)
        return self.data_dict

    def get_time(self):
        """
        Gets the current snapshot time.
        :return: Current time in the simulation (in AMRVAC units).
        """
        return self.header["time"]

    def get_extrema(self, var, print_result=False):
        """
        Returns the maximum and minimum value of the requested variable
        :param var: variable, see ds.known_fields which are available
        :param print_result: if True, prints the extrema to the console
        :return: min, max of the given variable
        """
        if var not in self.known_fields:
            raise KeyError("variable not known: {}".format(var))
        varmax = -1e99
        varmin = 1e99
        for offset in self.block_offsets:
            block = datfile_utilities.get_single_block_data(self.file, offset, self.block_shape)
            block = process_data.create_data_dict(block, self.header)
            varmax = np.maximum(varmax, np.max(block[var]))
            varmin = np.minimum(varmin, np.min(block[var]))
        if print_result:
            print(">> minimum {}: {:2.3e}   |   maximum {}: {:2.3e}".format(var, varmin, var, varmax))
        return varmin, varmax

    def set_adiabatic_constant(self, adiab_cte):
        if not self.isothermal:
            print("[INFO] dataset is not isothermal, 'e' is present. Setting the adiabatic constant has no effect.")
        self.adiab_constant = adiab_cte




    # 'protected' methods
    def _data_is_uniform(self):
        """
        Checks if the data is uniformely refined.
        :return: True if grid is uniform, False otherwise.
        """
        refined_nx = 2 ** (self.header['levmax'] - 1) * self.header['domain_nx']
        nleafs_uniform = np.prod(refined_nx / self.header['block_nx'])
        if not self.header["nleafs"] == nleafs_uniform:
            return False
        return True

    def _data_is_isothermal(self):
        """
        Checks if the data is isothermal
        :return: True if hd or mhd but no energy present
        """
        if (not 'e' in self.known_fields) and \
           (self.header['physics_type'] == 'hd' or self.header['physics_type'] == 'mhd'):
            return True
        else:
            return False

    def _check_regrid_directory(self, regriddir):
        """
        Checks if the specified save directory is present. Defaults to 'regridded_files', if this folder is not
        present it is created. If 'regriddir' is explicitly specified and the folder is not present, throws an error
        and terminates the program (to prevent creation of unwanted folders).
        :param regriddir: The directory to which to save the regridded files (and load them back in).
        """
        if regriddir is None:
            regriddir = self._regriddir
            if not os.path.isdir(regriddir):
                os.mkdir(regriddir)
                print("[INFO] Created directory: {}".format(regriddir))
        else:
            if not os.path.isdir(regriddir):
                sys.exit("Specified directory does not exist: {}".format(regriddir))
            self._regriddir = regriddir

    def _regrid_data(self, nbprocs, regriddir):
        """
        Regrids non-uniform data, which is saved as a Numpy file in the folder 'regridded_files'
        (created if not present). The regrid directory is searched first if the data is already present.
        If so it is loaded, if not regridding is started.
        :param nbprocs: Number of processors to do the regridding. Defaults to max number - 2
        :param regriddir: The directory to which the regridded file is saved. This directory is searched in order to
                          try and load the file if already present.
        """
        self._check_regrid_directory(regriddir)
        try:
            data = self._load_regridded_data()
        except FileNotFoundError:
            data = regridding.regrid_amr_data(self, nbprocs)
            self._save_regridded_data(data)
        return data

    def _save_regridded_data(self, data):
        """
        Saves the regridded data as a NumPy file.
        :param data: The raw, regridded data.
        """
        filename_regridded = self._regriddir + "/" + self._filename[:-4] + "_regridded.npy"
        np.save(filename_regridded, data)
        print("[INFO] Regridded data saved to {}".format(filename_regridded))

    def _load_regridded_data(self):
        """
        Loads the regridded data file back in, if present.
        :return: The data file if it is present.
        :exception: FileNotFoundError if the file is not present. The dataset will be regridded in this case.
        """
        filename_regridded = self._regriddir + "/" + self._filename[:-4] + "_regridded.npy"

        if os.path.isfile(filename_regridded):
            try:
                data = np.load(filename_regridded)
                print("[INFO] Regridded file found and loaded ({}).".format(filename_regridded))
                return data
            except OSError:
                print("[ERROR] File found but failed to load: {}".format(filename_regridded))
                sys.exit(1)
        else:
            raise FileNotFoundError

    def _check_datadict_exists(self):
        """
        Checks if the data dictionary is loaded.
        :raises AttributeError  if self.data_dict is None
        """
        if self.data_dict is None:
            print("[INFO] Dataset must be loaded to do this, call load_all_data() first.")
            raise AttributeError

    def _get_known_fields(self):
        """
        Retrieves the known fields for a given dataset.
        :return: Regular list with the known fields as strings
        """
        fields = copy.deepcopy(self.header['w_names'])
        if self.header['physics_type'] == 'hd' or self.header['physics_type'] == 'mhd':
            if self.header['ndim'] == 1:
                fields += ['v1', 'p']
            elif self.header['ndim'] == 2:
                fields += ['v1', 'v2', 'p']
            else:
                fields += ['v1', 'v2', 'v3', 'p']
            if 'e' in fields:
                fields += ['T']
        return fields
