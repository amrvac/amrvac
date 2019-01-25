"""
Created on 5 Dec 2018

@author: Niels Claes

Class to extract all data from header and raw data file.
Creates Instance variables for easy use later on.
Initialize by calling data = ProcessData(dat_file).
"""

from dataIO import dat_reader
from physics import conversions
import numpy as np
import sys, os
import settings
import print_tools

class ProcessData(object):

    def __init__(self, file):
        """
        Initializes Class instance.
        :param file: .dat file, opened in binary mode.
        """ 
        print("Reading %s" % settings.filename)
        hdr = dat_reader.get_header(file)
        
        # Obtain raw data
        try:
            # Mesh is uniformely refined
            raw_data = dat_reader.get_uniform_data(file)
        except IOError:
            print("    Data is not uniformely refined, performing regridding to finest level")
            # Check if data already regridded and saved:
            fn = "dat_files/" + print_tools.trim_filename(settings.filename) + "_regridded_dat.npy"
            if os.path.isfile(fn):
                print("    Regridded data already exists -- loading files.")
                raw_data = np.load(fn)
            else:
                # Data not found, initiate regridding
                if settings.multiple_procs:
                    print("    Regridding data using parallelization, number of procs = %s" % settings.nb_of_procs)
                    raw_data = dat_reader.get_amr_data_multiprocessing(file)
                else:
                    raw_data = dat_reader.get_amr_data(file)

                print("    Regridding done.")
        print("Processing data...")
        
        self._version = hdr["version"]
        self._nw = hdr["nw"]
        self._ndir = hdr["ndir"]
        self._ndim = hdr["ndim"]
        self._levmax = hdr["levmax"]
        self._nleafs = hdr["nleafs"]
        self._it = hdr["it"]
        self._time = hdr["time"]
        self._xmin = hdr["xmin"][0]
        self._xmax = hdr["xmax"][0]
        
        if self._ndim >= 2:
            self._ymin = hdr["xmin"][1]
            self._ymax = hdr["xmax"][1]

        if self._ndim == 3:
            self._zmin = hdr["xmin"][2]
            self._zmax = hdr["xmax"][2]

        self._wnames = hdr["w_names"]
        self._physics_type = hdr["physics_type"]
        if self._physics_type == "hd":
            settings.plot_Blines = False
        self._gamma = hdr["gamma"]
        
        #Initialize primitive and conservative variables.
        #In essence not needed, but done to flag usage before initialization.
        self.rho  = None
        self.momx = None
        self.momy = None
        self.momz = None
        self.e    = None
        self.b1   = None
        self.b2   = None
        self.b3   = None
        self.p    = None
        self.v1   = None
        self.v2   = None
        self.v3   = None
        self.T    = None

        #Ionization degree etc, prevent double loading.
        self.ion    = None
        self.fparam = None
        
        #Define conservative variables
        self._setConservativeVariables(raw_data)
        #Calculate primitive variables
        self._setPrimitiveVariables()
        #Calculate box sizes
        self._setupUniformGrid()
        
        
        print("    Processing done.\n")

        
        
    def _setConservativeVariables(self, raw_data):
        """
        Extracts conservative variables from raw data array and
        defines them as Class Instance variables. Works for 1D, 2D and 3D.
        :param raw_data: Raw data (type np.ndarray) from .dat file
        :raise ValueError: When number of dimensions is not 1, 2 or 3.
        """
        #Extract variables in 1D
        if self._ndim == 1:
            self.rho    = raw_data[:, self._wnames.index("rho")]
            self.momx   = raw_data[:, self._wnames.index("m1")]
            self.e      = raw_data[:, self._wnames.index("e")]
            if self._physics_type == "mhd":
                self.b1 = raw_data[:, self._wnames.index("b1")]
        #Extract variables in 2D
        elif self._ndim == 2:
            self.rho    = raw_data[:, :, self._wnames.index("rho")]
            self.momx   = raw_data[:, :, self._wnames.index("m1")]
            self.momy   = raw_data[:, :, self._wnames.index("m2")]
            self.e      = raw_data[:, :, self._wnames.index("e")]
            if self._physics_type == "mhd":
                self.b1 = raw_data[:, :, self._wnames.index("b1")]
                self.b2 = raw_data[:, :, self._wnames.index("b2")]
        #Extract variables in 3D
        elif self._ndim == 3:
            self.rho    = raw_data[:, :, :, self._wnames.index("rho")]
            self.momx   = raw_data[:, :, :, self._wnames.index("m1")]
            self.momy   = raw_data[:, :, :, self._wnames.index("m2")]
            self.momz   = raw_data[:, :, :, self._wnames.index("m3")]
            self.e      = raw_data[:, :, :, self._wnames.index("e")]
            if self._physics_type == "mhd":
                self.b1 = raw_data[:, :, :, self._wnames.index("b1")]
                self.b2 = raw_data[:, :, :, self._wnames.index("b2")]
                self.b3 = raw_data[:, :, :, self._wnames.index("b3")]
        else:
            raise ValueError("Number of dimensions must be 1, 2 or 3")

    def _setPrimitiveVariables(self):
        """
        Calculates primitive variables from conservative ones.
        Primitives are stored as Instance variables for easy use.
        """
        if self._physics_type == "hd":
            if self._ndim == 1:
                self.p, self.v1 = conversions.hd_to_primitive_1d(self.rho, self.momx, self.e, self._gamma)
            elif self._ndim == 2:
                self.p, self.v1, self.v2 = conversions.hd_to_primitive_2d(self.rho, self.momx, self.momy, self.e, self._gamma)
            else:
                self.p, self.v1, self.v2, self.v3 = conversions.hd_to_primitive_3d(self.rho, self.momx, self.momy, self.momz, self.e, self._gamma)
        elif self._physics_type == "mhd":
            if self._ndim == 1:
                self.p, self.v1 = conversions.mhd_to_primitive_1d(self.rho, self.momx, self.e, self.b1, self._gamma)
            elif self._ndim == 2:
                self.p, self.v1, self.v2 = conversions.mhd_to_primitive_2d(self.rho, self.momx, self.momy, self.e, self.b1, self.b2, self._gamma)
            else:
                self.p, self.v1, self.v2, self.v3 = conversions.mhd_to_primitive_3d(self.rho, self.momx, self.momy, self.momz,
                                                                                    self.e, self.b1, self.b2, self.b3, self._gamma)
        #Set temperature
        self.T = self.p / self.rho
        
    def _setupUniformGrid(self):
        """
        Calculates length and steps of axis arrays, after reduction to uniform grid in dat_reader.py.
        These are instance variables for easy access across modules.
        """
        self._nx = self.rho.shape[0]      
        self._dx = float(self._xmax - self._xmin) / self._nx
        self._x  = np.linspace(self._xmin, self._xmax, self._nx)
        
        if self._ndim >= 2:
            self._ny = self.rho.shape[1]
            self._dy = float(self._ymax - self._ymin) / self._ny
            self._y  = np.linspace(self._ymin, self._ymax, self._ny)
        if self._ndim == 3:
            self._nz = self.rho.shape[2]
            self._dz = float(self._zmax - self._zmin) / self._nz
            self._z  = np.linspace(self._zmin, self._zmax, self._nz)
            
        

        