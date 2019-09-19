import sys
import numpy as np
import matplotlib.colors
import matplotlib.pyplot as plt

from amrvac_tools.datfiles.reading import datfile_utilities
from amrvac_tools.datfiles.processing import regridding, process_data
from amrvac_tools.datfiles.physics import ionisation


class _syntheticmain():
    def __init__(self, dataset, **kwargs):
        self.dataset = dataset
        self.dataset.units.check_default_units()
        if self.dataset.header['ndim'] == 1:
            print("synthetic views can only be created for 2D/3D data")
            sys.exit(1)

        # initialise figure and axis
        fig = kwargs.get("fig", None)
        ax = kwargs.get("ax", None)
        if fig is None or ax is None:
            self.fig, self.ax = plt.subplots(1)
        else:
            self.fig = fig
            self.ax = ax
        # these are initialised in the subclasses
        self.cmap = None
        self.logscale = None
        self.colorbar = None

        # initialise variables
        self.line_of_sight = kwargs.get("line_of_sight", "x")
        self.altitude = kwargs.get("altitude", 20000)
        self.simulation_type = kwargs.get("simulation_type", "prominence")
        self.f23 = 0.6407  # oscillator strength of H-alpha line
        # initialise splines for interpolation
        ionisation.init_splines(self.altitude)
        self.integrated_block_list = []
        self.block_nx = self.dataset.header['block_nx']
        self.block_nx_int = self._reduce_list_to_2d(self.block_nx)

    def _get_ne(self, block, block_fields, block_ion):
        block_p = block[..., block_fields.index("p")] * self.dataset.units.unit_pressure
        block_T = block[..., block_fields.index("T")] * self.dataset.units.unit_temperature
        block_ne = block_p / ((1 + 1.1 / block_ion) * self.dataset.units.k_B * block_T)
        return block_ne

    def _integrate_block(self, block, l_edge, r_edge):
        """
        Integrates a given block along the given line of sight.
        :param block: the block to integrate
        :param l_edge: contains the left edge of the block, as a ndim length list [x0(, y0, z0)]
        :param r_edge: contains the right edge of the block, as an ndim length list [x1(, y1, z1)]
        :return: 2D numpy array, containing the integrated block. If the block is originally 2D the block itself
                 is returned.
        """
        if self.dataset.header["ndim"] == 2:
            return block

        if self.line_of_sight == 'x':
            x = np.linspace(l_edge[0], r_edge[0], self.block_nx[0])
            result = np.zeros_like(block[0, :, :])
            for i, j in np.ndindex(result.shape):
                col = block[:, i, j]
                integrated_col = np.trapz(col, x)
                result[i, j] = integrated_col
        elif self.line_of_sight == 'y':
            y = np.linspace(l_edge[1], r_edge[1], self.block_nx[1])
            result = np.zeros_like(block[:, 0, :])
            for i, j in np.ndindex(result.shape):
                col = block[i, :, j]
                integrated_col = np.trapz(col, y)
                result[i, j] = integrated_col
        else:
            z = np.linspace(l_edge[2], r_edge[2], self.block_nx[2])
            result = np.zeros_like(block[:, :, 0])
            for i, j in np.ndindex(result.shape):
                col = block[i, j, :]
                integrated_col = np.trapz(col, z)
                result[i, j] = integrated_col
        return result

    def _reduce_list_to_2d(self, list_in):
        if self.dataset.header['ndim'] == 2:
            return np.asarray(list_in)

        if self.line_of_sight == 'x':
            array2d = np.asarray([list_in[1], list_in[2]])
        elif self.line_of_sight == 'y':
            array2d = np.asarray([list_in[0], list_in[2]])
        else:
            array2d = np.asarray(list_in[:-1])
        return array2d

    def _merge_integrated_blocks(self):
        self.integrated_block_list = np.asarray(self.integrated_block_list)

        # Initialise merged 2D matrix
        maxlvl = np.max(self.dataset.block_lvls)
        refined_nx = 2**(maxlvl - 1) * self.dataset.header['domain_nx']
        refined_nx = self._reduce_list_to_2d(refined_nx)
        merged_result = np.zeros(tuple(refined_nx))

        # NOTE: one can also add all the blocks at a specific level together, and THEN regrid the resulting 2D
        #       matrix to the maximum level. This is faster, however, because the 2D matrix still contains zeros
        #       from the blocks not at the level at that moment considered, these are also accounted for during
        #       interpolation. This has unintended consequences near the block edges between two different levels,
        #       such as clearly visible grid lines and very small values not equal to zero (where they should be zero).

        # iterate over blocks and corresponding levels in dataset
        for ileaf, blocklvl in enumerate(self.dataset.block_lvls):
            grid_power_diff = 2 ** (maxlvl - blocklvl)          # power of 2 difference between level and max level
            regrid_width = self.block_nx_int * grid_power_diff  # index width of the block at max resolution
            # regrid current block to max level
            block = regridding.regrid_2dmatrix(self.integrated_block_list[ileaf], tuple(regrid_width))
            # retrieve current index in morton curve, reduce to 2D (due to integration)
            block_morton_idx = self._reduce_list_to_2d(self.dataset.block_ixs[ileaf])
            # calculate left and right index positions of the block using morton index
            idx0 = (block_morton_idx - 1) * grid_power_diff * self.block_nx_int
            idx1 = idx0 + tuple(regrid_width)
            # add this block to the integrated result
            merged_result[idx0[0]:idx1[0], idx0[1]:idx1[1]] += block
        return merged_result

    def _plot_synthetic_view(self, view):
        bounds_x, bounds_y = self._reduce_list_to_2d(self.dataset.get_bounds())
        norm = None
        if self.logscale:
            norm = matplotlib.colors.LogNorm()
        im = self.ax.imshow(np.rot90(view), extent=[*bounds_x, *bounds_y], cmap=self.cmap, norm=norm)
        self.colorbar = self.fig.colorbar(im)
        self.ax.set_title("integrated over {}".format(self.line_of_sight))


class h_alpha(_syntheticmain):
    def __init__(self, dataset, **kwargs):
        print(">> Creating H-alpha view...")

        super().__init__(dataset, **kwargs)
        self.cmap = kwargs.get("cmap", "Reds_r")
        self.logscale = kwargs.get("logscale", True)
        self._calculate_halpha_view()

    def _calculate_halpha_view(self):
        for ileaf, offset in enumerate(self.dataset.block_offsets):
            block = datfile_utilities.get_single_block_data(self.dataset.file, offset, self.dataset.block_shape)
            # this adds the temperature and pressure to the block
            block, block_fields = process_data.add_primitives_to_single_block(block, self.dataset)
            # interpolate ionisation and f parameter for each block
            block_ion, block_fpar = ionisation.block_interpolate_ionisation_f(block, block_fields, self.dataset)
            block_ne = super()._get_ne(block, block_fields, block_ion)
            n2 = block_ne**2 / (block_fpar * 1e16)              # parameter f is interpolated in units of 1e16 cm-3
            # calculate block opacity
            block_kappa = (np.pi * self.dataset.units.ec**2 / (self.dataset.units.m_e * self.dataset.units.c)) * \
                            self.f23 * n2 * self._gaussian(block, block_fields)
            # integrate block along line of sight to get opacity
            l_edge, r_edge = process_data.get_block_edges(ileaf, self.dataset)
            opacity = super()._integrate_block(block_kappa, l_edge, r_edge)

            S = self._source_function()
            intensity = S * (1 - np.exp(-opacity))
            if self.simulation_type == 'filament':
                Ibgr = 2.2 * S
                intensity += Ibgr * np.exp(-opacity)
            # add integrated block to list, this preserves block index
            self.integrated_block_list.append(intensity)

        # merge all integrated blocks into one single 2D array
        view = super()._merge_integrated_blocks()
        # plot final result
        super()._plot_synthetic_view(view)

    def _gaussian(self, block, block_fields):
        block_T = block[..., block_fields.index("T")] * self.dataset.units.unit_temperature
        ksi = 5 * 1e5  # microturbulence in cm/s
        nu_0 = self.dataset.units.c / (6562.8 * 1e-8)       # H-alpha wavelength is 6562.8 Angstrom
        delta_nu = 0
        delta_nuD = (nu_0 / self.dataset.units.c) * \
                    np.sqrt(2 * self.dataset.units.k_B * block_T / self.dataset.units.m_p + ksi ** 2)
        phi_nu = 1.0 / (np.sqrt(np.pi) * delta_nuD) * np.exp(-delta_nu / delta_nuD) ** 2
        return phi_nu

    def _source_function(self):
        H = self.altitude * 1e5     # altitude in cm
        # dilution factor W
        W = 0.5 * ( 1 - np.sqrt(1 - (self.dataset.units.Rsun**2 / (self.dataset.units.Rsun + H)**2)) )
        return W * 0.17 * 4.077 * 1e-5




class faraday(_syntheticmain):
    def __init__(self, dataset, **kwargs):
        print(">> Creating Faraday view...")
        if not dataset.header["physics_type"] == "mhd":
            print("calculating Faraday rotation measure is only possible for an MHD dataset")
            return
        if not dataset.header["ndim"] == 3:
            print("calculating Faraday rotation measure is only possible for a 3D dataset")
            return

        super().__init__(dataset, **kwargs)
        self.cmap = kwargs.get("cmap", "jet")
        self.logscale = kwargs.get("logscale", False)
        self._calculate_faraday_view()

    def _calculate_faraday_view(self):
        for ileaf, offset in enumerate(self.dataset.block_offsets):
            block = datfile_utilities.get_single_block_data(self.dataset.file, offset, self.dataset.block_shape)
            # add pressure and temperature to block
            block, block_fields = process_data.add_primitives_to_single_block(block, self.dataset)
            block_ion, block_fpar = ionisation.block_interpolate_ionisation_f(block, block_fields, self.dataset)
            block_ne = super()._get_ne(block, block_fields, block_ion)
            prefactor = self.dataset.units.ec**3 / (2*np.pi*self.dataset.units.m_e**2 * self.dataset.units.c**4)
            if self.line_of_sight == 'x':
                b_idx = 'b1'
            elif self.line_of_sight == 'y':
                b_idx = 'b2'
            else:
                b_idx = 'b3'
            b_para = block[..., block_fields.index(b_idx)] * self.dataset.units.unit_magneticfield
            l_edge, r_edge = process_data.get_block_edges(ileaf, self.dataset)
            fara_measure = super()._integrate_block(block_ne*b_para, l_edge, r_edge) * prefactor
            self.integrated_block_list.append(fara_measure)
        view = super()._merge_integrated_blocks()
        super()._plot_synthetic_view(view)
        self.colorbar.set_label("rad cm$^{-2}$")

