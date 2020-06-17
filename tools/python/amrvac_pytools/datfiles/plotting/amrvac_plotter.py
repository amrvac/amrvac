import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

from amrvac_pytools.datfiles.reading import datfile_utilities
from amrvac_pytools.datfiles.processing import process_data


class _plotsetup:
    """
    Parent class of amrplot and rgplot, contains the matplotlib figure initialisations and colormaps.
    """
    def __init__(self, dataset, **kwargs):
        self.dataset = dataset

        self.cmap = kwargs.get("cmap", "jet")
        self.logscale = kwargs.get("logscale", False)

        # initialise figure and axis
        fig = kwargs.get("fig", None)
        ax = kwargs.get("ax", None)
        if fig is None or ax is None:
            self.fig, self.ax = plt.subplots(1)
        else:
            self.fig = fig
            self.ax = ax
        self.colorbar = None
        self.varmin = kwargs.get('varmin', None)
        self.varmax = kwargs.get('varmax', None)

class amrplot(_plotsetup):
    """
    Plots AMR data for a variable in the known fields.
    """
    def __init__(self, dataset, var, **kwargs):
        super().__init__(dataset, **kwargs)
        self.var = var
        if not isinstance(var, str):
            raise TypeError("'amrplot' takes the variable name as argument")
        if self.var not in dataset.known_fields:
            raise KeyError("Variable '{}' not in the list of known fields".format(self.var))

        # mesh-related parameters
        self.draw_mesh = kwargs.get("draw_mesh", False)
        if not isinstance(self.draw_mesh, bool):
            raise ValueError("'draw_mesh' argument should be True or False")
        self.mesh_color = kwargs.get("mesh_color", "black")
        self.mesh_linestyle = kwargs.get("mesh_linestyle", "solid")
        self.mesh_linewidth = kwargs.get("mesh_linewidth", 2)
        self.mesh_opacity = kwargs.get("mesh_opacity", 1)

        if self.dataset.header["ndim"] == 1:
            self.plot_1d()
        elif self.dataset.header["ndim"] == 2:
            self.plot_2d()
        else:
            raise NotImplementedError("Plotting in 3D is not supported")

    def plot_1d(self):
        if self.varmin is None or self.varmax is None:
            self.varmin, self.varmax = self.dataset.get_extrema(self.var)

        for ileaf, offset in enumerate(self.dataset.block_offsets):
            l_edge, r_edge = process_data.get_block_edges(ileaf, self.dataset)
            block = datfile_utilities.get_single_block_data(self.dataset.file, offset, self.dataset.block_shape)
            block = process_data.create_data_dict(block, self.dataset.header)
            x = np.linspace(l_edge, r_edge, self.dataset.header['block_nx'][0])
            self.ax.plot(x, block[self.var], '-k')
            self.ax.set_title(self.var)
            self.ax.set_ylim([self.varmin,self.varmax])

    def plot_2d(self):
        if self.varmin is None or self.varmax is None:
            self.varmin, self.varmax = self.dataset.get_extrema(self.var)

        norm = None
        if self.logscale:
            norm = matplotlib.colors.LogNorm()
        # iterate over blocks in dataset
        for ileaf, offset in enumerate(self.dataset.block_offsets):
            # retrieve x and y coordinates for each block
            l_edge, r_edge = process_data.get_block_edges(ileaf, self.dataset)
            # read in block data (contains all variables)
            block = datfile_utilities.get_single_block_data(self.dataset.file, offset, self.dataset.block_shape)
            block = process_data.create_data_dict(block, self.dataset.header)
            x = np.linspace(l_edge[0], r_edge[0], self.dataset.header['block_nx'][0])
            y = np.linspace(l_edge[1], r_edge[1], self.dataset.header['block_nx'][1])
            im = self.ax.pcolormesh(x, y, block[self.var].T, cmap=self.cmap, vmin=self.varmin, vmax=self.varmax, norm=norm)

            # logic to draw the mesh
            if self.draw_mesh:
                if not r_edge[0] == self.dataset.header["xmax"][0]:
                    self.ax.vlines(x=r_edge[0], ymin=l_edge[1], ymax=r_edge[1], color=self.mesh_color,
                                   lw=self.mesh_linewidth, linestyle=self.mesh_linestyle, alpha=self.mesh_opacity)
                if not r_edge[1] == self.dataset.header["xmax"][1]:
                    self.ax.hlines(y=r_edge[1], xmin=l_edge[0], xmax=r_edge[0], color=self.mesh_color,
                                   lw=self.mesh_linewidth, linestyle=self.mesh_linestyle, alpha=self.mesh_opacity)
        self.ax.set_aspect('equal')
        divider = make_axes_locatable(self.ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        self.colorbar = self.fig.colorbar(im, cax=cax)
        self.ax.set_title(self.var)
        self.fig.tight_layout()


class rgplot(_plotsetup):
    """
    Plots a variable from a regridded dataset, requires that load_all_data() is called first.
    """
    def __init__(self, dataset, data, **kwargs):
        if dataset.data_dict is None:
            raise AttributeError("Make sure the regridded data is loaded when calling rgplot (ds.load_all_data)")
        if not isinstance(data, np.ndarray):
            raise Exception("Attribute 'data' passed should be a numpy array containing the data. "
                            "data = ds.load_all_data(), followed by eg. rgplot(ds, data['rho']")
        super().__init__(dataset, **kwargs)

        self.data = data
        if len(self.data.shape) == 1:
            self.plot_1d()
        elif len(self.data.shape) == 2:
            self.plot_2d()
        else:
            raise NotImplementedError("Plotting in 3D is not supported")


    def plot_1d(self):
        x = self.dataset.get_coordinate_arrays()[0]
        self.ax.plot(x, self.data, '-k')

    def plot_2d(self):
        if self.varmin is None or self.varmax is None:
            self.varmin, self.varmax = np.min(self.data), np.max(self.data)

        bounds_x, bounds_y = self.dataset.get_bounds()
        norm = None
        if self.logscale:
            norm = matplotlib.colors.LogNorm()
        im = self.ax.imshow(np.rot90(self.data), extent=[*bounds_x, *bounds_y], cmap=self.cmap, norm=norm,
                            vmin=self.varmin, vmax=self.varmax)
        self.ax.set_aspect('equal')
        divider = make_axes_locatable(self.ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        self.colorbar = self.fig.colorbar(im, cax=cax)
        self.fig.tight_layout()


