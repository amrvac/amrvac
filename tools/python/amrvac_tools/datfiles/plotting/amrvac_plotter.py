import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

from amrvac_tools.datfiles.reading import datfile_utilities
from amrvac_tools.datfiles.processing import process_data


class _plotsetup():
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

class amrplot(_plotsetup):
    def __init__(self, dataset, var, **kwargs):
        super().__init__(dataset, **kwargs)
        self.var = var

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
        for ileaf, offset in enumerate(self.dataset.block_offsets):
            l_edge, r_edge = process_data.get_block_edges(ileaf, self.dataset)
            block = datfile_utilities.get_single_block_data(self.dataset.file, offset, self.dataset.block_shape)
            block, block_fields = process_data.add_primitives_to_single_block(block, self.dataset, add_velocities=True)
            varidx = block_fields.index(self.var)
            block_data = block[:, varidx]
            x = np.linspace(l_edge, r_edge, self.dataset.header['block_nx'])
            self.ax.plot(x, block_data, '-k')

    def plot_2d(self):
        varmin, varmax = self.dataset.get_extrema(self.var)
        norm = None
        if self.logscale:
            norm = matplotlib.colors.LogNorm()
        # iterate over blocks in dataset
        for ileaf, offset in enumerate(self.dataset.block_offsets):
            # retrieve x and y coordinates for each block
            l_edge, r_edge = process_data.get_block_edges(ileaf, self.dataset)
            # read in block data (contains all variables)
            block = datfile_utilities.get_single_block_data(self.dataset.file, offset, self.dataset.block_shape)
            # add primitive variables to this block
            block, block_fields = process_data.add_primitives_to_single_block(block, self.dataset, add_velocities=True)
            # get index of variable to plot
            varidx = block_fields.index(self.var)
            block_data = block[:, :, varidx]
            x = np.linspace(l_edge[0], r_edge[0], self.dataset.header['block_nx'][0])
            y = np.linspace(l_edge[1], r_edge[1], self.dataset.header['block_nx'][1])
            im = self.ax.pcolormesh(x, y, block_data.T, cmap=self.cmap, vmin=varmin, vmax=varmax, norm=norm)

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
        self.fig.colorbar(im, cax=cax)
        self.fig.tight_layout()


class rgplot(_plotsetup):
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
        bounds_x, bounds_y = self.dataset.get_bounds()
        norm = None
        if self.logscale:
            norm = matplotlib.colors.LogNorm()
        im = self.ax.imshow(np.rot90(self.data), extent=[*bounds_x, *bounds_y], cmap=self.cmap, norm=norm)
        self.ax.set_aspect('equal')
        divider = make_axes_locatable(self.ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        self.fig.colorbar(im, cax=cax)
        self.fig.tight_layout()


