#!/usr/bin/env python
# coding: utf-8

# # Plot a single variable for chosen time stamps

# In[8]:


import yt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from pathlib import Path
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import distance_transform_edt
import matplotlib.colors as mcolors
from matplotlib.ticker import FuncFormatter
import argparse

parser = argparse.ArgumentParser(description="Plot UAWSoM solar-atmosphere output")
parser.add_argument("data_dir", nargs="?", default=".",
                    help="directory containing UAWSoM .dat files")
parser.add_argument("--prefix", default="uawsom_solar_",
                    help="snapshot filename prefix")
parser.add_argument("--snapshots", nargs="+", type=int,
                    default=[0, 5, 10, 15, 20])
args = parser.parse_args()

# [set what time stamps you want to show]
time_stamps = args.snapshots

units = dict(length_unit=(1e9, 'cm'),
             temperature_unit=(1e6, 'K'),
             numberdensity_unit=(1e9, 'cm**-3'))

path1 = Path(args.data_dir)

dataset = [
    yt.load(path1 / f'{args.prefix}{str(n).zfill(4)}.dat',
            units_override=units, unit_system='cgs')
    for n in time_stamps
]

# See what variables are available to plot

print("Available fields:")
print(dataset[0].field_list)

# Set the variable you want to plot

#available fields: [('amrvac', 'QAW'), ('amrvac', 'Qk'), ('amrvac', 'Te'), ('amrvac', 'b1'),
                #   ('amrvac', 'b1t'), ('amrvac', 'b2'), ('amrvac', 'b2t'), ('amrvac', 'b3'),
                #   ('amrvac', 'b3t'), ('amrvac', 'e'), ('amrvac', 'm1'), ('amrvac', 'm2'),
                #   ('amrvac', 'm3'), ('amrvac', 'rad'), ('amrvac', 'rho'), ('amrvac', 'wAminus'),
                #   ('amrvac', 'wAplus'), ('amrvac', 'wkminus'), ('amrvac', 'wkplus')]

# also can plot various 'gas' quantities in dimensional units

field_type = 'gas' # amrvac or gas only
field_name = 'number_density' # must match field
unit_value = r'n_H' # arbitrary name
cmap = 'viridis' # colour scheme for plotting


fontsize = 15
labelsize = 12

fig = plt.figure()

grid = AxesGrid(
    fig, (0.04, 0.05, 0.90, 0.94),
    axes_pad=0.,
    share_all=True,
    nrows_ncols=(1, len(time_stamps)),
    cbar_location="right",
    label_mode="L",     # axes get left side labels
    cbar_mode="single",
    cbar_size="5%",
    cbar_pad="5%",
)

for i, ds in enumerate(dataset):

    time = round(float(ds.current_time) * float(ds.time_unit) / 60) # 60 converts to minutes

    p = yt.plot_2d(ds, (field_type, field_name),
                   origin='native', fontsize=fontsize)

    p.set_axes_unit("Mm")
    p.set_cmap((field_type, field_name), cmap)
    p.set_colorbar_label((field_type, field_name), unit_value)
    p.set_log((field_type, field_name), True) # Plot the data log scaled or not

    # It might be useful to change the limits of data plotted so the chromosphere doesn't dominate the contrast
    p.set_zlim((field_type, field_name), 1e8, 1e11)

    # Sometimes you may wish to have different colour labels if the cmap colour contrast doesn't work nicely
    #if time in [1, 37, 74]:
    #    label_color = "black"
    #else:
    #    label_color = "white"

    # Can also set all labels to be the same colour
    label_color = "white"

    # Location of label in the top left of each panel
    p.annotate_text(
        (-45, 68),
        f"{time} min",
        coord_system='plot',
        text_args={"color": label_color}
    )

    # Attach to AxesGrid
    plot = p.plots[(field_type, field_name)]
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]
    plot.cax.tick_params(labelsize=labelsize)

    # yt finalizes the plot (this overwrites ticks)
    p._setup_plots()

    # --- NOW override ticks (yt won't overwrite them anymore) ---
    # These relate to the dimensions in Mm, your box is -50 to 50 Mm in x and 0 to 80 Mm in y (vertical direction)
    ax = grid[i].axes
    ax.set_xticks([-40, -20, 0, 20, 40])
    ax.set_xticklabels(['-40', '-20', '0', '20', '40'])

    ax.set_yticks([0, 20, 40, 60, 80])
    ax.set_yticklabels(['0', '20', '40', '60', '80'])

# Resize figure
fig.set_size_inches(14, 6)
#plt.show()

plt.savefig(path1 / 'fig1.png')

# Save figure
#fig.savefig(
#    "path_to_folder/plot.png",
#    dpi=300,
#    bbox_inches="tight"
#)


# # Plot a grid of variables at given time stamps

# In[2]:


# amrvac will save the conserved variables but you can derive new variables such as the velocity
# Here I use 116.4508 to convert to km/s. This is taken from unit_velocity = 11645084.295622544
# which you will see in your terminal, this value refers to cgs units and hence cm/s

def _v1(field, data):
    return data.ds.arr(
        116.4508*data[("amrvac", "m1")]/data[("amrvac", "rho")],
        "km/s"
    )

def _v2(field, data):
    return data.ds.arr(
        116.4508*data[("amrvac", "m2")]/data[("amrvac", "rho")],
        "km/s"
    )

yt.add_field(
    ("amrvac", "v1"),
    function=_v1,
    units="km/s",
    sampling_type="cell",
    force_override=True
)

yt.add_field(
    ("amrvac", "v2"),
    function=_v2,
    units="km/s",
    sampling_type="cell",
    force_override=True
)

def _vmag(field, data):
    return data.ds.arr(
        (data[("amrvac", "v2")]**2 + data[("amrvac", "v1")]**2)**0.5,
        "km/s"
    )

yt.add_field(
    ("amrvac", "vmag"),
    function=_vmag,
    units="km/s",
    sampling_type="cell",
    force_override=True
)

dataset = [
    yt.load(path1 / f'{args.prefix}{str(n).zfill(4)}.dat',
            units_override=units, unit_system='cgs')
    for n in time_stamps
]

# See what variables are available to plot

print("Available fields:")
print(dataset[0].field_list)

variables = [
    ("amrvac", "Te", r"T (MK)", "hot"),
    ("gas", "number_density", r"$n_{\mathrm{H}}$ (cm$^{-3}$)", "viridis"),
    ("amrvac", "v1", r"$v_{\mathrm{x}}$ (kms$^{-1}$)", "bwr"),
    ("amrvac", "v2", r"$v_{\mathrm{y}}$ (kms$^{-1}$)", "bwr")
]

fontsize = 15
labelsize = 12

fig = plt.figure(figsize=(16, 24))

grid = AxesGrid(
    fig, (0.05, 0.05, 0.90, 0.90),
    nrows_ncols=(len(variables), len(time_stamps)),
    axes_pad=0.0,
    share_all=False,
    cbar_location="right",
    cbar_mode="edge",      # separate colorbar per row
    cbar_size="5%",
    cbar_pad="5%",
    label_mode="L"
)
# Loop over rows (variables) and columns (snapshots)
for row, (field_type, field_name, unit_value, cmap) in enumerate(variables):

    for col, ds in enumerate(dataset):

        idx = row * len(time_stamps) + col   # grid index

        time = round(float(ds.current_time) * float(ds.time_unit) / 60)

        p = yt.plot_2d(ds, (field_type, field_name),
                       origin='native', fontsize=fontsize)

        p.set_axes_unit("Mm")
        p.set_cmap((field_type, field_name), cmap)
        p.set_colorbar_label((field_type, field_name), unit_value)

        # Annotate time
        if row == 1:
            label_color = "white"
        else:
            label_color = "black"

        p.annotate_text(
            (-45, 68),
            f"{time} min",
            coord_system='plot',
            text_args={"color": label_color}
        )
        #if row == 2:   # bottom row
        #    p.set_zlim((field_type, field_name), 1e-4, 1)
        #if row == 1:
        #    p.set_zlim((field_type, field_name), 1.1e8, 0.9e12)

        if row == 0:
            p.set_log((field_type, field_name), False)

        if row == 1:
            p.set_log((field_type, field_name), True)

        if row == 2:
            p.set_log((field_type, field_name), False)

        if row == 3:
            p.set_log((field_type, field_name), False)

        # Attach yt plot to AxesGrid
        plot = p.plots[(field_type, field_name)]
        plot.figure = fig
        plot.axes = grid[idx].axes
        plot.cax = grid[idx].cax
        plot.cax.tick_params(labelsize=labelsize)

        p._setup_plots()

        # Optional: custom ticks
        ax = grid[idx].axes
        ax.set_xticks([-40, -20, 0, 20, 40])
        ax.set_xticklabels(['-40', '-20', '0', '20', '40'])
        ax.set_yticks([0, 20, 40, 60, 80])
        ax.set_yticklabels(['0', '20', '40', '60', '80'])

fig.set_size_inches(14, 6)
#plt.show()
plt.savefig(path1 / 'fig2.png')

# # Print the minimum and maximum physical velocities

# In[9]:


for i in range(len(time_stamps)):
    ds = dataset[i]

    ad = ds.all_data()

    # Get coordinates
    y = ad[("index", "y")].to("Mm")

    # Select only y > 5 Mm
    mask = y > 0 #Change if you want to remove chromosphere

    vmag = ad[("amrvac", "vmag")]

    print("vmag above mask at time " + str(round(time_stamps[i]*85.87/60, 1)) + " mins :")
    print("min =", vmag[mask].min())
    print("max =", vmag[mask].max())
    print("---------------------------------")


    v1 = ad[("amrvac", "v1")]

    print("v1 above mask:")
    print("min =", v1[mask].min())
    print("max =", v1[mask].max())
    print("---------------------------------")

    v2 = ad[("amrvac", "v2")]

    print("v2 above mask:")
    print("min =", v2[mask].min())
    print("max =", v2[mask].max())
    print("---------------------------------")

    wkminus = ad[("amrvac", "wkminus")]

    print("wkminus above mask:")
    print("min =", wkminus[mask].min())
    print("max =", wkminus[mask].max())
    print("---------------------------------")


# # Visualise the data over a single field line

# In[4]:


# Change 10 to the time step you want to plot
ds = yt.load(path1 / f"{args.prefix}0010.dat")

footpoint = -4.0  # This means the footpoint of the field line is at x=-40Mm
footpoint = -4.0 # 40Mm will run faster than e.g., 47.5 Mm as the loop is much shorter

var_name = "Qk" # change this if you want to plot another variable
var_scale = 0.003697887481489527  # conversion factor to domensional heating rate
                                  # note: different variables have different conversion factors

level = ds.index.max_level
dims = ds.domain_dimensions * 2**level # if you turn on amr '*2**level' is important

cg = ds.covering_grid(
    level=level,
    left_edge=ds.domain_left_edge,
    dims=dims
)

x = cg["x"].to("code_length").v[:, 0, 0]
y = cg["y"].to("code_length").v[0, :, 0]

Bx = cg[("amrvac", "b1t")].v[:, :, 0]
By = cg[("amrvac", "b2t")].v[:, :, 0]
var = cg[("amrvac", var_name)].v[:, :, 0]

xmin, xmax = float(x.min()), float(x.max())
ymin, ymax = float(y.min()), float(y.max())

# Interpolators
Bx_i = RegularGridInterpolator((x, y), Bx)
By_i = RegularGridInterpolator((x, y), By)
var_i = RegularGridInterpolator((x, y), var)

def rk4_step(pos, ds_step):
    def B(p):
        return np.array([Bx_i(p).item(), By_i(p).item()])

    k1 = B(pos)
    k2 = B(pos + 0.5 * ds_step * k1)
    k3 = B(pos + 0.5 * ds_step * k2)
    k4 = B(pos + ds_step * k3)
    return pos + (ds_step / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

def trace_fieldline(seed, ds_step=1e-4, nsteps=200000,
                    y_target=None, tol=1e-3, min_dist=0.05): # if code takes too long reduce ds_step/tol/min_dist

    if y_target is None:
        y_target = ymin

    pts = [np.array(seed, dtype=float)]
    p = np.array(seed, dtype=float)
    arc = 0.0

    for _ in range(nsteps):
        p_new = rk4_step(p, ds_step)

        if not (xmin < p_new[0] < xmax):
            break
        if not (ymin <= p_new[1] <= ymax):
            break

        arc += np.linalg.norm(p_new - p)
        pts.append(p_new.copy())
        p = p_new

        if arc > min_dist and abs(p[1] - y_target) < tol:
            break

    return np.array(pts)

# Trace field line
seed = np.array([footpoint, ymin])
fl = trace_fieldline(seed)

# Sample heating along the line
var_line = var_i(fl)

# Arc length
ds_vals = np.sqrt(np.sum(np.diff(fl, axis=0)**2, axis=1))
s = np.concatenate([[0.0], np.cumsum(ds_vals)])

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)

# Left panel: profile along field line
ax1.plot(s * 10.0, var_scale * var_line, linewidth=2) #*10 to convert to Mm
ax1.set_xlabel("Arc length s (Mm)")
ax1.set_ylabel(r"$Q_{k}$ (erg cm$^{-3}$ s$^{-1}$)")
ax1.set_title("Kink wave heating along field line")
ax1.grid(True)
ax1.set_yscale("log") # you want this on for most variables other than velocity/momentum

# Right panel: 2D field with field line overlay
pcm = ax2.pcolormesh(
    x, y, var_scale*var.T,
    shading="auto",
    cmap="hot",
    norm=mcolors.LogNorm(vmin=np.nanmin(var_scale*var[var > 0]), vmax=np.nanmax(var_scale*var))
)

ax2.plot(
    fl[:, 0], fl[:, 1],
    color="black",
    linewidth=2,
    label=f'Footpoints = ±{-10 * seed[0]:.2f} Mm'
)

ax2.set_xlabel("x (Mm)", fontsize=16)
ax2.set_ylabel("y (Mm)", fontsize=16)

ax2.xaxis.set_major_formatter(FuncFormatter(lambda val, pos: f"{val*10:.0f}"))
ax2.yaxis.set_major_formatter(FuncFormatter(lambda val, pos: f"{val*10:.0f}"))

ax2.tick_params(axis="both", labelsize=14)
ax2.set_ylim(0, 8)
ax2.legend(loc="upper right", fontsize=11)

cbar = fig.colorbar(pcm, ax=ax2)
cbar.set_label(r"$Q_{k}$ (erg cm$^{-3}$ s$^{-1}$)", fontsize=14)
cbar.ax.tick_params(labelsize=12)

#plt.show()
plt.savefig(path1 / 'field_line_tracing.png')

# In[ ]:
