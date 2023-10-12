# This program generates an animation from AMRVAC data (dat format), and plot the magnetic field
#
# =============================================================================
#
#                                   Manual
#
# =============================================================================
#
#if you want to plot animation with denisty only in cylindrical coordinate:
#python gmunu_animation_rho_bfield.py --c cylindrical
#
#if you want to plot the density and the magnetic field: Note that the magnetic field plot does not support spherical coordinate
#python gmunu_animation_rho_bfield.py --c cylindrical --bfield

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
import yt
#import matplotlib
#matplotlib.use('Agg')
from matplotlib.animation import FuncAnimation
from matplotlib import rc_context
#import pylab # to show the plot
from optparse import OptionParser
import os
#import glob

#
# =============================================================================
#
#                                 Command Line
#
# =============================================================================
#
def parse_command_line():
	parser = OptionParser(description = __doc__)
	#extraction cooradinate
	parser.add_option("--c", "--coordinate-system", metavar = "coordinate system", type = 'string', help = 'Coordinate system used in the simulation.')

	parser.add_option("--data", default = './output*.dat', metavar = "filename", type = "string" , help = "Provide the input files. Example: --data './output000*dat' ")

	#plot related
	parser.add_option("--plane", metavar = "plane", type = int, help = "choose the plane to plot. choose 1 / 2 / 3")
	parser.add_option("--bfield", default = False, action = "store_true", help = "plot the magnetic field, default is False.")

	parser.add_option("--unit", default = 'code', metavar = "unit system", help = "Unit system of the animation.")
	
	#output data 
	parser.add_option("--o", "--output-dir", default = './', metavar = "directory", help = "Directory for ouputting the plots, default is ./ ")
	options, filenames = parser.parse_args()
	
	#check whether the input coordinate system is valid
	valid_coordinate_systems = ["cartesian", "cylindrical", "spherical"]
	valid_plane = [1, 2, 3]
	if options.c not in valid_coordinate_systems:
		raise ValueError("Invalid coordinate system, please input cartesian / cylindrical / spherical")

	if (options.plane not in  valid_plane) or (options.plane is None):
		raise ValueError("Invalid plane option, plase choose 1 / 2 / 3.")
	if options.bfield and (options.c == "spherical"):
		raise ValueError("Animation plotting does not support spherical coordinate.")
	return options, filenames



#
# =============================================================================
#
#                                  Utilities
#
# =============================================================================
#
def animate(i):
	ds = time_series[i]
	plot._switch_ds(ds)

#
# =============================================================================
#
#                                     Main
#
# =============================================================================
#
options, filenames = parse_command_line()
data_sets = options.data
if not os.path.exists(options.o):
	os.makedirs(options.o)
time_series = yt.load(data_sets, geometry_override=options.c, unit_system=options.unit)

plot = yt.SlicePlot(time_series[0], options.plane, 'rho', center=(50,0,0))
plot.set_cmap(field="rho", cmap='hot')
plot.set_zlim('rho', 1.4E-3, 1e-10)

#plot the magnetic field
if options.bfield:
	plot.annotate_streamlines('Bvec1', 'Bvec2', factor=16, density=1, plot_args={ "linewidth": 1.0, "color":'green'})

Radius = time_series[0].domain_width[0]
plot.set_center((Radius/2.0,Radius/2.0))
plot.set_width(Radius, Radius)
plot.annotate_timestamp(corner='upper_right', draw_inset_box=True)

fig = plot.plots['rho'].figure

animation = FuncAnimation(fig, animate, frames=len(time_series))

# Override matplotlib's defaults to get a nicer looking font
with rc_context({'mathtext.fontset': 'stix'}):
    animation.save('animation.mp4', fps=20)
