#This program plots the magnetic energy.
#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from matplotlib import rc_context
from matplotlib import rc
from optparse import OptionParser

#
# =============================================================================
#
#                                 Command Line
#
# =============================================================================
#

def parse_command_line():
	parser = OptionParser(description = __doc__)
	parser.add_option("--log", default = 'output.log', metavar = "filename", type = "string" , help = "Provide the input files. Example: --log output.log ")
	parser.add_option("--output", default = './B_energy.png', metavar = "filename", help = "Directory for ouput and the file name.")
	options, filenames = parser.parse_args()
	return options, filenames

#
# =============================================================================
#
#                                  Utilities
#
# =============================================================================
#

# =============================================================================
#
#                                     Main
#
# =============================================================================
#
options, filenames = parse_command_line()
df = pd.read_csv(options.log, sep = r'\s{1,}')
# convert code unit to milisecond
df['time_ms'] = df['global_time'] / 2.03001708e05 * 1000 
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':15})
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

plot_color = ["black", "r","g","b","fuchsia"]
ls = ["-","-.", "--", ":","-."]
lw = [1,1.5,1.5,1.5,1.5]

fig = plt.figure(1)
fig, (ax1) = plt.subplots(1, 1)
ax1.set_xlabel('$t$ $\\rm{(ms)}$')

ax1.plot(df['time_ms'] , df['EB_pol'], linestyle = '-' ,linewidth = 1.0, label='$\\mathcal{E}_{\\rm{pol}}(t) / \\mathcal{E}_{\\rm{mag}}(0)$' )
ax1.plot(df['time_ms'] , df['EB_tor'], linestyle = '-' ,linewidth = 1.0, label='$\\mathcal{E}_{\\rm{tor}}(t) / \\mathcal{E}_{\\rm{mag}}(0)$' )
ax1.legend(loc='best')
ax1.grid(True)
fig.savefig(options.output, bbox_inches="tight")

