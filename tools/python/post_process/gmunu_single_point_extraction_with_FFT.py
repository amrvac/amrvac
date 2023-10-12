#This program extracts information from a single point dat format time domain data and FFT, and interpolate frequency domain data to extract the peak freuqency.
#
# =============================================================================
#
#                                   Manual
#
# =============================================================================
#
#if you simply want to extract all spherical point data in the output directory with a specified coordinate, without FFT, simply use
#python gmunu_single_point_extraction_with_FFT.py --c spherical --x1-cut 4.0 --x2-cut 0.22 --x3-cut 0.22
#
#please be reminded that 'output*dat' includes special variable '*', so '' is needed to specifiy the filename / path when * is included.
# 
#if you want to input dat data / point_data.csv from a specific directory with path "stored_data", add --data stored_data/point_data.csv or --data 'stored_data/output*dat'
#python gmunu_single_point_extraction_with_FFT.py --c spherical --x1-cut 4.0 --x2-cut 0.22 --x3-cut 0.22 --data 'stored_data/output*dat'
#python gmunu_single_point_extraction_with_FFT.py --c spherical --x1-cut 4.0 --x2-cut 0.22 --x3-cut 0.22 --data stored_data/point_data.csv
#
#if you want to extraction the point data together with FFT, add --fft
#python gmunu_single_point_extraction_with_FFT.py --c spherical --x1-cut 4.0 --x2-cut 0.22 --x3-cut 0.22 --fft
#
#if you want to skip point extraction and FFT directly using point_data.csv, add --skip-extraction
#python gmunu_single_point_extraction_with_FFT.py --skip-extraction --c spherical --data point_data.csv --fft
#
#if you want to output the data and plots into a specific directory named plots_data, add --o plots_data
#python gmunu_single_point_extraction_with_FFT.py --c spherical --x1-cut 4.0 --x2-cut 0.22 --x3-cut 0.22 --fft --o plots_data
#
#Please check which arguments has default values, eg: x1-cut, frequency_upper_limit, fft etc

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
import yt
import os
import numpy as np
import pandas as pd
import glob
from optparse import OptionParser
from scipy.fftpack import fft, fftshift
from scipy.interpolate import interp1d
import scipy.signal as signal
from matplotlib import rc
import matplotlib.pyplot as plt

#
# =============================================================================
#
#                                 Command Line
#
# =============================================================================
#
def parse_command_line():
	parser = OptionParser(description = __doc__)
	#if you want to skip point extraction and perforom FFT directly, please activate --skip-extraction and 
	#specify the data as point.csv instead of output*.dat
	parser.add_option("--skip-extraction", default = False, action = "store_true", help = "choose to skip point extraction")
	
	#extraction cooradinate
	parser.add_option("--c", "--coordinate-system", metavar = "coordinate system", type = 'string', help = 'Coordinate system used in the simulation.')
	# coordinate system is stored as options.c
	parser.add_option("--x1-cut", default = 8.0, metavar = "x1 cutoff coordinate", type = 'float', help = 'Provide a cut in x1 coordinate. spherical: radial, cartesian: x, cylindrical: radial')
	parser.add_option("--x2-cut", default = 0.0, metavar = "x2 cutoff coordinate", type = 'float', help = 'Provide a cut in x2 coordinate. spherical: theta, cartesian: y, cylindrical: phi. ')
	parser.add_option("--x3-cut", default = 0.0, metavar = "x3 cutoff coordinate", type = 'float', help = 'Provide a cut in x3 coordinate. spherical: phi, cartesian: z, cylindrical: z.')
	
	#input data, if you choose skip point extraction, please input point.csv instead pf output*.dat
	parser.add_option("--data", default = './output*.dat', metavar = "filename", type = "string" , help = "Provide the input files. Example: --data './output000*dat' ")
	#parser.add_option("--log", default = './output.log', metavar = "filename", type = "string" , help = "Provide the input files. Example: --data './output.log'")


	#yt.load related 
	#parser.add_option("--geometry", default = 'spherical', metavar = "coordinate system", help = "Geometry override argument of yt.load.")
	parser.add_option("--unit", default = 'code', metavar = "unit system", help = "Unit system, default = code.")	

	#fft related
	parser.add_option("--fft", default = False, action = "store_true", help = "Perform fft, default is False.")
	parser.add_option("--window", default = True, action = "store_true", help = "multiply the time domain data with tukey window function before fft.")
	
	#plot related
	parser.add_option("--frequency-lower-limit", default = 300.0, metavar = "frequency (Hz)", type = 'float', help = 'frequency lower limit of the fft plot.')
	parser.add_option("--frequency-upper-limit", default = 8000.0, metavar = "frequency (Hz)", type = 'float', help = 'frequency lower limit of the fft plot.')
	
	#peak extraction related
	parser.add_option("--efl", "--extract-frequency-lower-limit", default = 0.0, metavar = "frequency (Hz)", type = 'float', help = 'frequency lower limit of peak extraction.')
	parser.add_option("--efu","--extract-frequency-upper-limit", default = 8000.0, metavar = "frequency (Hz)", type = 'float', help = 'frequency lower limit of peak extaction.')
	#output data
	parser.add_option("--o", "--output-dir", default = './', metavar = "directory", help = "Directory for ouputting the plots, default is ./ ")
	#output directory is stored as options.o
	options, filenames = parser.parse_args()
	
	#check whether the inputs are valid
	valid_coordinate_systems = ["cartesian", "cylindrical", "spherical"]
	if options.c not in valid_coordinate_systems:
		raise ValueError("Invalid coordinate system, please input cartesian / cylindrical / spherical")
	if options.skip_extraction and options.data is None:
		raise ValueError("Point_data.csv should be provided if you skip the extraction.")
	return options, filenames


#
# =============================================================================
#
#                                  Utilities
#
# =============================================================================
#
def initialize_coordinate_variables(coordinate_system):
	if coordinate_system == "cartesian":
		x1 = "x"
		x2 = "y"
		x3 = "z"
	elif coordinate_system == "cylindrical":
		x1 = "r"
		x2 = "z"
		x3 = "theta"
	elif coordinate_system == "spherical":
		x1 = "r" 
		x2 = "theta"
		x3 = "phi"
	return x1, x2, x3

def get_amplitude(vector1, vector2, vector3, coordinate_system):
	if coordinate_system == "cartesian":
		w = vector1**2 + vector2**2 + vector3**2
	elif coordinate_system == "cylindrical":
		w = vector1**2 + vector2**2 + options.x1_cut**2 * vector3**2
	elif coordinate_system == "spherical":
		w = vector1**2 + options.x1_cut**2 * vector2**2 + options.x1_cut**2 * np.sin(options.x2_cut)**2 * vector3**2
	return w

def fft_utils(df_TD, coordinate_system, use_window):
	print ('fft ...')
	if coordinate_system == "spherical":
		df_TD = df_TD[['t','rho','psi','alp','vel1','vel2','vel3']]
	else:
		df_TD = df_TD[['t','rho','psi','alp','vel1','vel2','vel3', 'velr_sphe', 'veltheta_sphe', 'velphi_sphe']]
	df_TD['t'] = df_TD['t'] / 2.03001708E05 # code unit to sec 
	t_min = 0.0
	t_max = float(df_TD.nlargest(1,'t')['t'])
	#print (t_max)
	time_step = 1.0/2.0**12 * t_max
	#print (time_step)
	time = np.arange(0.0, float(t_max), float(time_step))

	window = signal.tukey(len(time), alpha = 0.5)
	data_t0 = interp1d(df_TD['t'], df_TD['rho'], kind = 'cubic')
	data_t1 = interp1d(df_TD['t'], df_TD['vel1'], kind = 'cubic')
	data_t2 = interp1d(df_TD['t'], df_TD['vel2'], kind = 'cubic')
	data_t3 = interp1d(df_TD['t'], df_TD['vel3'], kind = 'cubic')
	data_temp0 = data_t0(time)
	data_temp1 = data_t1(time)# * 1.0e5
	data_temp2 = data_t2(time)# * 1.0e5
	data_temp3 = data_t3(time)
	# if there is a signal that is oscillating at some mean that is not 0, it will dominate the frequency domain psd - this is not something we want. to kill it off, we can remove the average
	data_temp0 -= np.mean(data_temp0)
	data_temp1 -= np.mean(data_temp1)
	data_temp2 -= np.mean(data_temp2)
	data_temp3 -= np.mean(data_temp3)
	if use_window:
		data_temp0 *= window
		data_temp1 *= window
		data_temp2 *= window
		data_temp3 *= window
	data_f0 = fft(data_temp0) * 2.0/time.size
	data_f1 = fft(data_temp1) * 2.0/time.size
	data_f2 = fft(data_temp2) * 2.0/time.size
	data_f3 = fft(data_temp3) * 2.0/time.size
	data_psd0 = np.abs(data_f0)**2
	data_psd1 = np.abs(data_f1)**2
	data_psd2 = np.abs(data_f2)**2
	data_psd3 = np.abs(data_f3)**2
	if coordinate_system != "spherical":
		data_t4 = interp1d(df_TD['t'], df_TD['velr_sphe'], kind = 'cubic')
		data_t5 = interp1d(df_TD['t'], df_TD['veltheta_sphe'], kind = 'cubic')
		data_t6 = interp1d(df_TD['t'], df_TD['velphi_sphe'], kind = 'cubic')
		data_temp4 = data_t4(time)
		data_temp5 = data_t5(time)
		data_temp6 = data_t6(time)
		data_temp4 -= np.mean(data_temp4)
		data_temp5 -= np.mean(data_temp5)
		data_temp6 -= np.mean(data_temp6)
		if use_window:
			data_temp4 *= window
			data_temp5 *= window
			data_temp6 *= window
		data_f4 = fft(data_temp4) * 2.0/time.size
		data_f5 = fft(data_temp5) * 2.0/time.size
		data_f6 = fft(data_temp6) * 2.0/time.size
		data_psd4 = np.abs(data_f4)**2
		data_psd5 = np.abs(data_f5)**2
		data_psd6 = np.abs(data_f6)**2
	freq = np.linspace(start=0.0, stop=1.0/(2.0 * time_step), num=int(time.size/2) )
	df_FD = pd.DataFrame( {'f': freq, 'eigen0': np.multiply(np.sign(np.real(data_f0[:time.size//2])),np.abs(data_f0[:time.size//2])), 'fft0': np.abs(data_f0[:time.size//2]), 'psd0': data_psd0[:time.size//2], 'eigen1': np.multiply(np.sign(np.real(data_f1[:time.size//2])),np.abs(data_f1[:time.size//2])),'fft1': np.abs(data_f1[:time.size//2]), 'psd1': data_psd1[:time.size//2], 'eigen2': np.multiply(np.sign(np.real(data_f2[:time.size//2])),np.abs(data_f2[:time.size//2])), 'fft2': np.abs(data_f2[:time.size//2]), 'psd2': data_psd2[:time.size//2], 'eigen3': np.multiply(np.sign(np.real(data_f3[:time.size//2])),np.abs(data_f3[:time.size//2])), 'fft3': np.abs(data_f3[:time.size//2]), 'psd3': data_psd3[:time.size//2] })
	df_FD = df_FD[['f','eigen0','fft0','psd0','eigen1','fft1','psd1','eigen2','fft2','psd2','eigen3','fft3','psd3']]
	if coordinate_system != "spherical":
		df_FD = pd.DataFrame( {'f': freq, 'eigen0': np.multiply(np.sign(np.real(data_f0[:time.size//2])),np.abs(data_f0[:time.size//2])), 'fft0': np.abs(data_f0[:time.size//2]), 'psd0': data_psd0[:time.size//2], 'eigen1': np.multiply(np.sign(np.real(data_f1[:time.size//2])),np.abs(data_f1[:time.size//2])),'fft1': np.abs(data_f1[:time.size//2]), 'psd1': data_psd1[:time.size//2], 'eigen2': np.multiply(np.sign(np.real(data_f2[:time.size//2])),np.abs(data_f2[:time.size//2])), 'fft2': np.abs(data_f2[:time.size//2]), 'psd2': data_psd2[:time.size//2], 'eigen3': np.multiply(np.sign(np.real(data_f3[:time.size//2])),np.abs(data_f3[:time.size//2])), 'fft3': np.abs(data_f3[:time.size//2]), 'psd3': data_psd3[:time.size//2], 'eigen4': np.multiply(np.sign(np.real(data_f4[:time.size//2])),np.abs(data_f4[:time.size//2])), 'fft4': np.abs(data_f4[:time.size//2]), 'psd4': data_psd4[:time.size//2], 'eigen5': np.multiply(np.sign(np.real(data_f5[:time.size//2])),np.abs(data_f5[:time.size//2])), 'fft5': np.abs(data_f5[:time.size//2]), 'psd5': data_psd5[:time.size//2],'eigen6': np.multiply(np.sign(np.real(data_f6[:time.size//2])),np.abs(data_f6[:time.size//2])), 'fft6': np.abs(data_f6[:time.size//2]), 'psd6': data_psd6[:time.size//2] })
		df_FD = df_FD[['f','eigen0','fft0','psd0','eigen1','fft1','psd1','eigen2','fft2','psd2','eigen3','fft3','psd3','eigen4','fft4','psd4','eigen5','fft5','psd5','eigen6','fft6','psd6']]
	
	return df_FD

#locate the index of the maximum amplitude and some points nearby 
def extract_max_frequency_index(df_FD, extra_points):
	#first filter the frequency rangg by the extract frequency limits
	#print ("number of freuqency elements before filtering is {}".format(len(df_FD)))
	#for index, row in df_FD.iterrows():
	#	if (row.f > extract_frequency_upper_limit) or (row.f < extract_frequency_lower_limit):
	#		df_FD = df_FD.drop(index)
	#print ("number of freuqency elements after filtering is {}".format(len(df_FD)))
	#assume the peak is the maximum amplitude of theta
	i = df_FD['fft2'].idxmax()
	if (i - extra_points) < 0 :
		i_min = 0
	else:
		i_min = i - extra_points
	i_max = i + extra_points
	return i_min, i, i_max

#select point around maximum amplitude for interpolation
def select_point(df_FD, extra_points):
	#extra_points: number of extra points from maximum amplitude
	frequency = []
	fft = []
	i_min, i, i_max = extract_max_frequency_index(df_FD, extra_points)
	max_fft = df_FD['fft2'][i]
	peak_frequency = df_FD['f'][i]
	for x in range(i_min, i_max+1):
		frequency.append(df_FD['f'][x])
		fft.append(df_FD['fft2'][x])
	frequency = np.array(frequency)
	fft = np.array(fft)
	return frequency, fft, max_fft, peak_frequency

# interpolate the fft array and frequency array 
def interpolate_peak_frequency(df_FD, frequency, fft, kind, extra_points, interval):
	print ("working on {} interpolation ...".format(kind))
	f = interp1d(frequency, fft, kind=kind)
	i_min, i, i_max = extract_max_frequency_index(df_FD, extra_points)
	xnew = np.arange(df_FD['f'][i_min], df_FD['f'][i_max], interval)
	ynew = f(xnew)
	df_interpolated = pd.DataFrame({'fft2_new' : xnew, 'psd2_new' : ynew})
	df_interpolated = df_interpolated[['fft2_new', 'psd2_new' ]]
	print ("writing interpolated_fft2_psd2.csv to {} ...".format(options.o))
	df_interpolated.to_csv("{}/interpolated_fft2_psd2.csv".format(options.o))
	#df1 stores the interpolated FD data
	df1_FD = pd.DataFrame(xnew, columns = ['frequency'])
	df1_FD['fft2'] = ynew.tolist()
	j = df1_FD['fft2'].idxmax()
	max_fft_new = df1_FD['fft2'][j]
	peak_frequency_new = df1_FD['frequency'][j]
	return max_fft_new, peak_frequency_new

#extract peak frequency from the interpolated frequency and fft array
def extract_peak(df_FD, extra_points, interval, kind, verbose):
	frequency, fft, max_fft, peak_frequency = select_point(df_FD, extra_points)
	max_fft_new, peak_frequency_new = interpolate_peak_frequency(df_FD, frequency, fft, kind, extra_points, interval)
	percentage_diff = (peak_frequency_new - peak_frequency) / peak_frequency * 100
	if verbose:
		print ("\t\tBefore interpolation\t\t\tAfter interpolation")
		print ("maximum fft:\t{}\t{}".format(max_fft, max_fft_new))
		print ("peak frequency:\t{} Hz\t{} Hz".format(peak_frequency, peak_frequency_new))
		print ("% difference of peak frquency = {} %".format(percentage_diff))
	return peak_frequency_new

def plot_all(df_TD, df_FD, coordinate_system, frequency_lower_limit, frequency_upper_limit):
	#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':10})
	print ("plotting...")
	df_TD['t'] = df_TD['t'] / 2.03001708E02
	#df_TD['rho'] = df_TD['rho'] / 1.6199E-18
	if coordinate_system == "spherical":
		df_FD = df_FD[['f','fft0', 'psd0', 'fft1', 'psd1', 'fft2', 'psd2', 'fft3', 'psd3']]
	else:
		df_FD = df_FD[['f','fft0', 'psd0', 'fft1', 'psd1', 'fft2', 'psd2', 'fft3', 'psd3', 'fft4', 'psd4', 'fft5', 'psd5', 'fft6', 'psd6']] 
	plot_color = ["black", "r","g","b","fuchsia"]
	ls = ["-","-.", "--", ":","-."]
	lw = [1,1.5,1.5,1.5,1.5]
	fig = plt.figure(1)
	fig, ((ax0, ax1, ax2, ax3), (ax4, ax5, ax6, ax7)) = plt.subplots(2, 4, figsize=(22, 9))
	#fig, ((ax0, ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8, ax9)) = plt.subplots(2, 4, figsize=(26, 9))
	#ax1: time domain plot, ax2: frequency domain plot
	ax0.set_ylabel('$\\rho(t)$')
	ax1.set_ylabel('$v^{r}(t)$')
	ax2.set_ylabel('$v^{\\theta}(t)$')
	ax3.set_ylabel('$v^{\\phi}(t)$')
	#ax4.set_ylabel('$|v(t)|$')
	ax0.set_xlabel('$t$ $\\rm{(ms)}$')
	ax1.set_xlabel('$t$ $\\rm{(ms)}$')
	ax2.set_xlabel('$t$ $\\rm{(ms)}$')
	ax3.set_xlabel('$t$ $\\rm{(ms)}$')
	#ax4.set_xlabel('$t$ $\\rm{(ms)}$')
	rho_TD = df_TD['rho']
	if coordinate_system == "spherical":
		v_r = df_TD['vel1']
		v_theta = df_TD['vel2']
		v_phi = df_TD['vel3']
	else:
		v_r = df_TD['velr_sphe']
		v_theta = df_TD['veltheta_sphe']
		v_phi = df_TD['velphi_sphe']
	#v_abs = np.sqrt(v_r**2 + 4.0**2 * v_theta**2 + 4.0**2 * np.sin(0.28)**2 * v_phi**2)
	ax0.plot(df_TD['t'], rho_TD, linestyle = '-' ,linewidth = 1.0 )
	ax1.plot(df_TD['t'], v_r, linestyle = '-' ,linewidth = 1.0 )
	ax2.plot(df_TD['t'], v_theta, linestyle = '-' ,linewidth = 1.0 )
	ax3.plot(df_TD['t'], v_phi, linestyle = '-' ,linewidth = 1.0 )
	#ax4.plot(df_TD['t'], v_abs, linestyle = '-' ,linewidth = 1.0 )
	ax4.set_ylabel('$\\rm{FFT\ of\ }$ $\\rho(t)$')
	ax5.set_ylabel('$\\rm{FFT\ of\ }$ $v^{r}(t)$')
	ax6.set_ylabel('$\\rm{FFT\ of\ }$ $v^{\\theta}(t)$')
	ax7.set_ylabel('$\\rm{FFT\ of\ }$ $v^{\\phi}(t)$')
	#ax8.set_ylabel('$\\rm{FFT\ of\ }$ $|v(t)|$')
	ax4.set_xlabel('$f$ $\\rm{(Hz)}$')
	ax5.set_xlabel('$f$ $\\rm{(Hz)}$')
	ax6.set_xlabel('$f$ $\\rm{(Hz)}$')
	ax7.set_xlabel('$f$ $\\rm{(Hz)}$')
	#ax8.set_xlabel('$f$ $\\rm{(Hz)}$')
	ax4.plot(df_FD['f'], df_FD['fft0'] ,linestyle = '-' ,linewidth = 1.0, label='$\\rm{FFT}$$(\\rho$' )
	if coordinate_system == "spherical":
		ax5.plot(df_FD['f'], df_FD['fft1'] ,linestyle = '-' ,linewidth = 1.0, label='$\\rm{FFT}$$(v^{r})$' )
		ax6.plot(df_FD['f'], df_FD['fft2'] ,linestyle = '-' ,linewidth = 1.0, label='$\\rm{FFT}$$(v^{\\theta})$' )
		ax7.plot(df_FD['f'], df_FD['fft3'] ,linestyle = '-' ,linewidth = 1.0, label='$\\rm{FFT}$$(v^{\\phi})$' )
	else:	
		ax5.plot(df_FD['f'], df_FD['fft4'] ,linestyle = '-' ,linewidth = 1.0, label='$\\rm{FFT}$$(v^{r})$' )
		ax6.plot(df_FD['f'], df_FD['fft5'] ,linestyle = '-' ,linewidth = 1.0, label='$\\rm{FFT}$$(v^{\\theta})$' )
		ax7.plot(df_FD['f'], df_FD['fft6'] ,linestyle = '-' ,linewidth = 1.0, label='$\\rm{FFT}$$(v^{\\phi})$' )
	#ax8.plot(df_FD['f'], df_FD['fft4'] ,linestyle = '-' ,linewidth = 1.0, label='$\\rm{FFT}$$(|v|$' )
	ax4.set_xlim([frequency_lower_limit, frequency_upper_limit])
	ax5.set_xlim([frequency_lower_limit, frequency_upper_limit])
	ax6.set_xlim([frequency_lower_limit, frequency_upper_limit])
	ax7.set_xlim([frequency_lower_limit, frequency_upper_limit])
	ax4.legend(loc='best')
	ax5.legend(loc='best')
	ax6.legend(loc='best')
	ax7.legend(loc='best')
	ax4.grid(True)
	ax5.grid(True)
	ax6.grid(True)
	ax7.grid(True)
	fig.tight_layout()
	plt.savefig("{}/rho_vel_TD_FD_combined.png".format(options.o), bbox_inches="tight", dpi= 300)
	return 
#
# =============================================================================
#
#                                     Main
#
# =============================================================================
#
options, filenames = parse_command_line()
data_sets = glob.glob(options.data)
if not os.path.exists(options.o):
	os.makedirs(options.o)

#extract point data or not 
if not options.skip_extraction:
	time = np.array([])
	rho_p = np.array([])
	psi_p = np.array([])
	alp_p = np.array([])
	W_vel1_p = np.array([])
	W_vel2_p = np.array([])
	W_vel3_p = np.array([])
	vel1_p = np.array([])
	vel2_p = np.array([])
	vel3_p = np.array([])
	velr_sphe = np.array([])
	veltheta_sphe = np.array([])
	velphi_sphe = np.array([])

	for i, data_set in enumerate(data_sets):
	  # Load the dataset
		ds = yt.load(data_set, geometry_override = options.c, unit_system = options.unit)
	  # get time
	  # get 1D data set
	  # locate where is the data
		x1, x2, x3 = initialize_coordinate_variables(options.c)
		slice_data = ds.ortho_ray( x1, (options.x2_cut, options.x3_cut) )
		data_loc = min(range(len(slice_data[x1])), key=lambda i: abs(float(slice_data[x1][i]) - options.x1_cut))
		x1_local = float(slice_data[x1][data_loc])
		x2_local = float(slice_data[x2][data_loc])
		x3_local = float(slice_data[x3][data_loc])
	  #
		time = np.append(time,[ds.current_time])
		rho_p = np.append(rho_p, [slice_data['rho'][data_loc]]) #extract density
		psi_p = np.append(psi_p, [slice_data['psi'][data_loc]]) 
		alp_p = np.append(alp_p, [slice_data['alp'][data_loc]]) 
		W_vel1_p = np.append(W_vel1_p, [slice_data['W_vel1'][data_loc]])
		W_vel2_p = np.append(W_vel2_p, [slice_data['W_vel2'][data_loc]])
		W_vel3_p = np.append(W_vel3_p, [slice_data['W_vel3'][data_loc]])
	  # now workout the veloc1 to 3
		wv1 = float(slice_data['W_vel1'][data_loc])
		wv2 = float(slice_data['W_vel2'][data_loc])
		wv3 = float(slice_data['W_vel3'][data_loc])
		psi4 = float(slice_data['psi'][data_loc]**4)
	  #w = wv1**2 + options.x1_cut**2 * wv2**2 + options.x1_cut**2 * np.sin(options.x2_cut)**2 * wv3**2
		w = get_amplitude(wv1, wv2, wv3, options.c)
		w = psi4 * w
		w = np.sqrt( 1.0 + w )
		v1 = wv1 / w
		v2 = wv2 / w
		v3 = wv3 / w
		vel1_p = np.append(vel1_p, [slice_data['W_vel1'][data_loc]/w])
		vel2_p = np.append(vel2_p, [slice_data['W_vel2'][data_loc]/w])
		vel3_p = np.append(vel3_p, [slice_data['W_vel3'][data_loc]/w])

	#transform vel_p into spherical coordinate
	if options.c == "cartesian":
		r_sphe = (x1_local**2 + x2_local**2 + x3_local**2)**0.5
		r_xy = (x1_local**2 + x2_local**2)**0.5
		sin_theta =  r_xy/r_sphe
		cos_theta =  x3_local/r_sphe
		sin_phi =  x2_local/r_xy
		cos_phi =  x1_local/r_xy
		velr_sphe     = np.append(velr_sphe,     [sin_theta*cos_phi*v1 + sin_theta*sin_phi*v2 + cos_theta*v3])
		veltheta_sphe = np.append(veltheta_sphe, [cos_theta*cos_phi*v1 + cos_theta*sin_phi*v2 - sin_theta*v3])
		veltheta_sphe /= r_sphe
		velphi_sphe   = np.append(velphi_sphe,   [-sin_phi*v1 + cos_phi*v2]) 
		velphi_sphe   /= r_xy
	elif options.c == "cylindrical":
		r_sphe = (x1_local**2 + x2_local**2)**0.5
		sin_theta = x1_local/r_sphe
		cos_theta = x2_local/r_sphe
		velr_sphe = np.append(velr_sphe, [sin_theta*v1 + cos_theta*v2])
		veltheta_sphe = np.append(veltheta_sphe, [cos_theta*v1 - sin_theta*v2])
		veltheta_sphe /= r_sphe
		velphi_sphe = np.append(velphi_sphe, [v3])

	# Sort the values by 't'
	srt = np.argsort(time)

	#save the dataframe according to the coordinate
	if options.c == "spherical":
		df_TD = pd.DataFrame({'t' : time[srt], 'rho' : rho_p[srt], 'psi' : psi_p[srt], 'alp' : alp_p[srt], 'W_vel1' : W_vel1_p[srt], 'W_vel2' : W_vel2_p[srt], 'W_vel3' : W_vel3_p[srt], 'vel1' : vel1_p[srt], 'vel2' : vel2_p[srt], 'vel3' : vel3_p[srt]   }  )
		df_TD = df_TD[['t', 'rho', 'psi', 'alp', 'W_vel1', 'W_vel2', 'W_vel3', 'vel1', 'vel2', 'vel3']]
	else:
		df_TD = pd.DataFrame({'t' : time[srt], 'rho' : rho_p[srt], 'psi' : psi_p[srt], 'alp' : alp_p[srt], 'W_vel1' : W_vel1_p[srt], 'W_vel2' : W_vel2_p[srt], 'W_vel3' : W_vel3_p[srt], 'vel1' : vel1_p[srt], 'vel2' : vel2_p[srt], 'vel3' : vel3_p[srt], 'velr_sphe' : velr_sphe[srt], 'veltheta_sphe' : veltheta_sphe[srt], 'velphi_sphe' : velphi_sphe[srt]   }  )
		df_TD = df_TD[['t','rho','psi','alp','W_vel1','W_vel2','W_vel3','vel1','vel2','vel3','velr_sphe','veltheta_sphe','velphi_sphe']]

	print ("writing point_data.csv to {} ...".format(options.o))
	df_TD.to_csv("{}/point_data.csv".format(options.o) , index = False, sep = '\t')
else:
	if options.c == "spherical":
		df_TD = pd.read_csv(options.data, delimiter = '\t', usecols = ['t', 'rho', 'psi', 'alp', 'vel1', 'vel2', 'vel3'])
	else:
		df_TD = pd.read_csv(options.data, delimiter = '\t', usecols = ['t','rho','psi','alp','W_vel1','W_vel2','W_vel3','vel1','vel2','vel3','velr_sphe','veltheta_sphe','velphi_sphe'])


if options.fft:
	df_FD = fft_utils(df_TD, options.c, options.window)
	print ("writing W_vel_FD.csv to {} ...".format(options.o))
	df_FD.to_csv("{}/W_vel_FD.csv".format(options.o), index = False, sep = '\t')
	df_FD_plot = df_FD
	for index, row in df_FD_plot.iterrows():
		if (row.f > options.frequency_upper_limit) or (row.f < options.frequency_lower_limit):
			df_FD_plot = df_FD_plot.drop(index)	
	plot_all(df_TD ,df_FD_plot, options.c, options.frequency_lower_limit, options.frequency_upper_limit)
	#interpolation setting
	extra_points = 2
	interval = 0.001
	kind = 'cubic'
	verbose = True
	#initializing interpolation
	peak_array = []
	#filtering the frequency range before extraction 
	print ("number of freuqency elements before filtering is {}".format(len(df_FD)))
	for index, row in df_FD.iterrows():
		if (row.f > options.efu) or (row.f < options.efl):
			df_FD = df_FD.drop(index)
	print ("number of freuqency elements after filtering is {}".format(len(df_FD)))
	peak_array.append(extract_peak(df_FD, extra_points, interval, kind, verbose))
	print ("extracted peak frequency is  {} Hz".format(peak_array[0]))
