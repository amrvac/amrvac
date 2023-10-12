#This program extracts data from multiple points defined by users and FFT.
#
# =============================================================================
#
#                                   Manual
#
# =============================================================================
#
#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
from matplotlib import rc
import matplotlib.pyplot as plt
import yt
from optparse import OptionParser
import pandas as pd
import glob
import scipy.signal as signal
from scipy.fftpack import fft
from scipy.interpolate import interp1d
from multiprocessing import Pool
import os
import sys
import numpy as np

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

	parser.add_option("--data", default = './output*.dat', metavar = "filename", type = "string" , help = "Provide the input files. Example: ==data './output000*dat' ")

	#radius
	parser.add_option("--radius", metavar = "radius", type = 'float', help = 'Radius in LogFile.dat')	
	
	#fft related
	parser.add_option("--fft", default = False, action = "store_true", help = "Perform fft, default is False.")
	#parser.add_option("--window", default = True, action = "store_true", help = "multiply the time domain data with tukey window function before fft.")
	
	#output data
	parser.add_option("--o", "--output-dir", default = './', metavar = "directory", help = "Directory for ouputting the plots, default is ./ ")
	#output directory is stored as options.o

	#yt.load related
	parser.add_option("--unit", default = 'code', metavar = "unit system", help = "Unit system, default = code.")

	#multiprocessing related
	parser.add_option("--cores", default = 32, metavar = "number of cores", type = 'int', help = "Number of cores being used")
	
	#averaging the plot the extracted and FFTed data only
	parser.add_option("--s", "--skip-extraction-and-FFT", default = False, action = "store_true", help = "Skip extraction and FFT, obtain the plots directly")
	
	#plot related
	parser.add_option("--frequency-lower-limit", default = 0.0, metavar = "frequency (Hz)", type = 'float', help = 'frequency lower limit of the fft plot.')
	parser.add_option("--frequency-upper-limit", default = 8000.0, metavar = "frequency (Hz)", type = 'float', help = 'frequency lower limit of the fft plot.')
	
	#peak extraction related
	parser.add_option("--efl", "--extract-frequency-lower-limit", default = 200.0, metavar = "frequency (Hz)", type = 'float', help = 'frequency lower limit of peak extraction.')
	parser.add_option("--efu","--extract-frequency-upper-limit", default = 8000.0, metavar = "frequency (Hz)", type = 'float', help = 'frequency lower limit of peak extaction.')

	options, filenames = parser.parse_args()


	valid_coordinate_systems = ["cartesian", "cylindrical", "spherical"]
	if options.c not in valid_coordinate_systems:
		raise ValueError("Invalid coordinate system, please input cartesian / cylindrical / spherical")
	if options.radius is None and options.s is False:
		raise ValueError("Radius cannot be empty")
	return options, filenames
#
# =============================================================================
#
#                                  Utilities
#
# =============================================================================
def initialize_polar_grid(nRadial, nTheta, radius):
	R = np.array([j * radius / nRadial for j in range(1, nRadial+1)])
	Theta = np.radians(np.linspace(0, 179, nTheta))
	work = np.zeros(((nRadial * nTheta)+1, 4))
	return R, Theta, work

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

def get_amplitude(vector1, vector2, vector3, coordinate_system, x1_cut, x2_cut, x3_cut):
	if coordinate_system == "cartesian":
		w = vector1**2 + vector2**2 + vector3**2
	elif coordinate_system == "cylindrical":
		w = vector1**2 + vector2**2 + x1_cut**2 * vector3**2
	elif coordinate_system == "spherical":
		w = vector1**2 + x1_cut**2 * vector2**2 + x1_cut**2 * np.sin(x2_cut)**2 * vector3**2
	return w

#transform spherical grid into the coordinate system 
#since the grid is defined on x-z plane in spherical coordinate
def coordinate_transform_xz_plane(coordinate_system, r_cut, theta_cut, phi_cut):
	if coordinate_system == "cartesian":
		x1_cut = r_cut * np.cos(theta_cut)
		x2_cut = 0.
		x3_cut = r_cut * np.sin(theta_cut)
	else:
		x1_cut = r_cut
		x2_cut = theta_cut
		x3_cut = phi_cut
	return x1_cut, x2_cut, x3_cut

#generate file path with 3 digits number tag
def generate_path(number_tag, directory, file_name):
	filepath = directory + file_name + '_%03d.csv'%(number_tag)
	return filepath

def FFT_time(t_min, t_max):
	time_step = 1.0/2.0**12 * t_max
	time = np.arange(0.0, float(t_max), float(time_step))
	return time_step, time

def denoising_and_FFT(data_t, t_max):
	t_min = 0.0
	time_step, time = FFT_time(t_min, t_max)
	window = signal.tukey(len(time), alpha = 0.5)
	data_temp = data_t(time)
	data_temp -= np.mean(data_temp)
	data_temp *= window
	data_f = fft(data_temp) * 2.0/time.size
	data_psd = np.abs(data_f)**2
	return data_temp, data_f, data_psd

def FFT_main(df_TD, loc_p, coordinate_system, counter):
	print ('fft point {}...'.format(int(counter)))
	if coordinate_system == "spherical":
		df_TD = df_TD[['t','rho','psi','alp','vel1','vel2','vel3']]
	else:
		df_TD = df_TD[['t','rho','psi','alp','vel1','vel2','vel3', 'velr_sphe', 'veltheta_sphe', 'velphi_sphe']]

	df_TD['t'] = df_TD['t'] / 2.03001708E05 # code unit to sec
	t_min = 0.0
	t_max = float(df_TD.nlargest(1,'t')['t'])
	time_step, time = FFT_time(t_min, t_max)
	data_t0 = interp1d(df_TD['t'], df_TD['rho'], kind = 'cubic')
	data_t1 = interp1d(df_TD['t'], df_TD['vel1'], kind = 'cubic')
	data_t2 = interp1d(df_TD['t'], df_TD['vel2'], kind = 'cubic')
	data_t3 = interp1d(df_TD['t'], df_TD['vel3'], kind = 'cubic')
	data_temp0, data_f0, data_psd0 = denoising_and_FFT(data_t0, t_max)
	data_temp1, data_f1, data_psd1 = denoising_and_FFT(data_t1, t_max)
	data_temp2, data_f2, data_psd2 = denoising_and_FFT(data_t2, t_max)
	data_temp3, data_f3, data_psd3 = denoising_and_FFT(data_t3, t_max)
	if coordinate_system != "spherical":
		data_t4 = interp1d(df_TD['t'], df_TD['velr_sphe'], kind = 'cubic')
		data_t5 = interp1d(df_TD['t'], df_TD['veltheta_sphe'], kind = 'cubic')
		data_t6 = interp1d(df_TD['t'], df_TD['velphi_sphe'], kind = 'cubic')
		data_temp4, data_f4, data_psd4 = denoising_and_FFT(data_t4, t_max)
		data_temp5, data_f5, data_psd5 = denoising_and_FFT(data_t5, t_max)
		data_temp6, data_f6, data_psd6 = denoising_and_FFT(data_t6, t_max)
	freq = np.linspace(start=0.0, stop=1.0/(2.0 * time_step), num=int(time.size/2) )
	FD_part1 = pd.DataFrame({'loc': loc_p})
	if coordinate_system == "spherical":
		FD_part2 = pd.DataFrame({'f': freq, 'eigen0': np.multiply(np.sign(np.real(data_f0[:time.size//2])),np.abs(data_f0[:time.size//2])), 'fft0': np.abs(data_f0[:time.size//2]), 'psd0': data_psd0[:time.size//2], 'eigen1': np.multiply(np.sign(np.real(data_f1[:time.size//2])),np.abs(data_f1[:time.size//2])),'fft1': np.abs(data_f1[:time.size//2]), 'psd1': data_psd1[:time.size//2], 'eigen2': np.multiply(np.sign(np.real(data_f2[:time.size//2])),np.abs(data_f2[:time.size//2])), 'fft2': np.abs(data_f2[:time.size//2]), 'psd2': data_psd2[:time.size//2], 'eigen3': np.multiply(np.sign(np.real(data_f3[:time.size//2])),np.abs(data_f3[:time.size//2])), 'fft3': np.abs(data_f3[:time.size//2]), 'psd3': data_psd3[:time.size//2] })	
		df_FD = pd.concat([FD_part1, FD_part2], axis=1)
		df_FD = df_FD[['loc', 'f','eigen0','fft0','psd0','eigen1','fft1','psd1','eigen2','fft2','psd2','eigen3','fft3','psd3']]
	else:
		FD_part2 = pd.DataFrame({'f': freq, 'eigen0': np.multiply(np.sign(np.real(data_f0[:time.size//2])),np.abs(data_f0[:time.size//2])), 'fft0': np.abs(data_f0[:time.size//2]), 'psd0': data_psd0[:time.size//2], 'eigen1': np.multiply(np.sign(np.real(data_f1[:time.size//2])),np.abs(data_f1[:time.size//2])),'fft1': np.abs(data_f1[:time.size//2]), 'psd1': data_psd1[:time.size//2], 'eigen2': np.multiply(np.sign(np.real(data_f2[:time.size//2])),np.abs(data_f2[:time.size//2])), 'fft2': np.abs(data_f2[:time.size//2]), 'psd2': data_psd2[:time.size//2], 'eigen3': np.multiply(np.sign(np.real(data_f3[:time.size//2])),np.abs(data_f3[:time.size//2])), 'fft3': np.abs(data_f3[:time.size//2]), 'psd3': data_psd3[:time.size//2], 'eigen4': np.multiply(np.sign(np.real(data_f4[:time.size//2])),np.abs(data_f4[:time.size//2])), 'fft4': np.abs(data_f4[:time.size//2]), 'psd4': data_psd4[:time.size//2], 'eigen5': np.multiply(np.sign(np.real(data_f5[:time.size//2])),np.abs(data_f5[:time.size//2])), 'fft5': np.abs(data_f5[:time.size//2]), 'psd5': data_psd5[:time.size//2],'eigen6': np.multiply(np.sign(np.real(data_f6[:time.size//2])),np.abs(data_f6[:time.size//2])), 'fft6': np.abs(data_f6[:time.size//2]), 'psd6': data_psd6[:time.size//2] })
		df_FD = pd.concat([FD_part1, FD_part2], axis=1)
		df_FD = df_FD[['loc', 'f','eigen0','fft0','psd0','eigen1','fft1','psd1','eigen2','fft2','psd2','eigen3','fft3','psd3','eigen4','fft4','psd4','eigen5','fft5','psd5','eigen6','fft6','psd6']]
	return df_FD

def average_data():
	##averging FFT data of 360 points
	print ("loading FD data ...")
	df_FD_average = pd.read_csv('FD_data_polar/W_vel_FD_000.csv', delimiter = '\t')
	print ("loading TD data ...")
	df_TD_average = pd.read_csv('point_data_polar/point_data_000.csv', delimiter = '\t')
	print ("averaging ...")
	for x in range(1, 361, 1):
		file_path_TD = generate_path(x, 'point_data_polar/', 'point_data')
		file_path_FD = generate_path(x, 'FD_data_polar/', 'W_vel_FD')
		df_TD_tmp = pd.read_csv(file_path_TD, delimiter = '\t')
		df_FD_tmp = pd.read_csv(file_path_FD, delimiter = '\t')
		df_FD_tmp.fft2
		df_TD_average = df_TD_average.add(df_TD_tmp)
		df_FD_average = df_FD_average.add(df_FD_tmp)
	df_TD_average = df_TD_average.div(361)
	df_FD_average = df_FD_average.div(361)
	print ("saving ...")
	df_TD_average.to_csv('average_point_data.csv', index = False, sep = '\t')
	df_FD_average.to_csv('average_FD_data.csv', index = False, sep = '\t')
	return df_TD_average, df_FD_average
	
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
	plt.savefig("{}/averaged_rho_vel_TD_FD_combined.png".format(options.o), bbox_inches="tight", dpi= 300)
	return 
#
# =============================================================================
#
#                                     Main
#
# =============================================================================
#
def point(argu):
	r_cut = argu[0]
	theta_cut = argu[1]
	phi_cut = argu[2]
	counter = argu[3]
	x1_cut, x2_cut, x3_cut = coordinate_transform_xz_plane(options.c, r_cut, theta_cut, phi_cut)
	loc_p = np.array([x1_cut, x2_cut, x3_cut])
	#print ("initial coordinate {} {} {}".format(r_cut, theta_cut, phi_cut))
	#print ("transformed to {} {} {}".format(x1_cut, x2_cut, x3_cut))
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
			x1, x2, x3 = initialize_coordinate_variables(options.c)
			slice_data = ds.ortho_ray( x1, (x2_cut, x3_cut) )
			# locate where is the data
			data_loc = min(range(len(slice_data[x1])), key=lambda i: abs(float(slice_data[x1][i])-x1_cut))
			x1_local = float(slice_data[x1][data_loc])
			x2_local = float(slice_data[x2][data_loc])
			x3_local = float(slice_data[x3][data_loc])
			time = np.append(time,[ds.current_time])
			rho_p = np.append(rho_p, [slice_data['rho'][data_loc]])
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
			w = get_amplitude(wv1, wv2, wv3, options.c, r_cut, theta_cut, phi_cut)
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
				velphi_sphe   = np.append(velphi_sphe,   [-sin_phi*v1 + cos_phi*v2])
			elif options.c == "cylindrical":
				r_sphe = (x1_local**2 + x2_local**2)**0.5
				sin_theta = x1_local/r_sphe
				cos_theta = x2_local/r_sphe
				velr_sphe     = np.append(velr_sphe,     [sin_theta*v1 + cos_theta*v2])
				veltheta_sphe = np.append(veltheta_sphe, [cos_theta*v1 - sin_theta*v2])
				veltheta_sphe /= r_sphe
				velphi_sphe   = np.append(velphi_sphe,   [v3])
			
		# Sort the values by 't'
		srt = np.argsort(time)
	#save the time domain data
		TD_part1 = pd.DataFrame({'loc' : loc_p})
		if options.c == "spherical":
			TD_part2 = pd.DataFrame({'t' : time[srt], 'rho' : rho_p[srt], 'psi' : psi_p[srt], 'alp' : alp_p[srt], 'W_vel1' : W_vel1_p[srt], 'W_vel2' : W_vel2_p[srt], 'W_vel3' : W_vel3_p[srt], 'vel1' : vel1_p[srt], 'vel2' : vel2_p[srt], 'vel3' : vel3_p[srt]   }  )
			df_TD = pd.concat([TD_part1, TD_part2], axis=1)
			df_TD = df_TD[['loc','t','rho','psi','alp','W_vel1','W_vel2','W_vel3','vel1','vel2','vel3']]
		else:
			TD_part2 = pd.DataFrame({'t' : time[srt], 'rho' : rho_p[srt], 'psi' : psi_p[srt], 'alp' : alp_p[srt], 'W_vel1' : W_vel1_p[srt], 'W_vel2' : W_vel2_p[srt], 'W_vel3' : W_vel3_p[srt], 'vel1' : vel1_p[srt], 'vel2' : vel2_p[srt], 'vel3' : vel3_p[srt], 'velr_sphe' : velr_sphe[srt], 'veltheta_sphe' : veltheta_sphe[srt], 'velphi_sphe' : velphi_sphe[srt]   }  )
			df_TD = pd.concat([TD_part1, TD_part2], axis=1)
			df_TD = df_TD[['t','rho','psi','alp','W_vel1','W_vel2','W_vel3','vel1','vel2','vel3','velr_sphe','veltheta_sphe','velphi_sphe']] 
		df_TD.csv_path = generate_path(counter, 'point_data_polar/', 'point_data')
		df_TD.to_csv(df_TD.csv_path, index = False, sep = '\t')
	#skip the extraction 
	else:
		df_TD_csv_path = generate_path(counter, 'point_data_polar/', 'point_data')
		if options.c == "spherical":
			df_TD = pd.read_csv(df_TD_csv_path, delimiter = '\t', usecols = ['t', 'rho', 'psi', 'alp', 'vel1', 'vel2', 'vel3'])
		else:
			df_TD = pd.read_csv(df_TD_csv_path, delimiter = '\t', usecols = ['t','rho','psi','alp','W_vel1','W_vel2','W_vel3','vel1','vel2','vel3','velr_sphe','veltheta_sphe','velphi_sphe'])
		
	if options.fft:
		df_FD = FFT_main(df_TD, loc_p, options.c, counter)
		df_FD.csv_path = generate_path(counter, 'FD_data_polar/', 'W_vel_FD')
		df_FD.to_csv(df_FD.csv_path, index = False, sep = '\t')

## main
options, filenames = parse_command_line()
poolsize = options.cores
#prepare output directory
if not os.path.exists(options.o):
	os.makedirs(options.o)
data_sets = glob.glob(options.data)
#current_directory = options.o
folder1 = os.path.join(options.o, 'point_data_polar')
folder2 = os.path.join(options.o, 'FD_data_polar')
if not os.path.exists(folder1):
	os.makedirs(folder1)
if not os.path.exists(folder2):
	os.makedirs(folder2)

if not options.s:
	#initialize array
	R, Theta, work = initialize_polar_grid(20, 18, options.radius)
	k = 1
	for r in R:
		for theta in Theta:
			work[k][0] = r
			work[k][1] = theta
			work[k][3] = k
			k = k + 1

	if __name__ == '__main__':
		p=Pool(processes = poolsize, maxtasksperchild = 1)
		p.map(point, work)
		p.close()
		p.join()

df_TD_average, df_FD_average = average_data()
df_FD_plot = df_FD_average
for index, row in df_FD_plot.iterrows():
	if (row.f > options.frequency_upper_limit) or (row.f < options.frequency_lower_limit):
		df_FD_plot = df_FD_plot.drop(index)
plot_all(df_TD_average ,df_FD_plot, options.c, options.frequency_lower_limit, options.frequency_upper_limit)
#interpolation setting
extra_points = 2
interval = 0.001
kind = 'cubic'
verbose = True
#initializing interpolation
peak_array = []
#filtering the frequency range before extraction
print ("number of freuqency elements before filtering is {}".format(len(df_FD_average)))
for index, row in df_FD_average.iterrows():
	if (row.f > options.efu) or (row.f < options.efl):
		df_FD_average = df_FD_average.drop(index)
print ("number of freuqency elements after filtering is {}".format(len(df_FD_average)))
peak_array.append(extract_peak(df_FD_average, extra_points, interval, kind, verbose))
print ("extracted peak frequency is  {} Hz".format(peak_array[0]))

















