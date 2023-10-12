#I hate yt
from matplotlib import rc
import matplotlib.pyplot as plt
import yt
import pandas as pd
import glob
import scipy.signal as signal
from scipy.fftpack import fft
from scipy.interpolate import interp1d
from multiprocessing import Pool
import os
import sys
import numpy as np
import configparser
from itertools import repeat
import re

def check_config_value(config):
	valid_coordinate_systems = ["cartesian", "cylindrical", "spherical"]
	workflow_extraction = config.getboolean("workflow", "extraction")
	valid_extraction_method = ["single", "multiple"]

	if config["extraction setting"]["coordinate_system"] not in valid_coordinate_systems:
		raise ValueError("Invalid coordinate system, please input cartesian / cylindrical / spherical")

	if workflow_extraction and config["IO setting"]["datafilepath"] == "":
		raise ValueError("outputxxx.dat should be provided if you want to extract.")
	elif not workflow_extraction and config["IO setting"]["datafilepath"]  == "":
		raise ValueError("point_data.csv should be provided if you skip the extraction.")

	#if config["extraction setting"]["coordinate_system"] == "cylindrical" and config["extraction setting"]["extraction_method"] == "multiple":
	#	raise NotImplementedError("multiple points extractions in cylindrical coordinate is not functional")#best line learn from yt

	if config["extraction setting"]["extraction_method"] == "multiple" not in valid_extraction_method:

		raise ValueError("please input single / multiple in extraction method ")

def create_directory(config, directory):
	output_data_path = config["IO setting"]["output_data_dir"]
	folder = os.path.join(output_data_path, directory)
	if not os.path.exists(folder):
		os.makedirs(folder)

def initialize_polar_grid(config):
	initial_profile = config["multiple points extraction setting"]['initial_profile_path']
	index = 0
	# check each line in xns logfile, move index up by 1, read only mode
	for line in open(initial_profile,"r") :
		index += 1
		# if string is found in line
		if 'Equatorial Radius' in line :
		# extract equatorial radius float, in code units (the first float found)
			radius =  float(re.findall(r'[\d]*[.][\d]+', line)[0])
	print ("radius = {}".format(radius))
	nRadial = config.getint("multiple points extraction setting", "nradial")
	nTheta = config.getint("multiple points extraction setting", "ntheta")
	R = np.array([j * radius / nRadial for j in range(1, nRadial+1)])
	Theta = np.radians(np.linspace(0, 179, nTheta))
	work = np.zeros(((nRadial * nTheta)+1, 4))
	k = 1
	for r in R:
		for theta in Theta:
			work[k][0] = r
			work[k][1] = theta
			work[k][3] = k
			k = k + 1
	return work

def generate_csv_path(path, counter, name):
	try:
		csv_path = path + name + "_%03d"%counter + ".csv"
	except TypeError:
		csv_path = path + name + ".csv"

	return csv_path

def filter_FD_data(df_FD,lower, upper):
	for index, row in df_FD.iterrows():
		if (row.frequency > upper) or (row.frequency < lower):
			df_FD = df_FD.drop(index)
	return df_FD
#
# =============================================================================
#
#                          Extraction Utilities
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

def extraction_coordinates(argu, config):
	coordinate_system = config["extraction setting"]["coordinate_system"]
	extraction_method = config["extraction setting"]["extraction_method"]
	if extraction_method == "single":
		x1_cut = config.getfloat("single point extraction setting", "x1_cut")
		x2_cut = config.getfloat("single point extraction setting", "x2_cut")
		x3_cut = config.getfloat("single point extraction setting", "x3_cut") 
	elif extraction_method == "multiple":	
		r_cut = argu[0]
		theta_cut = argu[1]
		phi_cut = argu[2]
		
	#coordinate transform the polar grid coordinate to xy plane
		if coordinate_system == "cartesian":
			x1_cut = r_cut * np.cos(theta_cut)
			x2_cut = 0.
			x3_cut = r_cut * np.sin(theta_cut)
		elif coordinate_system == "cylindrical":
			x1_cut = r_cut 
			x2_cut = 0.
			x3_cut = theta_cut
		elif coordinate_system == "spherical":
			x1_cut = r_cut
			x2_cut = theta_cut
			x3_cut = phi_cut
		print ("coordinate transform:{} {} {} --> {} {} {}".format(r_cut, theta_cut, phi_cut, x1_cut, x2_cut, x3_cut))
	return x1_cut, x2_cut, x3_cut

#transform spherical grid into the coordinate system 
#since the grid is defined on x-y plane in spherical coordinate
def get_amplitude(vector1, vector2, vector3, coordinate_system, x1_cut, x2_cut, x3_cut):
	if coordinate_system == "cartesian":
		w = vector1**2 + vector2**2 + vector3**2
	elif coordinate_system == "cylindrical":
		w = vector1**2 + vector2**2 + x1_cut**2 * vector3**2
	elif coordinate_system == "spherical":
		w = vector1**2 + x1_cut**2 * vector2**2 + x1_cut**2 * np.sin(x2_cut)**2 * vector3**2
	return w



def extract_vtu

def point(argu, config):
	time = rho_p = psi_p = alp_p = W_vel1_p = W_vel2_p = W_vel3_p = vel1_p = vel2_p = vel3_p = velr_sphe = veltheta_sphe = velphi_sphe = omega_r = omega_theta= omega_phi = np.array([])
	extraction_method = config["extraction setting"]["extraction_method"]
	data_sets = glob.glob(config["IO setting"]["datafilepath"])
	coordinate_system = config["extraction setting"]["coordinate_system"]

	#generate the extraction coordinates
	x1_cut, x2_cut, x3_cut = extraction_coordinates(argu, config)
	loc_p = np.array([x1_cut, x2_cut, x3_cut])	
	for i, data_set in enumerate(data_sets):

		#load the dataset
		print ("loading {} ...".format(data_set))
		if data_set.endwith('.dat'):
			input("Please press the Enter key to proceed")
			ds = yt.load(data_set, geometry_override = coordinate_system, unit_system = config["extraction setting"]["unit_system"])
			x1, x2, x3 = initialize_coordinate_variables(coordinate_system)
			swap_flag = 0
			if coordinate_system == "cylindrical":
				slice_data = ds.ortho_ray( x1, (x3_cut, x2_cut) )
				swap_flag = 1
			if len(slice_data[x1]) == 0:
				slice_data = ds.ortho_ray( x1, (x2_cut, x3_cut) )
				print ("axis order in ortho_ray swapped .........")
			else:
				slice_data = ds.ortho_ray( x1, (x2_cut, x3_cut) )
				print ("x2_cut = {}, x3_cut = {} ..........., length of slice_data = {} ".format(x2_cut, x3_cut, len(slice_data[x1])))
		elif data_set.endwith('.vtu'):
			ds = read.load_vtkfile(i, file = data_set, type='vtu')
			rho = ds.rho
			psi = ds.psi
			alp = ds.alp
			W_vel1 = ds.W_vel1
			W_vel2 = ds.W_vel2
			W_vel3 = ds.W_vel3
			coord = ds.getCenterPoints()
			diff_min = 1000
			for j in range(len(coord)):     
				x1_j = coord[j,0]
				y2_j = coord[j,1]
				diff = np.sqrt((x1_cut - x_j)**2 + (y2_cut - y_j)**2)
			if diff < diff_min:
				diff_min = diff
				data_loc = j     

		#find the closest point to our cut in the data
		data_loc = min(range(len(slice_data[x1])), key=lambda i: abs(float(slice_data[x1][i]) - x1_cut))	
		x1_local = float(slice_data[x1][data_loc])
		x2_local = float(slice_data[x2][data_loc])
		x3_local = float(slice_data[x3][data_loc])

		#append every data to the numpy array
		time = np.append(time,[ds.current_time])
		rho_p = np.append(rho_p, [slice_data['rho'][data_loc]])
		psi_p = np.append(psi_p, [slice_data['psi'][data_loc]])
		alp_p = np.append(alp_p, [slice_data['alp'][data_loc]])
		W_vel1_p = np.append(W_vel1_p, [slice_data['W_vel1'][data_loc]])
		W_vel2_p = np.append(W_vel2_p, [slice_data['W_vel2'][data_loc]])
		W_vel3_p = np.append(W_vel3_p, [slice_data['W_vel3'][data_loc]])
		wv1 = float(slice_data['W_vel1'][data_loc])
		wv2 = float(slice_data['W_vel2'][data_loc])
		wv3 = float(slice_data['W_vel3'][data_loc])
		psi4 = float(slice_data['psi'][data_loc]**4)
		w = get_amplitude(wv1, wv2, wv3, coordinate_system, x1_cut, x2_cut, x3_cut)
		w = psi4 * w
		w = np.sqrt( 1.0 + w )
		v1 = wv1 / w
		v2 = wv2 / w
		v3 = wv3 / w
		vel1_p = np.append(vel1_p, [slice_data['W_vel1'][data_loc]/w])
		vel2_p = np.append(vel2_p, [slice_data['W_vel2'][data_loc]/w])
		vel3_p = np.append(vel3_p, [slice_data['W_vel3'][data_loc]/w])
		#calculate omega according to the coordinate system
		if coordinate_system == "cartesian":
			r_sphe = (x1_local**2 + x2_local**2 + x3_local**2)**0.5
			r_xy = (x1_local**2 + x2_local**2)**0.5
			sin_theta =  r_xy/r_sphe
			cos_theta =  x3_local/r_sphe
			sin_phi =  x2_local/r_xy
			cos_phi =  x1_local/r_xy
			omega_r = np.append(omega_r, [sin_theta*cos_phi*v1 + sin_theta*sin_phi*v2 + cos_theta*v3])
			omega_theta = np.append(omega_theta, [cos_theta*cos_phi*v1 + cos_theta*sin_phi*v2 - sin_theta*v3])
			omega_phi = np.append(omega_phi, [-sin_phi*v1 + cos_phi*v2])
		elif coordinate_system == "cylindrical":
			r_sphe = (x1_local**2 + x2_local**2)**0.5
			sin_theta = x1_local/r_sphe
			cos_theta = x2_local/r_sphe
			omega_r = np.append(omega_r, [cos_theta*v1 + sin_theta*v2])
			#ofcoz don't forget the swap issue due to ortho_ray
			if swap_flag == 0:
				omega_theta = np.append(omega_phi,[v3])
				omega_phi = np.append(omega_theta, [sin_theta*v1 - cos_theta*v2])
			elif swap_flag == 1:
				omega_theta = np.append(omega_theta, [sin_theta*v1 - cos_theta*v2])
				omega_phi = np.append(omega_phi,[v3])
	# transform vel into spherical
	if coordinate_system == "cartesian":
		velr_sphe = (vel1_p ** 2 + vel2_p **2 + vel3_p **2) ** 1/2
		veltheta_sphe = np.arctan(vel3_p / vel1_p)
		veltheta_sphe =  np.nan_to_num(veltheta_sphe, copy=True, nan=0.0, posinf=None, neginf=None)
		velphi_sphe = np.arctan( (( vel1_p ** 2 + vel2_p ** 2) ** 1/2) / vel3_p)
		velphi_sphe = np.nan_to_num(velphi_sphe, copy=True, nan=0.0, posinf=None, neginf=None)
	elif coordinate_system == "cylindrical":
		velr_sphe = (vel1_p ** 2 + vel2_p ** 2 ) ** 1/2
		#ofcoz don't forget the swap issue due to ortho_ray
		if swap_flag == 0:
			veltheta_sphe = vel3_p
			velphi_sphe = np.arctan(vel1_p / vel2_p)
		elif swap_flag == 1:
			veltheta_sphe = np.arctan(vel1_p / vel2_p)
			veltheta_sphe =  np.nan_to_num(veltheta_sphe, copy=True, nan=0.0, posinf=None, neginf=None)
			velphi_sphe = vel3_p

	# Sort the values by 't' before saving
	srt = np.argsort(time)
	#save the dataframe according to the coordinate
	TD_part1 = pd.DataFrame({'loc' : loc_p})
	if coordinate_system == "spherical":
		TD_part2 = pd.DataFrame({'t' : time[srt], 
					'rho' : rho_p[srt],
					'psi' : psi_p[srt],
					'alp' : alp_p[srt],
					'W_vel1' : W_vel1_p[srt],
					'W_vel2' : W_vel2_p[srt],
					'W_vel3' : W_vel3_p[srt],
					'vel1' : vel1_p[srt],
					'vel2' : vel2_p[srt],
					'vel3' : vel3_p[srt]   }  )
		TD_part2 = TD_part2[['t', 'rho', 'psi', 'alp', 'W_vel1', 'W_vel2', 'W_vel3', 'vel1', 'vel2', 'vel3']]
		df_TD = pd.concat([TD_part1, TD_part2], axis=1)
	else:
		TD_part2 = pd.DataFrame({'t' : time[srt],
					'rho' : rho_p[srt],
					'psi' : psi_p[srt],
					'alp' : alp_p[srt],
					'W_vel1' : W_vel1_p[srt],
					'W_vel2' : W_vel2_p[srt],
					'W_vel3' : W_vel3_p[srt],
					'vel1' : vel1_p[srt],
					'vel2' : vel2_p[srt],
					'vel3' : vel3_p[srt],
					'velr_sphe' : velr_sphe[srt],
					'veltheta_sphe' : veltheta_sphe[srt],
					'velphi_sphe' : velphi_sphe[srt],
					'omega_r': omega_r[srt],
					'omega_theta': omega_theta[srt],
					'omega_phi': omega_phi[srt]}  )
		TD_part2 = TD_part2[['t','rho','psi','alp','W_vel1','W_vel2','W_vel3','vel1','vel2','vel3','velr_sphe','veltheta_sphe','velphi_sphe', 'omega_r', 'omega_theta', 'omega_phi']]
		df_TD = pd.concat([TD_part1, TD_part2], axis=1)
		df_TD.fillna(0,inplace = True)

	output_data_path = config["IO setting"]["output_data_dir"]
	output_data_path = output_data_path + "point_data_polar/" if extraction_method == "multiple" else output_data_path

	counter = None if argu is None else argu[3]
	TD_csv_path = generate_csv_path(output_data_path, counter, "point_data")
	df_TD.to_csv(TD_csv_path, index = False , sep = '\t')
	return 

def extraction(config):
	poolsize = config.getint("multiple points extraction setting", "poolsize")
	extraction_setting = config["extraction setting"]
	if extraction_setting["extraction_method"] == "single":
		point(None, config)
	elif extraction_setting["extraction_method"] == "multiple":
		create_directory(config, "point_data_polar")
		work = initialize_polar_grid(config)

		if __name__ == '__main__':
			p=Pool(processes = poolsize, maxtasksperchild = 1)
			p.starmap(point, zip(work, repeat(config)))
			p.close()
			p.join()

#
# =============================================================================
#
#                               FFT Utilities
#
# =============================================================================
#
def average_data(work):
	max_num = len(work)
	for x in range(1, max_num, 1):
		TD_csv_path = generate_csv_path('point_data_polar/', x, 'point_data')
		FD_csv_path = generate_csv_path('FD_data_polar/', x, 'W_vel_FD')
		try:
			average_TD
			average_FD
		except NameError:
			average_TD = pd.read_csv(TD_csv_path, delimiter = '\t')
			average_FD = pd.read_csv(FD_csv_path, delimiter = '\t')
		else:
			TD_tmp = pd.read_csv(TD_csv_path, delimiter = '\t')
			FD_tmp = pd.read_csv(FD_csv_path, delimiter = '\t')
			average_TD = average_TD.add(TD_tmp)
			average_FD = average_FD.add(FD_tmp)
	average_TD = average_TD.div(max_num)
	average_FD = average_FD.div(max_num)
	average_TD.to_csv('average_point_data.csv', index = False, sep = '\t')
	average_FD.to_csv('average_FD_data.csv', index = False, sep = '\t')
	return 

def FFT_main(argu, config):
	coordinate_system = config["extraction setting"]["coordinate_system"]
	output_data_path = config["IO setting"]["output_data_dir"]
	extraction_method = config["extraction setting"]["extraction_method"]
	#import time domain data
	point_data_path = output_data_path + "point_data_polar/" if extraction_method == "multiple" else output_data_path
	counter = None if argu is None else argu[3]
	name = "point_data" if argu is None else "point_data_polar/point_data"
	TD_data_path = generate_csv_path(output_data_path, counter, name)
	df_TD = pd.read_csv(TD_data_path, delimiter = '\t')

	#select specific columns for FFT
	if coordinate_system == "spherical":
		df_for_fft = df_TD[['t', 'rho', 'vel1', 'vel2', 'vel3']]
	else:
		df_for_fft = df_TD[['t', 'rho', 'vel1', 'vel2', 'vel3', 'velr_sphe', 'veltheta_sphe', 'velphi_sphe', 'omega_r', 'omega_theta', 'omega_phi']]
	
	# code unit to sec
	df_TD['t'] = df_TD['t'] / 2.03001708E05 
	
	#FFT paramters
	t_min = 0.0
	t_max = float(df_TD.nlargest(1,'t')['t'])
	time_step = 1.0/2.0**12 * t_max
	time = np.arange(0.0, float(t_max), float(time_step))
	
	#denoising window
	window = signal.tukey(len(time), alpha = 0.5)

	#resultant dataframe: df_FD
	df_FD = pd.DataFrame()	
	freq = np.linspace(start=0.0, stop=1.0/(2.0 * time_step), num=int(time.size/2) )	
	df_FD["frequency"] = freq.tolist()
	#iterate the coloumn for FFT
	for name, values in df_for_fft.iteritems():
		data = interp1d(df_TD['t'], values, kind = 'cubic')
		data_temp = data(time)
		data_temp -= np.mean(data_temp)
		data_temp *= window
		data_FD = fft(data_temp) * 2.0/time.size
		data_fft_output = np.abs(data_FD[:time.size//2])
		data_psd = np.abs(data_FD)**2
		data_psd_output = data_psd[:time.size//2]
		data_fft_output = np.abs(data_FD[:time.size//2])
		data_eigen_output = np.multiply(np.sign(np.real(data_FD[:time.size//2])),np.abs(data_FD[:time.size//2]))
		df_FD["{}_fft".format(name)] = data_fft_output.tolist()
		df_FD["{}_psd".format(name)] = data_psd_output.tolist()
		df_FD["{}_eigen".format(name)] = data_eigen_output.tolist()


	#save
	output_data_path = output_data_path + "FD_data_polar/" if (extraction_method == "multiple") else output_data_path
	counter = None if argu is None else argu[3]

	FD_csv_path = generate_csv_path(output_data_path, counter, "W_vel_FD")	
	df_FD.to_csv(FD_csv_path, index = False , sep = '\t')
	
	return 

def FFT_execute(config):
	poolsize = config.getint("multiple points extraction setting", "poolsize")
	extraction_setting = config["extraction setting"]
	if extraction_setting["extraction_method"] == "single":
		FFT_main(None, config)
	elif extraction_setting["extraction_method"] == "multiple":
		create_directory(config, "FD_data_polar")
		work = initialize_polar_grid(config)

		if __name__ == '__main__':
			p=Pool(processes = poolsize, maxtasksperchild = 1)
			p.starmap(FFT_main, zip(work, repeat(config)))
			p.close()
			p.join()

		average_data(work)
	return 

#
# =============================================================================
#
#                               Peak Extraction Utilities
#
# =============================================================================
#

#locate the index of the maximum amplitude and some points nearby 
def extract_max_frequency_index(df_FD, extra_points, config):
	#first filter the frequency rangg by the extract frequency limits
	#print ("number of freuqency elements before filtering is {}".format(len(df_FD)))
	#for index, row in df_FD.iterrows():
	#	if (row.f > extract_frequency_upper_limit) or (row.f < extract_frequency_lower_limit):
	#		df_FD = df_FD.drop(index)
	#print ("number of freuqency elements after filtering is {}".format(len(df_FD)))
	#assume the peak is the maximum amplitude of theta
	coordinate_system = config["extraction setting"]["coordinate_system"]
	i = df_FD['vel2_fft'].idxmax() if coordinate_system == "spherical" else df_FD['omega_theta_fft'].idxmax()
	if (i - extra_points) < df_FD.index[0] or (i - extra_points) < 0:
		raise ValueError("Your extraction frequency lower limit hit the lower boundary of the FD data, no extra points can be selected for interpolation. Please adjust your extraction frequency lower limit in extraction.ini.")
	else:
		i_min = i - extra_points
	i_max = i + extra_points
	return i_min, i, i_max

#select point around maximum amplitude for interpolation
def select_point(df_FD, extra_points, config):
	#extra_points: number of extra points from maximum amplitude
	frequency = []
	fft = []
	i_min, i, i_max = extract_max_frequency_index(df_FD, extra_points, config)
	coordinate_system = config["extraction setting"]["coordinate_system"]
	max_fft = df_FD['vel2_fft'][i] if coordinate_system == "spherical" else df_FD['omega_theta_fft'][i]
	peak_frequency = df_FD['frequency'][i]
	for x in range(i_min, i_max+1):
		frequency.append(df_FD['frequency'][x])
		if coordinate_system == "spherical":
			fft.append(df_FD['vel2_fft'][x])
		else:
			fft.append(df_FD['omega_theta_fft'][x])
	frequency = np.array(frequency)
	fft = np.array(fft)
	return frequency, fft, max_fft, peak_frequency

# interpolate the fft array and frequency array 
def interpolate_peak_frequency(df_FD, frequency, fft, kind, extra_points, interval, config):
	print ("working on {} interpolation ...".format(kind))
	coordinate_system = config["extraction setting"]["coordinate_system"]
	f = interp1d(frequency, fft, kind=kind)
	i_min, i, i_max = extract_max_frequency_index(df_FD, extra_points, config)
	xnew = np.arange(df_FD['frequency'][i_min], df_FD['frequency'][i_max], interval)
	ynew = f(xnew)
	if coordinate_system == "spherical":
		df_interpolated = pd.DataFrame({'vel2_fft_new' : xnew, 'vel2_psd_new' : ynew})
		df_interpolated = df_interpolated[['vel2_fft_new', 'vel2_psd_new' ]]
		print ("writing interpolated_vel2_fft_psd.csv to {} ...".format(config['IO setting']["output_data_dir"]))
		df_interpolated.to_csv("{}/interpolated_vel2_fft_psd.csv".format(config['IO setting']["output_data_dir"]))
	else:
		df_interpolated = pd.DataFrame({'omega_theta_fft_new' : xnew, 'omega_theta_psd_new' : ynew})
		df_interpolated = df_interpolated[['omega_theta_fft_new', 'omega_theta_psd_new' ]]
		print ("writing interpolated_veltheta_sphe_fft_psd.csv to {} ...".format(config['IO setting']["output_data_dir"]))
		df_interpolated.to_csv("{}/interpolated_veltheta_sphe_fft_psd.csv".format(config['IO setting']["output_data_dir"]))
	#df1 stores the interpolated FD data
	df1_FD = pd.DataFrame(xnew, columns = ['frequency'])
	df1_FD['veltheta_fft'] = ynew.tolist()
	j = df1_FD['veltheta_fft'].idxmax()
	max_fft_new = df1_FD['veltheta_fft'][j]
	peak_frequency_new = df1_FD['frequency'][j]
	return max_fft_new, peak_frequency_new

#extract peak frequency from the interpolated frequency and fft array
def extract_peak(config):
	extraction_method = config["extraction setting"]["extraction_method"]
	frequency_lower_limit = config.getfloat("frequency peak extraction setting", "extraction_frequency_lower_limit")
	frequency_upper_limit = config.getfloat("frequency peak extraction setting", "extraction_frequency_upper_limit")
	output_data_path = config["IO setting"]["output_data_dir"]
		
	if extraction_method == "single":
		FD_data_path = generate_csv_path(output_data_path, None, "W_vel_FD")
	elif extraction_method == "multiple":
		FD_data_path = generate_csv_path(output_data_path, None, "average_FD_data")
	df_FD = pd.read_csv(FD_data_path, delimiter = '\t')
	df_FD = filter_FD_data(df_FD, frequency_lower_limit, frequency_upper_limit)
	extra_points = 2
	interval = 0.001
	kind = 'cubic'
	verbose = True
	frequency, fft, max_fft, peak_frequency = select_point(df_FD, extra_points, config)
	max_fft_new, peak_frequency_new = interpolate_peak_frequency(df_FD, frequency, fft, kind, extra_points, interval, config)
	percentage_diff = (peak_frequency_new - peak_frequency) / peak_frequency * 100
	if verbose:
		print ("Interpolation of vel_theta:")
		print ("\t\tBefore interpolation\tAfter interpolation")
		print ("maximum fft:\t{}\t{}".format(max_fft, max_fft_new))
		print ("peak frequency:\t{} Hz\t{} Hz".format(peak_frequency, peak_frequency_new))
		print ("% difference of peak frquency = {} %".format(percentage_diff))
	return peak_frequency_new

#
# =============================================================================
#
#                               Plot Utilities
#
# =============================================================================
#
def plot_all(config):
	rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':10})

	coordinate_system = config["extraction setting"]["coordinate_system"]
	output_data_path = config["IO setting"]["output_data_dir"]
	extraction_method = config["extraction setting"]["extraction_method"]
	frequency_lower_limit = config.getfloat("plot setting", "frequency_lower_limit")
	frequency_upper_limit = config.getfloat("plot setting", "frequency_upper_limit")
	velr_peak_frequency = config.getfloat("plot setting", "velr_peak_frequency") 
	veltheta_peak_frequency = config.getfloat("plot setting", "veltheta_peak_frequency") 
	velphi_peak_frequency = config.getfloat("plot setting", "velphi_peak_frequency")
	plot_avline = config.getboolean("plot setting", "plot_avline")
	#import TD and FD data
	if extraction_method == "single":
		TD_data_path = generate_csv_path(output_data_path, None, "point_data")
		FD_data_path = generate_csv_path(output_data_path, None, "W_vel_FD")

	elif extraction_method == "multiple":
		TD_data_path = generate_csv_path(output_data_path, None, "average_point_data")
		FD_data_path = generate_csv_path(output_data_path, None, "average_FD_data")
	
	df_TD = pd.read_csv(TD_data_path, delimiter = '\t')
	df_FD = pd.read_csv(FD_data_path, delimiter = '\t')
	

	df_FD = filter_FD_data(df_FD, frequency_lower_limit, frequency_upper_limit)

	df_TD['t'] = df_TD['t'] / 2.03001708E02
	#df_TD['rho'] = df_TD['rho'] / 1.6199E-18
	if coordinate_system == "spherical":
		df_FD = df_FD[['frequency',
				'rho_fft', 
				'rho_psd', 
				'vel1_fft', 
				'vel1_psd', 
				'vel2_fft', 
				'vel2_psd', 
				'vel3_fft', 
				'vel3_psd']]
	else:
		df_FD = df_FD[['frequency',
				'rho_fft', 
				'rho_psd', 
				'vel1_fft', 
				'vel1_psd', 
				'vel2_fft', 
				'vel2_psd', 
				'vel3_fft', 
				'vel3_psd',
				'velr_sphe_fft',
				'velr_sphe_psd',
				'veltheta_sphe_fft',
				'veltheta_sphe_psd',
				'velphi_sphe_fft',
				'velphi_sphe_psd',
				'omega_theta_fft']]
		
	plot_color = ["black", "r","g","b","fuchsia"]
	ls = ["-","-.", "--", ":","-."]
	lw = [1,1.5,1.5,1.5,1.5]
	fig = plt.figure(1)
	fig, ((ax0, ax1, ax2, ax3), (ax4, ax5, ax6, ax7)) = plt.subplots(2, 4, figsize=(22, 12))
	ax0.set_ylabel('$\\rho(t)$')
	ax1.set_ylabel('$v^{r}(t)$')
	ax2.set_ylabel('$v^{\\theta}(t)$')
	ax3.set_ylabel('$v^{\\phi}(t)$')
	ax0.set_xlabel('$t$ $\\rm{(ms)}$')
	ax1.set_xlabel('$t$ $\\rm{(ms)}$')
	ax2.set_xlabel('$t$ $\\rm{(ms)}$')
	ax3.set_xlabel('$t$ $\\rm{(ms)}$')
	rho_TD = df_TD['rho']
	if coordinate_system == "spherical":
		v_r = df_TD['vel1']
		v_theta = df_TD['vel2']
		v_phi = df_TD['vel3']
	else:
		v_r = df_TD['velr_sphe']
		v_theta = df_TD['omega_theta']
		v_phi = df_TD['velphi_sphe']
	ax0.plot(df_TD['t'], rho_TD, linestyle = '-' ,linewidth = 1.0 )
	ax1.plot(df_TD['t'], v_r, linestyle = '-' ,linewidth = 1.0 )
	ax2.plot(df_TD['t'], v_theta, linestyle = '-' ,linewidth = 1.0 )
	ax3.plot(df_TD['t'], v_phi, linestyle = '-' ,linewidth = 1.0 )
	rho_TD = df_TD['rho']
	ax4.set_ylabel('$\\rm{FFT\ of\ }$ $\\rho(t)$')
	ax5.set_ylabel('$\\rm{FFT\ of\ }$ $v^{r}(t)$')
	ax6.set_ylabel('$\\rm{FFT\ of\ }$ $v^{\\theta}(t)$')
	ax7.set_ylabel('$\\rm{FFT\ of\ }$ $v^{\\phi}(t)$')
	ax4.set_xlabel('$f$ $\\rm{(Hz)}$')
	ax5.set_xlabel('$f$ $\\rm{(Hz)}$')
	ax6.set_xlabel('$f$ $\\rm{(Hz)}$')
	ax7.set_xlabel('$f$ $\\rm{(Hz)}$')
	ax4.plot(df_FD['frequency'], df_FD['rho_fft'] ,linestyle = '-' ,linewidth = 1.0, label='$\\rm{FFT}$$(\\rho$' )
	if coordinate_system == "spherical": 
		ax5.plot(df_FD['frequency'], df_FD['vel1_fft'] ,linestyle = '-' ,linewidth = 1.0, label='$\\rm{FFT}$$(v^{r})$' )
		ax6.plot(df_FD['frequency'], df_FD['vel2_fft'] ,linestyle = '-' ,linewidth = 1.0, label='$\\rm{FFT}$$(v^{\\theta})$' )
		ax7.plot(df_FD['frequency'], df_FD['vel3_fft'] ,linestyle = '-' ,linewidth = 1.0, label='$\\rm{FFT}$$(v^{\\phi})$' )
	else:
		ax5.plot(df_FD['frequency'], df_FD['velr_sphe_fft'] ,linestyle = '-' ,linewidth = 1.0, label='$\\rm{FFT}$$(v^{r})$' )
		ax6.plot(df_FD['frequency'], df_FD['omega_theta_fft'] ,linestyle = '-' ,linewidth = 1.0, label='$\\rm{FFT}$$(v^{\\theta})$' )
		ax7.plot(df_FD['frequency'], df_FD['velphi_sphe_fft'] ,linestyle = '-' ,linewidth = 1.0, label='$\\rm{FFT}$$(v^{\\phi})$' )
		if plot_avline:
			ax5.axvline(x=velr_peak_frequency, color='g',label = '$^2f$', linestyle='--',linewidth = 0.6)
			ax6.axvline(x=veltheta_peak_frequency, color='g',label = '$^2f$', linestyle='--',linewidth = 0.6)
			ax7.axvline(x=velphi_peak_frequency, color='g',label = '$^2f$', linestyle='--',linewidth = 0.6)
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
	fig.suptitle('density and fluid velocity in time and frequency domain', fontsize=20)
	fig.tight_layout()
	plt.savefig("{}/{}_rho_vel_TD_FD_combined_{}.png".format(output_data_path, coordinate_system, extraction_method), bbox_inches="tight", dpi= 300)
	return 

def workflow_execute(config):
	workflow_extraction = config.getboolean("workflow", "extraction")
	workflow_fft = config.getboolean("workflow", "fft")
	workflow_frequency_peak_extraction = config.getboolean("workflow", "frequency_peak_extraction")
	workflow_plot = config.getboolean("workflow", "plot") 

	if workflow_extraction:
		print ("extracting data ...")
		extraction(config)

	if workflow_fft:
		print ("fft ...")
		FFT_execute(config)

	if workflow_frequency_peak_extraction:
		print ("peak extraction ...")
		extract_peak(config)

	if workflow_plot:
		print ("plotting ...")
		plot_all(config)


#config variables
config = configparser.ConfigParser()
config.read("./extraction.ini")

check_config_value(config)
workflow_execute(config)
#I hate yt very much
