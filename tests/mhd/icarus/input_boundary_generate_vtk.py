#!/usr/bin/env python

import numpy as np
import sys
from astropy import constants as const
import scipy.constants
import matplotlib.pyplot as plt
import glob

from datetime import datetime




######################## User Defined Variables ##########################


dataset =  sorted(glob.glob('./solar_wind_bc_used_in_paper.in'))

filename_boundary_output = './solar_wind_2015_single_gong.vtk'
generated_boundary_images_flag = 'yes'
generated_boundary_images_directory ='./'

######################## Units for AMPRVAC ##########################

unit_length = 6.957e+8
unit_time= 3600.0
unit_velocity = unit_length/unit_time
unit_density = 1.6726e-19
unit_b = np.sqrt( scipy.constants.mu_0 * unit_velocity**2 * unit_density)
unit_pressure = unit_velocity**2 * unit_density


def retrieve_first_last_magnetogram_times(dataset):
    for i, data in enumerate(dataset):
        if i == 0:
            with open(data) as f:
                lines = f.readlines()
            magnetogram_time = lines[1][:-1]
            first_magnetogram = datetime.strptime(magnetogram_time,'%Y-%m-%dT%H:%M:%S')
        if i == len(dataset)-1:
            with open(data) as f:
                lines = f.readlines()
            magnetogram_time = lines[1][:-1]
            last_magnetogram = datetime.strptime(magnetogram_time,'%Y-%m-%dT%H:%M:%S')
    return first_magnetogram, last_magnetogram

first_magnetogram, last_magnetogram = retrieve_first_last_magnetogram_times(dataset)
############ Print infromation for generation of the input .vtk file #############
# first_magnetogram = datetime.strptime(first_magnetogram,'%Y-%m-%dT%H:%M:%S')
# last_magnetogram = datetime.strptime(last_magnetogram,'%Y-%m-%dT%H:%M:%S')
n_time = len(dataset)
if n_time >1:
    average_cadence = (last_magnetogram-first_magnetogram).total_seconds()/3600./(n_time-1)
else:
    average_cadence = (last_magnetogram-first_magnetogram).total_seconds()/3600./(n_time)

print ("Information about the dataset and input boundary file for Icarus")
print ("\n")
print ("Total number of magentograms used for creating a boundary file: ", n_time)
print ("Start of the magnetgorams: ", first_magnetogram)
print ("End of the magnetograms: ", last_magnetogram)
print ("Total time in hours: ", (last_magnetogram-first_magnetogram).total_seconds()/3600.)
print ("Average timestep between magnetograms [hour]: ", str(average_cadence))
print ("\n")
print ("Units for generation the boundary file for Icarus")
print ("\n")
print ("Velocity: ", unit_velocity)
print ("Density: ", unit_density)
print ("Magentic field: ", unit_b)
print ("Pressure: ", unit_pressure)


def retrieve_first_last_magnetogram_times(dataset):
    for i, data in enumerate(dataset):
        if i == 0:
            with open(data) as f:
                lines = f.readlines()
            magnetogram_time = lines[1][:-1]
            first_magnetogram = datetime.strptime(magnetogram_time,'%Y-%m-%dT%H:%M:%S')
        if i == len(dataset)-1:
            with open(data) as f:
                lines = f.readlines()
            magnetogram_time = lines[1][:-1]
            last_magnetogram = datetime.strptime(magnetogram_time,'%Y-%m-%dT%H:%M:%S')
    return first_magnetogram, last_magnetogram


def shift_data_by_longitude(data, longitude_shift, lons):
    longitude_resolution = 1/len(lons)
    # Determine the index offset
    lon_shift_index = int(longitude_shift / longitude_resolution)  # Assuming a constant resolution

    # Create a new array for shifted data
    shifted_data = [[0] * len(data[0]) for _ in range(len(data))]

    # Shift data by longitude
    for i in range(len(data)):
        for j in range(len(data[i])):
            new_j = (j + lon_shift_index) % len(data[i])  # Circular shifting
            shifted_data[i][new_j] = data[i][j]

    return shifted_data


density_all = []
vr_all = []
pressure_all = []
br_all = []
prev_magnetogram = datetime(2999,1,1)
current_magnetgoram = datetime.now()
shift_accumulated = 0.

for i, data in enumerate(dataset):

    with open(data) as f:
        lines = f.readlines()

    magnetogram_time = lines[1][:-1]
    if (prev_magnetogram == datetime(2999,1,1)):
        prev_magnetogram = datetime.strptime(magnetogram_time,'%Y-%m-%dT%H:%M:%S')
        current_magnetgoram = datetime.strptime(magnetogram_time,'%Y-%m-%dT%H:%M:%S')
    if (magnetogram_time != current_magnetgoram):
        current_magnetgoram = datetime.strptime(magnetogram_time,'%Y-%m-%dT%H:%M:%S')
    time_step = (current_magnetgoram-prev_magnetogram).total_seconds()/3600.
    prev_magnetogram = datetime.strptime(magnetogram_time,'%Y-%m-%dT%H:%M:%S')

    radius_inner_bc = float(lines[3][:-1])
    n_clt_points = int(lines[5][:-1])
    lon_num_index = 5+1+n_clt_points+1+1
    n_lon_points = int(lines[lon_num_index][:-1])

    for j, line in enumerate(lines):
    # Check if "vr\n" is in the line
        if "vr\n" in line:
            vr_index = j
        elif "number_density\n" in line:
            n_index=j
        elif "temperature\n" in line:
            temp_index=j
        elif "Br\n" in line:
            br_index=j


    vr_string = lon_num_index+1+n_lon_points+1
    clts = np.array(lines[7:lon_num_index-1])
    lons = np.array(lines[lon_num_index+2:vr_string])
    if (n_index > br_index):
        vr_array = np.array(lines[vr_index+1:br_index])
        num_density_array = np.array(lines[n_index+1:temp_index])
        temperature_array = np.array(lines[temp_index+1:])
        br_array = np.array(lines[br_index+1:n_index])
    else:
        vr_array = np.array(lines[vr_index+1:n_index])
        num_density_array = np.array(lines[n_index+1:temp_index])
        temperature_array = np.array(lines[temp_index+1:br_index])
        br_array = np.array(lines[br_index+1:])


    lons_original = lons.astype(np.float64)
    clts_original=clts.astype(np.float64)

    print (time_step, current_magnetgoram)
    clts = np.linspace(1.745329251994309772e-02, 3.124139361069850018e+00, len(clts))
    clt_min = 8.14606741573e-2*np.pi*2.0
    clt_max = 4.18539325843e-1*np.pi*2.0
    lon_correct = abs(float(lons_original[0]))
    lons = np.linspace(float(lons_original[0]) + lon_correct, float(lons_original[-1])+lon_correct, len(lons))
    delta_lon=-time_step*2.*np.pi/(24.47*24)
    shift_accumulated+=delta_lon


    vr_array = vr_array.astype(np.float64)
    num_density_array = num_density_array.astype(np.float64)
    temperature_array = temperature_array.astype(np.float64)
    pressure_array = scipy.constants.k*num_density_array*temperature_array
    br_array =br_array.astype(np.float64)
    num_density_array = num_density_array.astype(np.float64)* 1.6726217770000001e-027 * 0.5

    if (min(pressure_array) < 0):
        print ("There is negative pressure in " + data)

    vr_nondim = vr_array/unit_velocity
    n_nondim = num_density_array/unit_density
    p_nondim = pressure_array/unit_pressure
    br_nondim = br_array/unit_b

    generated_lons = []
    generated_lats = []
    temp_lons = []
    temp_lats = []

    for lon_1 in lons:
        for clt_1 in clts:
            temp_lons.append(lon_1)
            temp_lats.append(clt_1)


    for clt_0 in clts:
        for lon_0 in lons:
            generated_lons.append(lon_0)
            generated_lats.append(clt_0)



    generated_lats=np.array(generated_lats)
    generated_lons=np.array(generated_lons)


    nlon=n_lon_points
    nlat=n_clt_points
    lon=generated_lons.reshape((nlat,nlon))
    lat=generated_lats.reshape((nlat,nlon))

    vr=vr_nondim.reshape((nlon,nlat)).T
    n_dens=n_nondim.reshape((nlon,nlat)).T
    pressure=p_nondim.reshape((nlon,nlat)).T
    br=br_nondim.reshape((nlon,nlat)).T

    longitude_resolution = (max(lons)-min(lons))/len(lons)
    num_cols_shift = int(shift_accumulated / longitude_resolution)

    vr = np.roll(vr, num_cols_shift, axis=1)
    n_dens = np.roll(n_dens, num_cols_shift, axis=1)
    pressure = np.roll(pressure, num_cols_shift, axis=1)
    br = np.roll(br, num_cols_shift, axis=1)


    vr_all.append(vr.transpose().flatten())
    density_all.append(n_dens.transpose().flatten())
    pressure_all.append(pressure.transpose().flatten())
    br_all.append(br.transpose().flatten())



    if (generated_boundary_images_flag == 'yes'):
        plt.pcolormesh(lon, lat, vr, cmap='rainbow')
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(8,10),sharex=True)

        im1 = ax1.pcolormesh(lon, lat, np.flipud(vr*unit_velocity), cmap='jet', shading='auto')
        fig.colorbar(im1, ax=ax1)
        ax1.set_ylabel(r'$V_r$')
        ax1.yaxis.set_label_position("right")

        im2 = ax2.pcolormesh(lon, lat, np.flipud(n_dens*unit_density), cmap='hot_r', shading='auto')
        fig.colorbar(im2, ax=ax2)
        ax2.set_ylabel(r'Number Desnity')
        ax2.yaxis.set_label_position("right")

        im3 = ax3.pcolormesh(lon, lat, np.flipud(pressure*unit_pressure* 1.6726217770000001e-027 * 0.5)/(np.flipud(n_dens*unit_density)*scipy.constants.k), cmap='hot', shading='auto')
        fig.colorbar(im3, ax=ax3)
        ax3.set_ylabel(r'Temperture')
        ax3.yaxis.set_label_position("right")

        im4 = ax4.pcolormesh(lon, lat, np.flipud(br*unit_b), cmap='seismic', shading='auto')
        fig.colorbar(im4, ax=ax4)
        ax4.set_ylabel(r'$B_r$')
        ax4.yaxis.set_label_position("right")

        if (i < 10):
            file_name = '0000'+str(i)
        elif (i < 100):
            file_name = '000'+str(i)
        elif (i < 1000):
            file_name = '00'+str(i)
        elif (i < 10000):
            file_name = '0' +str(i)
        plt.savefig(generated_boundary_images_directory+"boundary_"+file_name)

    if (i == 0):

        dlon = (lon[0,1]-lon[0,0])
        dlat = (lat[1,0]-lat[0,0])



        left=generated_lons[0]
        right=generated_lons[-1]
        bottom=generated_lats[0]
        top=generated_lats[-1]

vtk_header="""\
# vtk DataFile Version 2.0
File generated by picture_to_vtk.py script
ASCII
DATASET STRUCTURED_POINTS
DIMENSIONS {} {} {}
SPACING {} {} {}
ORIGIN {} {}  0
POINT_DATA {}""".format(nlat,nlon,n_time,dlat,dlon,average_cadence, bottom, left, nlat*nlon*n_time)

file1 = open(filename_boundary_output, "w")
file1.write(vtk_header)
headerdataset="""
SCALARS bc_data float 1
LOOKUP_TABLE default
"""
vr_all = np.concatenate((vr_all[:]))
density_all = np.concatenate((density_all[:]))
pressure_all = np.concatenate((pressure_all[:]))
br_all = np.concatenate((br_all[:]))


##WRITE number density
file1.write(headerdataset)
file1.write(" ".join(map(lambda x : str(x), density_all.transpose().flatten())))


#Write vr
file1.write(headerdataset)
file1.write(" ".join(map(lambda x : str(x), vr_all.transpose().flatten())))

##WRITE pressure
file1.write(headerdataset)
file1.write(" ".join(map(lambda x : str(x), pressure_all.transpose().flatten())))


##WRITE Br
file1.write(headerdataset)
file1.write(" ".join(map(lambda x : str(x), br_all.transpose().flatten())))




file1.close()
print ("The delta_phi variable to modify in the .par file is: ", float(lons_original[0]))
print ("The generated boundary file for Icarus is saved tos ", filename_boundary_output)
if (generated_boundary_images_flag == 'yes'):
    print ("The generated boundary images are stored in ", generated_boundary_images_directory[:-1])
