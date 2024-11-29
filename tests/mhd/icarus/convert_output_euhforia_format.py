import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import interp1d
import scipy.constants
import glob
import dateutil.parser
from datetime import datetime, timedelta
import csv
import sys

# Usage : python3 convert_output_euhforia_format.py forecast 14 0 2015-06-25T01:04:00


trajectory_data = sorted(glob.glob('./icarus_testcase_particle_00000?.csv'))


##### Reading in the passed arguments
arguments = sys.argv
readout_start = str(arguments[1])
relaxation = float(arguments[2])
cme_insertion = int(arguments[3])
magnetogram_timestamp = datetime.strptime(arguments[4], '%Y-%m-%dT%H:%M:%S')
timestamp = relaxation+cme_insertion
r_earth = []


save_from_hrs = (relaxation+cme_insertion)*24.0
save_from_index = 0
if (readout_start == 'forecast'):
    save_from_hrs = (relaxation+cme_insertion)*24.0
elif (readout_start == 'cme_insertion'):
    save_from_hrs = relaxation*24.0
elif (readout_start == 'relaxation'):
    save_from_hrs = 0.0

start_date = magnetogram_timestamp-timedelta(days=(relaxation+cme_insertion))
# start_date = magnetogram_timestamp-timedelta(days=(cme_insertion))


##### Constants used in the script
au_in_solar_radius = 215.032
r_sun = 6.957e+8
omega = 2.97e-6
m_p = scipy.constants.proton_mass


def get_v_lon(r, v_lon):
    v_lon = v_lon + r * r_sun* omega/100.0
    return v_lon

for i, data in enumerate(trajectory_data):
    csv_file = open(data, 'r')
    file = csv.DictReader(csv_file)
    date_csv1 = []
    r_csv1 = []
    clt_csv1 = []
    lon_csv1 = []
    lon_heeq = []
    n_csv1 = []
    p_csv1 = []
    vr_csv1 = []
    v_clt_csv1 = []
    v_lon_csv1 = []
    v_lon_heeq = []
    br_csv1 = []
    b_clt_csv1 = []
    b_lon_csv1 = []


    magnetogram_time= timestamp*24.0
    magnetogram_index = 0


    for col in file:
        date_csv1.append(float(col[' time']))
        r_csv1.append(float(col[' x1'])/au_in_solar_radius)
        clt_csv1.append(float(col[' x2']))
        lon_csv1.append(float(col[' x3']))
        n_csv1.append(float(col[' rho'])/((10**6)*0.5*m_p))
        p_csv1.append(float(col[' p']))
        vr_csv1.append(float(col[' v1']))
        v_clt_csv1.append(float(col[' v2']))
        v_lon_csv1.append(float(col[' v3']))
        br_csv1.append(float(col[' b1']))
        b_clt_csv1.append(float(col[' b2']))
        b_lon_csv1.append(float(col[' b3']))

    if (i == i):
        r_earth = r_csv1

    for j in range(len(date_csv1)):
        if (date_csv1[j] >= magnetogram_time):
            magnetogram_index = j
            break

    for j in range(len(date_csv1)):
        if (date_csv1[j] >= save_from_hrs):
            save_from_index = j
            break
    for k in range(len(date_csv1)):
        v_lon_heeq.append(get_v_lon(r_earth[k], v_lon_csv1[k]))
        lon_heeq.append(np.mod(lon_csv1[k]+omega*timedelta(hours=(date_csv1[k]-date_csv1[magnetogram_index])).total_seconds(), 2*np.pi))


    ###### Generate time array from simulation times
    freq = (timedelta(hours=(date_csv1[-1]-date_csv1[0]))).total_seconds()/len(date_csv1)
    time_array = np.arange(start_date, start_date+timedelta(hours=(date_csv1[-1]-date_csv1[0])), timedelta(seconds=freq)).astype(datetime)
    for n in range(len(time_array)):
        time_array[n] = time_array[n].strftime('%Y-%m-%dT%H:%M:%S')

    ####### Generate set of lists to write to the output file
    data_to_write = list(zip(time_array[save_from_index:], r_csv1[save_from_index:], clt_csv1[save_from_index:], lon_csv1[save_from_index:],n_csv1[save_from_index:], p_csv1[save_from_index:], vr_csv1[save_from_index:], v_clt_csv1[save_from_index:], v_lon_heeq[save_from_index:], br_csv1[save_from_index:], b_clt_csv1[save_from_index:], b_lon_csv1[save_from_index:]))


    if (int(data[-5]) == 1) :
        file_path = data[:-20]+'_Earth.dsv'
    if (int(data[-5]) == 2) :
        file_path = data[:-20]+'_Mars.dsv'
    if (int(data[-5]) == 3) :
        file_path = data[:-20]+'_Mercury.dsv'
    if (int(data[-5]) == 4) :
        file_path = data[:-20]+'_Venus.dsv'
    if (int(data[-5]) == 5) :
        file_path = data[:-20]+'_STA.dsv'
    if (int(data[-5]) == 6) :
        file_path = data[:-20]+'_STB.dsv'
    if (int(data[-5]) == 7) :
        file_path = data[:-20]+'_PSP.dsv'
    if (int(data[-5]) == 8) :
        file_path = data[:-20]+'_SolO.dsv'


    # Write the array to a text file
    with open(file_path, 'w') as file:
        file.write('date r[AU] clt[rad] lon[rad] n[1/cm^3] P[Pa] vr[km/s] vclt[km/s] vlon[km/s] Br[nT] Bclt[nT] Blon[nT]'+'\n')
        for row in data_to_write:
            file.write(' '.join([str(element) for element in row]) + '\n')
