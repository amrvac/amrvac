from amrvac_pytools.vtkfiles import read, amrplot
import numpy as np
import glob
import pandas as pd

#extraction location
x_cut = 3.0
y_cut = 3.0
z_cut = 0.0

data_sets = glob.glob('./output/output_d3_x+0.00D+00_n*.vtu')

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

for i, data_sets in enumerate(data_sets):
 # Load the content
 ds = read.load_vtkfile(i, file='output_d3_x+0.00D+00_n', type='vtu')	
 rho = ds.rho
 psi = ds.psi
 alp = ds.alp
 W_vel1 = ds.W_vel1
 W_vel2 = ds.W_vel2
 W_vel3 = ds.W_vel3
 coord = ds.getCenterPoints()

 # locate where is the data
 diff_min = 1000
 for j in range(len(coord)):     
  x_j = coord[j,0]
  y_j = coord[j,1]
  diff = np.sqrt((x_cut - x_j)**2 + (y_cut - y_j)**2)
  if diff < diff_min:
   diff_min = diff
   data_loc = j     

 # fill the arrays 
 time = np.append(time, ds.time) 
 rho_p = np.append(rho_p, rho[data_loc])  
 psi_p = np.append(psi_p, psi[data_loc])  
 alp_p = np.append(alp_p, alp[data_loc])  
 W_vel1_p = np.append(W_vel1_p, W_vel1[data_loc])  
 W_vel2_p = np.append(W_vel2_p, W_vel2[data_loc])  
 W_vel3_p = np.append(W_vel3_p, W_vel3[data_loc])  
 
 # now workout veloc1 to 3 
 wv1 = float(W_vel1[data_loc])
 wv2 = float(W_vel2[data_loc])
 wv3 = float(W_vel3[data_loc]) 
 psi4 = float(psi[data_loc]**4)
 w = wv1**2 + wv2**2 + wv3**3
 w = psi4 * w
 w = np.sqrt( 1.0 + w )
 vel1_p = np.append(vel1_p, W_vel1[data_loc]/w)
 vel2_p = np.append(vel2_p, W_vel2[data_loc]/w)
 vel3_p = np.append(vel3_p, W_vel3[data_loc]/w)


# Sort the values by 't'
srt = np.argsort(time)

df_TD = pd.DataFrame({'t' : time[srt], 
                      'rho' : rho_p[srt],
                      'psi' : psi_p[srt],
                      'alp' : alp_p[srt],
                      'W_vel1' : W_vel1_p[srt],
                      'W_vel2' : W_vel2_p[srt],
                      'W_vel3' : W_vel3_p[srt],
                      'vel1' : vel1_p[srt],
                      'vel2' : vel2_p[srt],
                      'vel3' : vel3_p[srt]   }  )

df_TD = df_TD[['t','rho','psi','alp','W_vel1','W_vel2','W_vel3','vel1','vel2','vel3']]
df_TD.to_csv("point_data.csv", index = False, sep = '\t')
