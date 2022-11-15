import numpy as np


#read input file
fin = open("init_struc_total", "rt")
#read file contents to string
data = fin.read()
#replace all occurrences of the required string
data = data.replace('D', 'e')
#close the input file
fin.close()
#open the input file in write mode
fin = open("init_struc_total", "wt")
#overrite the input file with the resulting data
fin.write(data)
#close the file
fin.close()

A = np.loadtxt('init_struc_total',skiprows=1)
i,   rho,   t,   pg,   tau,   mc,   gammar,   er,   r,   v = A.T

np.savetxt('1D_stable.txt',np.transpose([r, rho, v, pg, er]))
