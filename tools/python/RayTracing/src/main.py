"""
Program with raytracing algorithms for a given MPI-AMRVAC .dat file.
Currently implemented functions: H-alpha views and the Faraday rotation effect.
Basic steps:
    1) Read in .dat file
    2) Check if mesh is uniformely refined, if not perform regridding.
       Parallelization across multiple cores can be used for large data sets.
    3) Interpolate the ionization degree and parameter f across the entire mesh using pressure and temperature.
       This is used to derive the electron density at each grid point.
    4) Save the interpolated matrices as .npy files (Numpy files), as for large matrices the interpolation can
       take a while. This way the data can be loaded easily if a re-run is required.
    5) Use the ionization degree and parameter f to calculate required parameters across the mesh
    6) Integrate the mesh along a given line of sight (x, y or z axis) in 3D, or, in 2D,
       simply return the mesh.
    7) Plot the required results.
    
Usage: 2 options
    1) Change the filename to use in settings.py and run from a compiler or command line without arguments
    2) Run from command line via python main.py -d DAT_FILE -ns -np x
       with:
         - d DAT_FILE: the .dat file to read in
         - ns: If this option is passed, the code will NOT save the interpolated data and/or
               the regridded data as .npy (Numpy) files. Every re-run will require a new calculation from scratch.
         - np: Optional, if passed Python's multiprocessing module will be activated to perform the regridding.
               The argument that follows equals the number of CPUs to use, eg. -np 6 will use 6 cores.
       If no arguments are given, everything defaults to the values in settings.py
       
REMARK: Before using, change the parameters in settings.py to your own desired values.
        Make sure the normalizations are correct, as unit_pressure and unit_temperature are used
        to re-dimensionalize the variables for a correct calculation of the H-alpha view.

Created on 4 Dec 2018

Last update on 09 Jan 2019

@author: Niels Claes
@contact: niels.claes@kuleuven.be
"""

from dataIO.reduce_data import ProcessData
from views import h_alpha_view, faraday_view
import settings
import os, time
import matplotlib.pyplot as plt
import argparse
import multiprocessing


if __name__ == '__main__':
    # ====================================================
    # ==============  COMMAND LINE ARGUMENTS  ============
    # ====================================================
    parser = argparse.ArgumentParser(description='.dat file from MPI-AMRVAC simulation')
    parser.add_argument('-d', '--dat_file', type=argparse.FileType('rb'), help='MPI-AMRVAC .dat file to read')
    parser.add_argument('-ns', '--disable_save', default=False, action='store_true', 
                        help="If given, disables saving of interpolated and/or regridded files.")
    parser.add_argument('-np', '--number_of_procs', type=int, help='Setting a number of processors different from'
                                                                   'one activates parallelization of computationally'
                                                                   'expensive tasks (interpolation and regridding).')
    args = parser.parse_args()

    # ====================================================
    # ===============   INITIALIZATIONS   ================
    # ====================================================
    # Check if .dat file is given as command line argument    
    if args.dat_file is not None:
        dat_file = args.dat_file
        settings.filename = dat_file.name
    else:
        dat_file = open(settings.filename, 'rb')
        
    # Check for optional save argument
    if args.disable_save:
        settings.saveFiles = False
        print("Saving of files is disabled.")
    else:
        # Create save directory
        if not os.path.isdir("interpolated_files"):
            os.mkdir("interpolated_files")

    # Check if multiprocessing needs to be activated
    if args.number_of_procs is not None:
        if not args.number_of_procs == 1:
            print("entered")
            settings.multiple_procs = True
            # Check if supplied numer is greater than available processors
            if args.number_of_procs > multiprocessing.cpu_count():
                settings.nb_of_procs = multiprocessing.cpu_count()
            else:
                settings.nb_of_procs = args.number_of_procs
            print("Multiprocessing enabled.")
    time_start = time.time()




    # ====================================================
    # ================   PROCESSING   ====================
    # ====================================================
    # Create data object
    data = ProcessData(dat_file)




    # ====================================================
    # =================   PLOTTING   =====================
    # ====================================================
    if settings.halpha:
        h_alpha_view.plot_h_alpha(data)
    if settings.faraday:
        faraday_view.plot_faraday(data)


    time_end = time.time() - time_start
    print("Total running time: %.0f s" % time_end)

    plt.show()
    
    print("Done.")


    
    
    
