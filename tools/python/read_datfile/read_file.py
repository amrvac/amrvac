
import reduce_data
import argparse
import time


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='.dat file from MPI-AMRVAC simulation')
    parser.add_argument('-d', '--dat_file', type=str, help='MPI-AMRVAC .dat file to read')
    args = parser.parse_args()

    time_start = time.time()

    filename = "PATH_TO_FILE"

    # Check if .dat file is given as command line argument
    if args.dat_file is not None:
        filename = args.dat_file

    data = reduce_data.ProcessData(filename=filename)

    temperature = data.T

    # Now call the required data, for temperature use eg. data.T
    # Possible calls are rho, momx, momy, momz, e, b1, b2, b3, p, v1, v2, v3, T
    # Load previously regridded files by doing data = numpy.load("path_to_data.npy")
