import numpy as np
from scipy.interpolate import RectBivariateSpline

from amrvac_pytools.datfiles.physics import ionisation_tables

# Create bivariate spline approximation over rectangular mesh
# Do this outside of methods, it only has to be done once when importing the module.
spline_ion = None
spline_f = None

def init_splines(altitude):
    """
    Initialises the ionisation and f splines for a certain altitude (in km).
    :param altitude: Height above the solar surface (in km)
    """
    global spline_ion, spline_f
    itable, ftable = ionisation_tables.get_ionisation_table(altitude)

    spline_ion = RectBivariateSpline(ionisation_tables.T_table, ionisation_tables.pg_table, itable)
    spline_f = RectBivariateSpline(ionisation_tables.T_table, ionisation_tables.pg_table, ftable)

def block_interpolate_ionisation_f(block, dataset):
    """
    Interpolates the ionisation and f parameter for every point in a given block.
    :param block: a datfile block as a dictionary, containing the known_variables fields as keys
    :param dataset: instance of the load_datfile class
    :return: ionisation and f parameters as two Numpy arrays of same size as the input block variables
    """
    # Re-dimensionalise temperature and pressure
    temp = block["T"]*dataset.units.unit_temperature
    pg = block["p"]*dataset.units.unit_pressure

    # create empty matrices
    ion = np.zeros_like(temp)
    fpar = np.zeros_like(temp)

    it = np.nditer(temp, flags=['multi_index'])
    while not it.finished:
        idx = it.multi_index
        ion[idx] = spline_ion.ev(temp[idx], pg[idx])
        fpar[idx] = spline_f.ev(temp[idx], pg[idx])
        it.iternext()
    return ion, fpar





