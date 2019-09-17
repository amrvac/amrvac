import numpy as np
from scipy.interpolate import RectBivariateSpline

from amrvac_tools.datfiles.physics import ionisation_tables

# Create bivariate spline approximation over rectangular mesh
# Do this outside of methods, it only has to be done once when importing the module.
spline_ion = None
spline_f = None

def init_splines(altitude):
    global spline_ion, spline_f
    itable, ftable = ionisation_tables.get_ionisation_table(altitude)

    spline_ion = RectBivariateSpline(ionisation_tables.T_table, ionisation_tables.pg_table, itable)
    spline_f = RectBivariateSpline(ionisation_tables.T_table, ionisation_tables.pg_table, ftable)

def block_interpolate_ionisation_f(block, block_fields, dataset):
    # Re-dimensionalise temperature and pressure. Use Ellipsis object (...) to fill in missing dimensions
    temp = block[..., block_fields.index("T")]*dataset.units.unit_temperature
    pg = block[..., block_fields.index("p")]*dataset.units.unit_pressure

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





