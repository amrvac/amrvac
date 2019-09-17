"""
Module to calculate the (M)HD energies and conversions between primitive and conservative variables.
"""

# ====== HYDRODYNAMICS ======
def _hd_ekin(rho, mom1, mom2, mom3):
    return 0.5 * (mom1 ** 2 + mom2 ** 2 + mom3 ** 2) / rho


def hd_to_primitive_1d(rho, mom1, e, gamma):
    v1 = mom1 / rho
    p = (gamma - 1) * (e - _hd_ekin(rho, mom1, 0, 0))
    return v1, p


def hd_to_conserved_1d(rho, v1, p, gamma):
    mom1 = v1 * rho
    e    = p / (gamma - 1) + 0.5 * v1**2 * rho
    return mom1, e


def hd_to_primitive_2d(rho, mom1, mom2, e, gamma):
    v1 = mom1 / rho
    v2 = mom2 / rho
    p = (gamma - 1) * (e - _hd_ekin(rho, mom1, mom2, 0))
    return v1, v2, p


def hd_to_conserved_2d(rho, v1, v2, p, gamma):
    mom1 = v1 * rho
    mom2 = v2 * rho
    e    = p / (gamma - 1) + 0.5 * (v1**2 + v2**2) * rho
    return mom1, mom2, e


def hd_to_primitive_3d(rho, mom1, mom2, mom3, e, gamma):
    v1 = mom1 / rho
    v2 = mom2 / rho
    v3 = mom3 / rho
    p = (gamma - 1) * (e - _hd_ekin(rho, mom1, mom2, mom3))
    return v1, v2, v3, p


def hd_to_conserved_3d(rho, v1, v2, v3, p, gamma):
    mom1 = v1 * rho
    mom2 = v2 * rho
    mom3 = v3 * rho
    e    = p / (gamma - 1) + 0.5 * (v1**2 + v2**2 + v3**3) * rho
    return mom1, mom2, mom3, e


# ====== MAGNETOHYDRODYNAMICS ======
def _mhd_ekin(rho, mom1, mom2, mom3):
    return 0.5 * (mom1 ** 2 + mom2 ** 2 + mom3 ** 2) / rho


def _mhd_emag(b1, b2, b3):
    return 0.5 * (b1 ** 2 + b2 ** 2 + b3 ** 2)


def mhd_to_primitive_1d(rho, mom1, e, b1, gamma):
    v1 = mom1 / rho
    p = (gamma - 1) * (e - _mhd_ekin(rho, mom1, 0, 0) - _mhd_emag(b1, 0, 0))
    return v1, p


def mhd_to_conserved_1d(rho, v1, p, b1, gamma):
    mom1 = v1 * rho
    e    = p / (gamma - 1) + 0.5 * v1**2 * rho + _mhd_emag(b1, 0, 0)
    return mom1, e


def mhd_to_primitive_2d(rho, mom1, mom2, e, b1, b2, gamma):
    v1 = mom1 / rho
    v2 = mom2 / rho
    p = (gamma - 1) * (e - _mhd_ekin(rho, mom1, mom2, 0) - _mhd_emag(b1, b2, 0))
    return v1, v2, p


def mhd_to_conserved_2d(rho, v1, v2, p, b1, b2, gamma):
    mom1 = v1 * rho
    mom2 = v2 * rho
    e    = p / (gamma - 1) + 0.5 * (v1**2 + v2**2) * rho + _mhd_emag(b1, b2, 0)
    return mom1, mom2, e


def mhd_to_primitive_3d(rho, mom1, mom2, mom3, e, b1, b2, b3, gamma):
    v1 = mom1 / rho
    v2 = mom2 / rho
    v3 = mom3 / rho
    p = (gamma - 1) * (e - _mhd_ekin(rho, mom1, mom2, mom3) - _mhd_emag(b1, b2, b3))
    return v1, v2, v3, p


def mhd_to_conserved_3d(rho, v1, v2, v3, p, b1, b2, b3, gamma):
    mom1 = v1 * rho
    mom2 = v2 * rho
    mom3 = v3 * rho
    e    = p / (gamma - 1) + 0.5 * (v1**2 + v2**2 + v3**2) * rho + _mhd_emag(b1, b2, b3)
    return mom1, mom2, mom3, e