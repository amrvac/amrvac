"""
Created on 5 Dec 2018

Module for conversion between conservative and primitive variables.
MPI-AMRVAC 2.0 .dat snapshots contain conservative data.

@author: Niels Claes
"""

# ===== HYDRODYNAMICS =====
def _hd_ekin(rho, momx, momy, momz):
    return 0.5 * (momx**2 + momy**2 + momz**2) / rho

def hd_to_primitive_1d(rho, momx, e, gamma):
    p  = (gamma - 1) * (e - _hd_ekin(rho, momx, 0, 0))
    v1 = momx / rho
    return p, v1

def hd_to_primitive_2d(rho, momx, momy, e, gamma):
    p  = (gamma - 1) * (e - _hd_ekin(rho, momx, momy, 0))
    v1 = momx / rho
    v2 = momy / rho
    return p, v1, v2

def hd_to_primitive_3d(rho, momx, momy, momz, e, gamma):
    p  = (gamma - 1) * (e - _hd_ekin(rho, momx, momy, momz))
    v1 = momx / rho
    v2 = momy / rho
    v3 = momz / rho
    return p, v1, v2, v3
    
# ===== MAGNETOHYDRODYNAMICS =====
def _mhd_ekin(rho, momx, momy, momz):
    return 0.5 * (momx**2 + momy**2 + momz**2) / rho

def _mhd_emag(b1, b2, b3):
    return 0.5 * (b1**2 + b2**2 + b3**2)

def mhd_to_primitive_1d(rho, momx, e, b1, gamma):
    p  = (gamma - 1) * (e - _mhd_ekin(rho, momx, 0, 0) - _mhd_emag(b1, 0, 0))
    v1 = momx / rho
    return p, v1

def mhd_to_primitive_2d(rho, momx, momy, e, b1, b2, gamma):
    p  = (gamma - 1) * (e - _mhd_ekin(rho, momx, momy, 0) - _mhd_emag(b1, b2, 0))
    v1 = momx / rho
    v2 = momy / rho
    return p, v1, v2

def mhd_to_primitive_3d(rho, momx, momy, momz, e, b1, b2, b3, gamma):
    p  = (gamma - 1) * (e - _mhd_ekin(rho, momx, momy, momz) - _mhd_emag(b1, b2, b3))
    v1 = momx / rho
    v2 = momy / rho
    v3 = momz / rho
    return p, v1, v2, v3
    
    


