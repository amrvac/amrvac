"""
Created on 7 Dec 2018

Module to integrate the absorption coefficient kappa along the line of sight.
Works in 1d, 2d and 3d, along any of the given coordinate axes.
@note: At the moment, only integration along one of the coordinate axes is supported.

@author: Niels Claes
"""

from physics import absorption_coeff
import numpy as np

def integrate_los(data, kappa, axis):
    """
    Method to integrate the absorption coefficient matrix along a given axis.
    @param data: Reduced data object.
    @param kappa: The absorption coefficient.
                  Type is np.ndarray of dimension ndim
    @param axis:  (String) Axis along which kappa is integrated, this is the line of sight.
                  Default is x.
    @return: Integrated absorption coefficient (opacity) along the required line of sight.
             Type is np.ndarray of dimension ndim
    @note: For a 3D matrix in numpy, axes are labelled as:
           - axis = 0 means depth  (x)
           - axis = 1 means height (y), (0,0) is top left corner
           - axis = 2 means width  (z)
           This can be checked through
           a = np.ones((3,5,7))
           x0 = np.sum(a,axis=0)
           x1 = np.sum(a,axis=1)
           x2 = np.sum(a,axis=2)
    """
    print("Integrating matrix along line of sight (%s-axis)" % axis)
    if axis == "x":
        tau = np.zeros_like(kappa[0, :, :])
        #Iterate over indices of y-z plane
        for i, j in np.ndindex(tau.shape):
            col = kappa[:, i, j]
            integrated_col = np.trapz(col, data._x)
            coords = (i, j)
            tau[coords] = integrated_col
        return tau
    elif axis == "y":
        tau = np.zeros_like(kappa[:, 0, :])
        for i, j in np.ndindex(tau.shape):
            col = kappa[i, :, j]
            integrated_col = np.trapz(col, data._y)
            coords = (i, j)
            tau[coords] = integrated_col
        return tau
    else:
        tau   = np.zeros_like(kappa[:, :, 0])
        for i, j in np.ndindex(tau.shape):
            col = kappa[i, j, :]
            integrated_col = np.trapz(col, data._z)
            coords = (i, j)
            tau[coords] = integrated_col
        return tau
    

def get_opacity(data, axis):
    """
    Returns the opacity along the line of sight.
    @param data: Reduced data object
    @param axis: (String) Line of sight, default = x
    @return: - Opacity interpolated from input if ndim < 3
             - Integrated opacity along line of sight if ndim == 3
             Type is np.ndarray of dimension 2
    """
    #Retrieve absorption coefficient
    kappa = absorption_coeff.get_kappa(data)
    
    #If 1D or 2D, simply return the grid itself
    if not data._ndim == 3:
        return kappa
    else:
        return integrate_los(data, kappa, axis)
        
        
