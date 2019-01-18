"""
Module to integrate a given 3D matrix along an axis.

Created on 7 Dec 2018

@author: Niels Claes
"""

import numpy as np

def integrate_los(data, matrix, axis):
    """
    Method to integrate the absorption coefficient matrix along a given axis.
    @param data: Reduced data object.
    @param matrix: Matrix to integrate.
                   Type is np.ndarray of dimension 3
    @param axis:  (String) Axis along which kappa is integrated, this is the line of sight.
                  Default is x.
    @return: Integrated matrix along the required line of sight.
             Type is np.ndarray of dimension 2
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
        tau = np.zeros_like(matrix[0, :, :])
        #Iterate over indices of y-z plane
        for i, j in np.ndindex(tau.shape):
            col = matrix[:, i, j]
            integrated_col = np.trapz(col, data._x)
            coords = (i, j)
            tau[coords] = integrated_col
        return tau
    elif axis == "y":
        tau = np.zeros_like(matrix[:, 0, :])
        for i, j in np.ndindex(tau.shape):
            col = matrix[i, :, j]
            integrated_col = np.trapz(col, data._y)
            coords = (i, j)
            tau[coords] = integrated_col
        return tau
    else:
        tau   = np.zeros_like(matrix[:, :, 0])
        for i, j in np.ndindex(tau.shape):
            col = matrix[i, j, :]
            integrated_col = np.trapz(col, data._z)
            coords = (i, j)
            tau[coords] = integrated_col
        return tau