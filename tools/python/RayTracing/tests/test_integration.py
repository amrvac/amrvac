'''
Created on 10 Dec 2018

@author: Niels
'''
import unittest
import numpy as np
from physics import opacity

class DataObject():
    def __init__(self):
        self._x = np.arange(3)
        self._y = np.arange(3)
        self._z = np.arange(3)
        self._ndim = 3


class Test(unittest.TestCase):
    
    def get_test_array(self):
        a = np.asarray([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        return np.asarray([a, a, a])
        
    
    def test_integration_x(self):
        data = DataObject()
        matrix_test = self.get_test_array()
        
        expected = np.asarray([[2,4,6],[8,10,12], [14,16,18]])
        integrated = opacity.integrate_los(data, matrix_test, "x")
        
        np.testing.assert_array_almost_equal(integrated, expected, decimal=6)
        
    def test_integration_y(self):
        data = DataObject()
        matrix_test = self.get_test_array()
        
        expected = np.asarray([[8,10,12], [8,10,12], [8,10,12]])
        integrated = opacity.integrate_los(data, matrix_test, "y")
        
        np.testing.assert_array_almost_equal(integrated, expected, decimal=6)
        
    def test_integration_z(self):
        data = DataObject()
        matrix_test = self.get_test_array()
        
        expected = np.asarray([[4,10,16], [4,10,16], [4,10,16]])
        integrated = opacity.integrate_los(data, matrix_test, "z")
        
        np.testing.assert_array_almost_equal(integrated, expected, decimal=6)


if __name__ == "__main__":
    unittest.main()