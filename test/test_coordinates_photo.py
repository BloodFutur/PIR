# Import required libraries
import unittest

import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/../src")

# Import the functions to test
import coordinates_photo as cp

class TestCameraCalibrationWithAstrometryNet(unittest.TestCase):
    def test_calibrate_camera(self):
        photo = "images/test/test_12_mars_21h46mn36s.jpg"
        reference_pixel, reference_point, transformation_matrix = cp.calibrate_camera(photo)
        self.assertEqual(reference_pixel, (955.153951009, 437.333496094))
        self.assertEqual(reference_point, (79.8535963865, 10.0878946534))
        self.assertEqual(transformation_matrix, ([-0.0358029908004, -0.027435406908],
                                                 [0.0270064129609, -0.034174545714]))

class TestCameraCalibrationWithoutAstrometryNet(unittest.TestCase):
    def test_get_calibrated_camera_parameters(self):
        photo = "images/test/test_12_mars_21h46mn36s.jpg"
        reference_pixel, reference_point, transformation_matrix = cp.get_calibrated_camera_parameters('output/test/test_wcs_header.txt')
        self.assertEqual(reference_pixel, (955.153951009, 437.333496094))
        self.assertEqual(reference_point, (79.8535963865, 10.0878946534))
        self.assertEqual(transformation_matrix, ([-0.0358029908004, -0.027435406908],
                                                 [0.0270064129609, -0.034174545714]))
# Example test case
class TestCoordinatesPhoto(unittest.TestCase):
    def test_convert_coordinates(self):
        observer_location = (35.6895, 139.6917, 0)
        observer_time = "2021-06-01 12:00:00"
        object_coordinates = (0, 0)
        object_location = cp.convert_coordinates(observer_location, observer_time, object_coordinates)
        self.assertEqual(object_location, (35.6895, 139.6917))

if __name__ == "__main__":
    unittest.main()