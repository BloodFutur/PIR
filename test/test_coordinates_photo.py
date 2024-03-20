# Import required libraries
import unittest

import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/../src")

# Import the functions to test
import coordinates_photo as cp

# Example test case
class TestCoordinatesPhoto(unittest.TestCase):
    def test_convert_coordinates(self):
        observer_location = (35.6895, 139.6917, 0)
        observer_time = "2021-06-01 12:00:00"
        object_coordinates = (0, 0)
        object_location = cp.convert_coordinates(observer_location, observer_time, object_coordinates)
        self.assertEqual(object_location, (35.6895, 139.6917))

if __name__ == "__main__":
    print("test")
    unittest.main()