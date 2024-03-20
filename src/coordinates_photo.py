from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
from astropy.time import Time


def calibrate_camera(photo, observer_location):
    """
    This function calibrates the camera and returns the camera parameters

    Parameters
    ----------
    photo : str (path to the photo)
        The photo to calibrate the camera
    observer_location : tuple
        The observer location on Earth in (latitude, longitude, altitude) format

    Returns
    -------
    reference_pixel : tuple
        The reference pixel of the camera in (x, y) format
    reference_point: tuple
        The reference point of the camera in (RA, dec) format
    transformation_matrix : matrix
        The transformation matrix of the camera
    """
    photo = None
    observer_location = None

    reference_pixel = None
    reference_point = None
    transformation_matrix = None
    return reference_pixel, reference_point, transformation_matrix




def convert_coordinates(observer_location, observer_time, object_coordinates):
    """
    This function converts the celestial coordinates of an object to Earth-based coordinates and then to GPS coordinates
    We assume that the object flying at 32000 feet altitude

    Parameters
    ----------
    observer_location : tuple
        The observer location on Earth in (latitude, longitude, altitude) format
    observer_time : str
        The time of observation in the format "YYYY-MM-DD HH:MM:SS"
    object_coordinates : tuple
        The celestial coordinates of the object in (RA, Dec) format

    Returns
    -------
    object_location : tuple
        The GPS coordinates of the object in (latitude, longitude) format
    """
    observer_lat, observer_lon, observer_alt = observer_location
    observer_time = Time(observer_time)
    object_ra, object_dec = object_coordinates

    object_lat = observer_lat
    object_lon = observer_lon
    return object_lat, object_lon

def main():
    observer_location = (35.6895, 51.3896, 0)
    observer_time = "2021-06-21 00:00:00"
    object_coordinates = (0, 0)
    object_location = convert_coordinates(observer_location, observer_time, object_coordinates)
    print(f"The GPS coordinates of the object are: {object_location}")
    return 0