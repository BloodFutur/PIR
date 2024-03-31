from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
from astropy.time import Time
import astropy.io.fits

from astroquery.astrometry_net import AstrometryNet

import os
from dotenv import load_dotenv
from pathlib import Path
from dotenv import dotenv_values


def calibrate_camera(photo: str):
    """
    This function calibrates the camera and returns the reference pixel, reference point, and transformation matrix
    Some other camera parameters are stored in /output/wcs_header.txt

    If the calibration fails, the function returns None for all the parameters

    Parameters
    ----------
    photo : str (path to the photo)
        The photo to calibrate the camera
    Returns
    -------
    reference_pixel : tuple
        The reference pixel of the camera in (x, y) format
    reference_point: tuple
        The reference point of the camera in (RA, dec) format
    transformation_matrix : matrix
        The transformation matrix of the camera
    """
    ast = AstrometryNet()

    # Load the API key from the .env file
    config = dotenv_values(".env")  
    API_KEY_ASTROMETRY = config['API_KEY_ASTROMETRY']
    # print(API_KEY_ASTROMETRY)
    ast.api_key = API_KEY_ASTROMETRY

    # Set up variables processing the image
    try_again = True
    submission_id = None

    # While processing the image, keep trying until the solve is successful
    while try_again:
        try:
            # Send the image to Astrometry.net to solve
            if not submission_id:
                wcs_header = ast.solve_from_image(photo, submission_id=submission_id)
            # Monitor the submission
            else:
                wcs_header = ast.monitor_submission(submission_id, solve_timeout=120)
        except TimeoutError as e: # handle the timeout error
            submission_id = e.args[1]
        else: # got a result, so terminate
            print('Solve succeeded')
            try_again = False

    # If the solve is successful
    if wcs_header:
        print('Writing WCS solution to file')

        # Save the WCS solution to a file
        with open('output/wcs_header.txt', 'w') as f:
            f.write(wcs_header.tostring())

        # Get the reference pixel
        reference_pixel = (wcs_header['CRPIX1'], wcs_header['CRPIX2'])

        # Get the reference point
        reference_point = (wcs_header['CRVAL1'], wcs_header['CRVAL2'])

        # Get the transformation matrix
        transformation_matrix = ([wcs_header['CD1_1'], wcs_header['CD1_2']],
                                 [wcs_header['CD2_1'], wcs_header['CD2_2']])
        
    else:
        # Code to execute when solve fails
        print('Solve failed')
        reference_pixel = None
        reference_point = None
        transformation_matrix = None

    return reference_pixel, reference_point, transformation_matrix


def get_calibrated_camera_parameters(path: str='output/wcs_header.txt'):
    """
    This function gets the calibrated camera parameters from the cache.
    The calibrated camera parameters are stored in /output/wcs_header.txt
    If the file is not found, the function returns None for all the parameters

    Parameters
    ----------
    None

    Returns
    -------
    reference_pixel : tuple
        The reference pixel of the camera in (x, y) format
    reference_point: tuple
        The reference point of the camera in (RA, dec) format
    transformation_matrix : matrix
        The transformation matrix of the camera
    """
    # Read the WCS solution from the file
    try:
        with open(path, 'r') as f:
            # Read the WCS solution from the file
            wcs_header_raw = f.read()
            
            # Parse the WCS solution
            wcs_header = astropy.io.fits.Header.fromstring(wcs_header_raw)

            # Get the reference pixel
            reference_pixel = (wcs_header['CRPIX1'], wcs_header['CRPIX2'])

            # Get the reference point
            reference_point = (wcs_header['CRVAL1'], wcs_header['CRVAL2'])

            # Get the transformation matrix
            transformation_matrix = ([wcs_header['CD1_1'], wcs_header['CD1_2']],
                                    [wcs_header['CD2_1'], wcs_header['CD2_2']])

    except FileNotFoundError:
        # Code to execute when the file is not found
        print('File not found')
        reference_pixel = None
        reference_point = None
        transformation_matrix = None

    return reference_pixel, reference_point, transformation_matrix



def local_sidereal_time(observer_lon, observer_time):
    #la date
    dt = datetime.strptime(observer_time, "%Y-%m-%d %H:%M:%S")
    #la date julienne 
    jd = 2451545.0 + (dt - datetime(2000, 1, 1)).total_seconds() / 86400.0
    #temps en siècles juliens 
    t = (jd - 2451545.0) / 36525.0
    #temps Moyen de Greenwich Sidéral
    gmst = 280.46061837 + 360.98564736629 * (jd - 2451545) + 0.000387933 * t**2 - t**3 / 38710000
    #Ajustement au Temps Sidéral Local
    lst = (gmst + observer_lon) % 360 # Le modulo 360 assure que le résultat reste dans l'intervalle de 0° à 360°.
    return lst

def equatorial_to_azimuthal(ra, dec, lst, observer_lat):
    ra_rad = math.radians(ra)
    dec_rad = math.radians(dec)
    lst_rad = math.radians(lst)
    observer_lat_rad = math.radians(observer_lat)

    sin_altitude = math.sin(dec_rad) * math.sin(observer_lat_rad) + \
                   math.cos(dec_rad) * math.cos(observer_lat_rad) * math.cos(hour_angle)
    altitude = math.asin(sin_altitude)
    
    cos_azimuth = (math.sin(dec_rad) - math.sin(altitude) * math.sin(observer_lat_rad)) / \
                  (math.cos(altitude) * math.cos(observer_lat_rad))
    azimuth = math.acos(cos_azimuth)
   
    return math.degrees(azimuth), math.degrees(altitude)

def polar_to_cartesian(altitude, azimuth):
   
    radius = math.cos(altitude)*10000
    azimuth_rad = math.radians(azimuth)

    # Conversion to Cartesian coordinates
    x = radius * math.cos(azimuth_rad)
    y = radius * math.sin(azimuth_rad)
    
    return x, y

def convert_coordinates(observer_location, observer_time, object_coordinates):
    observer_lat, observer_lon, observer_alt = observer_location
    
    observer_time_dt = datetime.strptime(observer_time, "%Y-%m-%d %H:%M:%S")
    
    object_ra, object_dec = object_coordinates

    lst = local_sidereal_time(observer_lon, observer_time_dt)
    
    azimuth, altitude = equatorial_to_azimuthal(object_ra, object_dec, lst, observer_lat)


    x, y = polar_to_cartesian(radius, azimuth)

    # Convert Cartesian coordinates to GPS coordinates 
    λ = math.atan2(y, x)
    φ = math.atan2(z, math.sqrt(x**2 + y**2))

    object_lon = math.degrees(λ)
    object_lat = math.degrees(φ)

    return object_lat, object_lon



def main():

    # observer_location = (35.6895, 51.3896, 0)
    # observer_time = "2021-06-21 00:00:00"
    # object_coordinates = (0, 0)
    # object_location = convert_coordinates(observer_location, observer_time, object_coordinates)
    # print(f"The GPS coordinates of the object are: {object_location}")

    # calibrate_camera('images/Image_test_12_mars_21h46mn36s.jpg')
    return 0

if __name__ == "__main__":
    dotenv_path = Path('../.env')
    load_dotenv(dotenv_path=dotenv_path)
    main()
