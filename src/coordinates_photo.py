from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
from astropy.time import Time
import astropy.io.fits

from astroquery.astrometry_net import AstrometryNet

import os
from dotenv import load_dotenv
from pathlib import Path
from dotenv import dotenv_values
import pytz
from datetime import datetime

import math
import re


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
            
            # Parse the WCSD solution
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


def polar_to_cartesian(altitude, azimuth):
   
    radius = math.cos(altitude)*10000
    azimuth_rad = math.radians(azimuth)

    # Conversion to Cartesian coordinates
    x = radius * math.cos(azimuth_rad)
    y = radius * math.sin(azimuth_rad)
    
    return x, y

def convert_coordinates(observer_location, observer_time, object_coordinates):
    observer_lat, observer_lon, observer_alt = observer_location
    observing_location = EarthLocation(lat=observer_lat*u.deg, lon=observer_lon*u.deg)


    observer_time = Time(observer_time)
    # local_time = datetime(2024, 4, 20, 0, 37, 38)  # Example: May 13, 2024, 12:00 PM local time
    
    local_time = datetime(  observer_time.datetime.year, 
                            observer_time.datetime.month, 
                            observer_time.datetime.day, 
                            observer_time.datetime.hour, 
                            observer_time.datetime.minute, 
                            observer_time.datetime.second)  

    # Specify the local timezone
    local_timezone = pytz.timezone('Europe/Paris')  # Example: 'America/New_York'

    # Localize the time to the specified timezone
    localized_time = local_timezone.localize(local_time)

    # Convert the localized time to GMT/UTC
    gmt_time = localized_time.astimezone(pytz.utc)

    print("gmt time: ",gmt_time)

    # Celestial coordinates (RA and Dec)
    object_ra, object_dec = object_coordinates

    # Convert to Zenith and Azimuth
    celestial_coord = SkyCoord(ra=object_ra, dec=object_dec, unit="deg")
    alt_az = celestial_coord.transform_to(AltAz(obstime=gmt_time, location=observing_location))
    
    print(f"Azimuth: {alt_az.az.deg} deg\nZenith: {alt_az.alt.deg} deg")
    
    r, theta = altaz_to_polar(alt_az.alt.rad, alt_az.az.rad) 
    print(f"r: {r} m\ntheta: {theta} rad /{math.degrees(theta)} deg \n")

    # en m
    x,y = polar_to_cartesian(r,theta)
    print(f"x: {x} m, y:{y} m")

    # en rad
    object_lat, object_lon = cartesian_to_lat_lon(x,y)
    object_lat, object_lon = math.degrees(object_lat) , math.degrees(object_lon)
    # print(f"lon: {lon} rad / {math.degrees(lon)}deg, lat: {lat} rad/{math.degrees(lat)}deg")

    # object_lat, object_lon = method1(alt_az.alt.rad,alt_az.az.rad)

    return object_lat, object_lon

def altaz_to_polar(alt, az):
    """
    Submethod to convert altitude (radians), azimuth (radians) of an object to polar coordinates (r,theta)
    Returns
    -------
    (r, theta) : tuple
    """
    alt_objet = 10000 # hypothesis: flying at 10km altitude
    # r = math.cos(alt) * alt_objet # en mètres
    r = 10000 / math.cos(alt)
    theta = az
    return r, theta

def polar_to_cartesian(r, theta):
    """
    Submethod to convert radius (radians), theta (radians) of an object 
    to cartesian coordinates (x,y,z) in meters
    Returns
    -------
    (x, y, z) : tuple
    """
    x = r * math.cos(theta) # axis: 
    y = r * math.sin(theta)

    return x,y

def cartesian_to_lat_lon(x,y,z=10000):
    """
    Submethod to convert cartesian coordinates (x,y,z) in meters (we assume z = 10km) 
    to longitude, latitude in radians
    Returns
    -------
    (lon, lat) : tuple
    """
    lon = math.atan2(y,x)
    print("lon en rad:",lon)
    lat = math.atan(z/(math.sqrt(x*x+y*y)))
    return lat, lon

def get_celestial_coordinates(x,y,path):
    """
    Submethod to compute celestial coordinates of a given pixel x,y in the image 
    Returns 
    -------
    (Ra,Dec) : tuple
    """

    ref_pixel,ref_celestial,matrix = get_calibrated_camera_parameters(path)

    delta_Ra = matrix[0][0] * (x - ref_pixel[0]) + matrix[0][1] * (y - ref_pixel[1])
    delta_Dec = matrix[1][0] * (x - ref_pixel[0]) + matrix[1][1] * (y - ref_pixel[1])

    Ra = ref_celestial[0] + delta_Ra
    Dec = ref_celestial[1] + delta_Dec

    return Ra,Dec


def main():
    # Get observer location and time from file 
    with open('src/params', 'r') as file:
        lines = file.readlines()

    # Get observer location in lon, lat, alt
    observer_location = lines[0].strip()
    observer_location = tuple(map(float, observer_location.split(',')))
    
    # Get observer time using format (AAAA--MM-DD hh:mm:ss)
    observer_time = lines[1].strip()

    # Printing out the entered values
    print("Observer Location:", observer_location)
    print("Observer Time:", observer_time)

    # Get ref_pixel (x,y) in image, ref_celestial_coord (RA,Dec) of ref pixel and transformation matrix
    ref_pixel,ref_celestial_coord,matrix = get_calibrated_camera_parameters("output/wcs_header.txt")
    
    print("get celestial coord ref pixel:",get_celestial_coordinates(ref_pixel[0],ref_pixel[1],"output/wcs_header.txt"))

    # Compute Celestial coordinates of a given point (x,y) in the photo
    object_coordinates = get_celestial_coordinates(ref_pixel[0],ref_pixel[1],"output/wcs_header.txt")

    # Convert celestial coordinates of the object to GPS coordinates (longitude, latitude)
    object_location = convert_coordinates(observer_location, observer_time, object_coordinates)
    print(f"The GPS coordinates of the object are: \n{object_location}")

    return 0

if __name__ == "__main__":
    dotenv_path = Path('../.env')
    load_dotenv(dotenv_path=dotenv_path)
    main()
