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
import numpy as np
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


def convert_coordinates(observer_location, observer_time, object_coordinates):
    # Define observer location in the Astropy format
    observer_lat, observer_lon, observer_alt = observer_location
    observing_location = EarthLocation(lat=observer_lat*u.deg, lon=observer_lon*u.deg)

    # Create time object
    observer_time = Time(observer_time)
    
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
    print("observing location:",observer_location)
    alt_az = celestial_coord.transform_to(AltAz(obstime=gmt_time, location=observing_location))

    print(f"Azimuth: {alt_az.az.deg} deg\nZenith: {alt_az.alt.deg} deg")

    # Calculate latitude and longitude
    latitude, longitude = calculate_lon_lat(observer_lat, observer_lon,
                                                    alt_az.alt.deg, alt_az.az.deg)

    return latitude, longitude

def calculate_lon_lat(observer_lat, observer_lon, elevation, azimuth,aircraft_alt = 10):
    """
    assume aircrafts fly object_alt = 10km
    """
    # Earth's radius in kilometers
    Earth_radius = 6371.1363
    
    # Convert angles from degrees to radians
    observer_lat_rad = np.radians(observer_lat)
    observer_lon_rad = np.radians(observer_lon)
    elevation_rad = np.radians(elevation)
    azimuth_rad = np.radians(azimuth)

    # Calculate intermediate angle
    psi_alt = np.pi/2 - elevation_rad - np.arcsin(Earth_radius/(Earth_radius+aircraft_alt)*np.cos(elevation_rad))
    
    #calculate lat
    lat = np.arcsin(np.sin(observer_lat_rad) * np.cos(psi_alt) + 
                        np.cos(observer_lat_rad) * np.sin(psi_alt) * np.cos(azimuth_rad))
    
    # Calclute longitude based on latitude's val
    if (observer_lat>70   and  np.tan(psi_alt)*np.cos(azimuth_rad) > np.tan(np.pi/2-observer_lat_rad)) or \
       (observer_lat<-70  and  np.tan(psi_alt)*np.cos(azimuth_rad + np.pi) > np.tan(np.pi/2+observer_lat_rad)):
      lon = observer_lon_rad +np.pi - np.arcsin(np.sin(psi_alt)*np.sin(azimuth_rad)/np.cos(lat))
      
    else:
      lon = observer_lon_rad + np.arcsin(np.sin(psi_alt)*np.sin(azimuth_rad)/np.cos(lat))
    
    return np.degrees(lat), np.degrees(lon)


def get_celestial_coordinates(x,y,path):
    """
    Submethod to compute celestial coordinates of a given pixel x,y in the image using data calibration (Astrometry)
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

    # Define the expected format
    expected_format = r"-?\d+(\.\d+)?,\s*-?\d+(\.\d+)?,\s*-?\d+(\.\d+)?"

    # Prompting the user for observer location
    observer_location_input = input("Enter observer location (latitude, longitude, altitude): ")

    # Check if input matches the expected format
    if re.match(expected_format, observer_location_input):
        print("Input is in the correct format.")
    else:
        print("Input is not in the correct format. Please enter in latitude, longitude, altitude format.")
        return -1

    # Splitting the input string into latitude, longitude, and altitude
    observer_location = tuple(map(float, observer_location_input.split(',')))


    # Define the expected format
    expected_format = r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}"

    # Prompting the user for observer time
    observer_time = input("Enter observer time (YYYY-MM-DD HH:MM:SS): ")

    # Check if input matches the expected format
    if re.match(expected_format, observer_time):
        print("Input is in the correct format.")
    else:
        print("Input is not in the correct format. Please enter in YYYY-MM-DD HH:MM:SS format.")
        return -1

    # Printing out the entered values
    print("Observer Location:", observer_location)
    print("Observer Time:", observer_time)

    # Get ref_pixel (x,y) in image, ref_celestial_coord (RA,Dec) of ref pixel and transformation matrix
    print("get celestial coord ref pixel:")


    # Compute Celestial coordinates of a given point (x,y) in the photo
    object_coordinates = get_celestial_coordinates(0,0,"output/wcs_header.txt")

    # Convert celestial coordinates of the object to GPS coordinates (latitude, longitude)
    object_location = convert_coordinates(observer_location, observer_time, object_coordinates)
    print(f"The GPS coordinates of the object are: \n{object_location}")

    return 0

if __name__ == "__main__":
    dotenv_path = Path('../.env')
    load_dotenv(dotenv_path=dotenv_path)
    main()
