from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
from astropy.time import Time
import astropy.io.fits

from astroquery.astrometry_net import AstrometryNet

import os
from dotenv import load_dotenv
from pathlib import Path
from dotenv import dotenv_values

import math


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
            
            # Parse the WCSD. Les mécanismes à l’origine d’un accident du travail
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
    observing_location = EarthLocation(lat=observer_lat*u.deg, lon=observer_lon*u.deg)


    observer_time = Time(observer_time)
    # Compute Julian Date (JD)
    jd = observer_time.jd

    # Compute Greenwich Sidereal Time (GST) in degrees
    gst = 100.46 + 0.985647 * jd
    gst = gst % 360  # Ensure GST is within 0 to 360 degrees

    # Convert GST to hours, minutes, and seconds
    gst_hour = gst / 15
    gst_min = (gst_hour % 1) * 60
    gst_sec = (gst_min % 1) * 60

    # Compute Local Sidereal Time (LST) in degrees
    lst = gst + observing_location.lon.deg
    lst = lst % 360  # Ensure LST is within 0 to 360 degrees

    lst_hour = lst / 15
    lst_min = (lst_hour % 1) * 60
    lst_sec = (lst_min % 1) * 60

    # Celestial coordinates (RA and Dec)
    object_ra, object_dec = object_coordinates
    celestial_coord = SkyCoord(ra=object_ra, dec=object_dec)

    # Convert to Zenith and Azimuth
    alt_az = celestial_coord.transform_to(AltAz(obstime=observer_time, location=observing_location))

    #call method1 with zenith (zenith = alt_az.alt) and azimut in degrees
    object_lat, object_lon = method1(alt_az.alt.deg,alt_az.az.deg)

    return object_lat, object_lon


def method1(zenith, azimuth):
    """
    Submethod to convert zenith (degrees), azimuth (degrees) of an object to longitude, latitude
    Approach : zenith azimuth -> polar coordinates (r,theta) -> cartesian coordinates -> longitude, latitude 

    Returns
    -------
    (longitude, latitude) : tuple
    """
    #assume altitude (cartesian coordinates) equals 10000
    z = 10000

    #Convert to Polar coordinates 
    r = math.cos(zenith) * z

    #coord cartesiennes
    x = r * math.sin(azimuth)
    y = r * math.cos(azimuth)
    
    lon = math.atan(y/x)
    lat = math.atan(z/(math.sqrt(x*x+y*y)))
    return lon, lat

def method2(zenith, azimuth):
    """
    Submethod to convert zenith (degrees), azimuth (degrees) of an object to longitude, latitude
    Approach : 
    https://fr.wikipedia.org/wiki/Coordonn%C3%A9es_sph%C3%A9riques#/media/Fichier:Spherical_Coordinates_(Colatitude,_Longitude)_(b).svg

    From the spherical coordinates in the image, determine zenith azimuth 
    and convert them to longitude latitude
    
    x = p * sin(zenith) * cos(azimuth)
    y = p * sin(zenith) * sin(azimuth)
    z = p * cos(zenith) = 10000 (approximation)

    Returns
    -------
    (longitude, latitude) : tuple
    """
    #assume altitude (cartesian coordinates) equals 10000
    z = 10000

    zenith = 90 -zenith

    #Cartesian coordinates
    p = z / math.cos(zenith)
    x = p * math.sin(zenith) * math.cos(azimuth)
    y = p * math.sin(zenith) * math.sin(azimuth)

    lon = math.atan(y/x)
    lat = math.atan(z/(math.sqrt(x*x+y*y)))
    return lon, lat

def main():
    #43.59332974210737, 1.4568521462595372 
    observer_location = (43.59332974210737, 1.4568521462595372, 0)
    observer_time = "2024-03-12 21:46:36"
    object_coordinates = ("4h56m33s", "3d18m15.545s")
    object_location = convert_coordinates(observer_location, observer_time, object_coordinates)
    print(f"The GPS coordinates of the object are: {object_location}")
    print(convert_coordinates(observer_location,observer_time,object_coordinates))
    # calibrate_camera('images/Image_test_12_mars_21h46mn36s.jpg')
    return 0

if __name__ == "__main__":
    dotenv_path = Path('../.env')
    load_dotenv(dotenv_path=dotenv_path)
    main()