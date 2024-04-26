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
    celestial_coord = SkyCoord(ra=object_ra, dec=object_dec, unit="deg")

    # Convert to Zenith and Azimuth
    alt_az = celestial_coord.transform_to(AltAz(obstime=observer_time, location=observing_location))

    print(f"Azimuth: {alt_az.az.deg} deg\nZenith: {alt_az.alt.deg} deg")
    
    r, theta = altaz_to_polar(alt_az.alt.rad, alt_az.az.rad) 
    print(f"r: {r} m\ntheta: {theta} rad /{math.degrees(theta)} deg \n")

    # en m
    x,y = polar_to_cartesian(r,theta)
    print(f"x: {x} m, y:{y} m")

    # en rad
    object_lat, object_lon = cartesian_to_lon_lat(x,y)
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
    r = math.cos(alt) * alt_objet # en mètres
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
    x = r * math.sin(theta) # axis: 
    y = r * math.cos(theta)

    return x,y

def cartesian_to_lon_lat(x,y):
    """
    Submethod to convert cartesian coordinates (x,y,z) (we assume z = 10km) 
    to longitude, latitude in radians
    Returns
    -------
    (lon, lat) : tuple
    """
    z = 10000 # en m
    lon = math.atan(y/x)
    lat = math.atan(z/(math.sqrt(x*x+y*y)))
    return lon, lat

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
    
    # Prompting the user for observer location
    observer_location_input = input("Enter observer location (latitude, longitude, altitude): ")
    # Splitting the input string into latitude, longitude, and altitude
    observer_location = tuple(map(float, observer_location_input.split(',')))

    # Prompting the user for observer time
    observer_time = input("Enter observer time (YYYY-MM-DD HH:MM:SS): ")

    # Printing out the entered values
    print("Observer Location:", observer_location)
    print("Observer Time:", observer_time)

    # Get ref_pixel (x,y) in image, ref_celestial_coord (RA,Dec) of ref pixel and transformation matrix
    ref_pixel,ref_celestial_coord,matrix = get_calibrated_camera_parameters("../output/wcs_header.txt")
    
    print(get_celestial_coordinates(ref_pixel[0],ref_pixel[1],"../output/wcs_header.txt"))
    
    # Compute Celestial coordinates of a given point (x,y) in the photo
    object_coordinates = get_celestial_coordinates(ref_pixel[0],ref_pixel[1],"../output/wcs_header.txt")
    
    # Convert celestial coordinates of the object to GPS coordinates (longitude, latitude)
    object_location = convert_coordinates(observer_location, observer_time, object_coordinates)
    print(f"The GPS coordinates of the object are: \n{object_location}")

    return 0

if __name__ == "__main__":
    dotenv_path = Path('../.env')
    load_dotenv(dotenv_path=dotenv_path)
    main()