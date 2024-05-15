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
import pyproj


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


def convert_to_gps_coordinates(observer_location, observer_time, object_coordinates):
    """
    This function converts celestial coordinates to GPS coordinates.

    First we convert the observer's location to cartesian coordinates.
    Then we convert the altitude and azimuth of the object to cartesian coordinates.
    Finally, we convert the cartesian coordinates of the object to longitude and latitude.

    Parameters
    ----------
    observer_location : tuple
        The observer location in (latitude, longitude, altitude) format
    observer_time : str
        The observer time in "YYYY-MM-DD HH:MM:SS" format
    object_coordinates : tuple
        The celestial coordinates of the object in (RA, Dec) format

    Returns
    -------
    object_location : tuple
        The GPS coordinates of the object in (longitude, latitude) format
    """
    # Unpack observer data
    observer_lon, observer_lat, observer_alt = observer_location
    observing_location = EarthLocation(lat=observer_lat*u.deg, lon=observer_lon*u.deg)

    # Specify the observer time to be in the UTC timezone
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

    print(f"Observer Time (GMT): {observer_time}")

    # Unpack object data
    object_ra, object_dec = object_coordinates
    celestial_coord = SkyCoord(ra=object_ra, dec=object_dec, unit="deg", frame='icrs')

    # Convert to Zenith and Azimuth observing_location
    alt_az = celestial_coord.transform_to(AltAz(obstime=gmt_time, location=observing_location))

    # Get Altitude and Azimuth of the object
    E = alt_az.alt.rad # in radians (A) (Altitude)
    A = alt_az.az.rad # in radians (E) (Azimuth)
    
    print('\n')
    print(f"{'Obj Property':<20} {'Degrees':<15}")
    print(f"{'-'*35}")
    print(f"{'Object RA':<20} {object_ra:<15.4f}")
    print(f"{'Object Dec':<20} {object_dec:<15.4f}")
    print()
    print(f"{'Airplane Altitude':<20} {alt_az.alt.deg:<15.4f}")
    print(f"{'Airplane Azimuth':<20} {alt_az.az.deg:<15.4f}")
    print(f"{'-'*35}")

    # Convert longitude, latitude of the observer to cartesian coordinates
    # To convert longitude, latitude to cartesian coordinates, we use the Mercator projection as an approximation
    # To compute the x, y coordinates manually, we use the following formulas:
    # xc = R * math.cos(obs_lat) * math.cos(obs_lon) # in m
    # yc = R * math.cos(obs_lat) * math.sin(obs_lon) # in m
    # Here, we use the pyproj library to convert longitude, latitude to x, y coordinates

    obs_lat = math.radians(observer_lat)
    obs_lon = math.radians(observer_lon)
   
    transformer_mercator = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:3857")
    yc, xc = transformer_mercator.transform(observer_lon, observer_lat)

    # Convert radians to degrees
    obs_lon_deg = math.degrees(obs_lon)
    obs_lat_deg = math.degrees(obs_lat)

    print('\n')
    print(f"{'Observer Property':<20} {'Value':<19} {'Value2':<15}")
    print(f"{'-'*58}")
    print(f"{'Observer X':<20} {xc/1000.0:<15.4f} km")
    print(f"{'Observer Y':<20} {yc/1000.0:<15.4f} km")
    print(f"{'Observer Longitude':<20} {obs_lon:<15.4f} rad {obs_lon_deg:<10.4f} deg")
    print(f"{'Observer Latitude':<20} {obs_lat:<15.4f} rad {obs_lat_deg:<10.4f} deg")
    print(f"{'-'*58}")

    # Hypothesis: flying at 10km altitude
    alt_objet = 10000 # in m
    R = 6371000 # Earth radius in m

    # We make a triangle to convert altitude, azimuth to coordinates of a horizontal plane x-y 
    # tangent to the sphere of the Earth at the observer's location
    
    Rprime = R + observer_alt
    alpha = E + math.pi/2
    beta = math.asin((Rprime * math.sin(A)) / (Rprime + alt_objet))
    psi =  alpha + beta
    d = psi * R
    x = xc + d * math.cos(A)
    y = yc + d * math.sin(A)

    # Convert radians to degrees
    alpha_deg = math.degrees(alpha)
    beta_deg = math.degrees(beta)
    psi_deg = math.degrees(psi)
    E_deg = math.degrees(E)
    A_deg = math.degrees(A)

    # Print as a table with fixed width columns
    print("\n")
    print(f"{'Obj Property':<12} {'Radians':<12} {'Degrees':<15}")
    print(f"{'-'*40}")
    print(f"{'Alpha':<12} {alpha:<12.4f} {alpha_deg:<15.4f}")
    print(f"{'Beta':<12} {beta:<12.4f} {beta_deg:<15.4f}")
    print(f"{'Psi':<12} {psi:<12.4f} {psi_deg:<15.4f}")
    print(f"{'Altitude':<12} {E:<12.4f} {E_deg:<15.4f}")
    print(f"{'Azimuth':<12} {A:<12.4f} {A_deg:<15.4f}")
    print()
    print(f"{'Distance':<12} {d/1000.0:<12.4f} {'km':<15}")
    print(f"{'X':<12} {x/1000.0:<12.4f} {'km':<15}")
    print(f"{'Y':<12} {y/1000.0:<12.4f} {'km':<15}")
    print(f"{'-'*40}\n")

    # Convert x, y to longitude, latitude
    # lon = math.atan2(y,x)
    # lat = math.atan2(z,(math.sqrt(x*x+y*y)))
    transformer_wgs64 = pyproj.Transformer.from_crs("EPSG:3857", "EPSG:4326")
    lon, lat = transformer_wgs64.transform(y, x)

    return lon, lat

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


    wcs_header_path = "output/wcs_header_20avril.txt"

    # Get ref_pixel (x,y) in image, ref_celestial_coord (RA,Dec) of ref pixel and transformation matrix
    ref_pixel,ref_celestial_coord,matrix = get_calibrated_camera_parameters(wcs_header_path)
    
    print(get_celestial_coordinates(ref_pixel[0],ref_pixel[1],wcs_header_path))
    
    # Compute Celestial coordinates of a given point (x,y) in the photo
    object_coordinates = get_celestial_coordinates(ref_pixel[0],ref_pixel[1],wcs_header_path)
    
    # Convert celestial coordinates of the object to GPS coordinates (longitude, latitude)
    # object_location = convert_coordinates(observer_location, observer_time, object_coordinates)
    # object_coordinates = (743,543)
    object_coordinates = get_celestial_coordinates(0,0,wcs_header_path)
    object_location = convert_to_gps_coordinates(observer_location, observer_time, object_coordinates)
    print(f"The GPS coordinates of the object are: \n{object_location}")

    image_width = 1600
    image_height = 1066
    image_left_top = (0,0)
    image_right_top = (image_width,0)
    image_right_bottom = (image_width,image_height)
    image_left_bottom = (0,image_height)

    left_top = get_celestial_coordinates(image_left_top[0],image_left_top[1],wcs_header_path)
    right_top = get_celestial_coordinates(image_right_top[0],image_right_top[1],wcs_header_path)
    right_bottom = get_celestial_coordinates(image_right_bottom[0],image_right_bottom[1],wcs_header_path)
    left_bottom = get_celestial_coordinates(image_left_bottom[0],image_left_bottom[1],wcs_header_path)

    left_top_gps = convert_to_gps_coordinates(observer_location, observer_time, left_top)
    right_top_gps = convert_to_gps_coordinates(observer_location, observer_time, right_top) 
    right_bottom_gps = convert_to_gps_coordinates(observer_location, observer_time, right_bottom)
    left_bottom_gps = convert_to_gps_coordinates(observer_location, observer_time, left_bottom)

    print("\n")
    print(f"Left Top GPS: {left_top_gps}")
    print(f"Right Top GPS: {right_top_gps}")
    print(f"Right Bottom GPS: {right_bottom_gps}")
    print(f"Left Bottom GPS: {left_bottom_gps}")

    return 0

if __name__ == "__main__":
    dotenv_path = Path('../.env')
    load_dotenv(dotenv_path=dotenv_path)
    main()
