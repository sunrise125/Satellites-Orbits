# ---- azel_TLE.py ---------- May.25, 2019 ------
# Azimuth, elevation, distance froma TLE's elements
# Source: https://bwinkel.github.io/pycraf/satellite/index.html
import datetime
import numpy as np
from astropy.coordinates import EarthLocation
from astropy import time
from pycraf import satellite

# -------- 3-line TLE elements input ------------- 
tle_string = '''Object A (Starlink)
1 44235U 19029A   19147.55729408  .00003298  00000-0  89486-4 0  9992
2 44235  52.9978 158.9310 0001303 181.8828 296.9785 15.40623568  1490'''

satname, sat = satellite.get_sat(tle_string)
print ('\n Sat Name = ', satname,'\n Sat Epoch= ', sat.epoch, '\n Sat Numb = ', sat.satnum)

# Datetime input
year=2019; month=5; day=30; utch=2; utcm=22; utcs=0
print ('\n Prediction date: ', year, month, day )
print (' Prediction time: ', utch, utcm, utcs, ' UTC' )
dt = datetime.datetime(year, month, day, utch, utcm, utcs) 
obstime = time.Time(dt)

# define observer location
Lambda=15.0650; Phi=37.0328; H=370.0
location = EarthLocation(Lambda, Phi, H)
# create a SatelliteObserver instance
sat_obs = satellite.SatelliteObserver(location)
#
az, el, dist = sat_obs.azel_from_sat(tle_string, obstime)  
print('\n Observer location: Lambda=', Lambda, 'deg  Phi=', Phi,  'deg  H=', H,' m')
print(' Azimuth  : {:.5f}'.format(az))
print(' Elevation: {:.5f}'.format(el))
print(' Distance : {:.3f}'.format(dist))
# EOF azel_TLE.py ----------
'''
 Sat Name =  Object A (Starlink)
 Sat Epoch=  2019-05-27 13:22:30.208511
 Sat Numb =  44235

 Prediction date:  2019 5 30
 Prediction time:  2 22 0  UTC

 Observer location: Lambda= 15.065 deg  Phi= 37.0328 deg  H= 370.0  m
 Azimuth  : 137.01560 deg
 Elevation: -19.42666 deg
 Distance : 5328.111 km
'''
