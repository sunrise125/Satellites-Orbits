# ---- azel_TLE.py ---------- May.25, 2019 ------
# Azimuth, elevation, distance froma TLE's elements
# Source: https://bwinkel.github.io/pycraf/satellite/index.html
import datetime
import numpy as np
from astropy.coordinates import EarthLocation
from astropy import time
from pycraf import satellite

# -------- 3-line TLE elements input ------------- 
tle_string = '''Object E (Starlink)
1 44239U 19029E   19147.55806232  .00002267  00000-0  70150-4 0  9990
2 44239  52.9963 158.9546 0002970 202.5743 258.8860 15.37972393  1497'''

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
 Sat Name =  Object E (Starlink)
 Sat Epoch=  2019-05-27 13:23:36.584448
 Sat Numb =  44239

 Prediction date:  2019 5 30
 Prediction time:  2 22 0  UTC

 Observer location: Lambda= 15.065 deg  Phi= 37.0328 deg  H= 370.0  m
 Azimuth  : -161.99410 deg
 Elevation: 48.90299 deg
 Distance : 587.296 km
'''
