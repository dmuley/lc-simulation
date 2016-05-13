if __package__ is None:
	import sys;
	from os import path;
	sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) );
	print sys.path[-1];
	from gravsims.basicKepler import *;
	from processing_pipeline.pipeline import import_idl;	
else:
	from ..gravsims.basicKepler import *;
	from ..processing_pipeline.pipeline import import_idl;

import numpy as np;
from itertools import chain;
import matplotlib.pyplot as plt;
import general_transits as gt;
import time;
from scipy.stats import chisquare;

################################################
######## PLACE ALL PLANETARY DATA HERE  ########
################################################
#FOR RETROGRADE ORBITS: Add pi radians to inclination, add pi radians to arg_periastron

steps = 15001;
revs = 0.25;

starPlanet = OrbitingSystem();
starPlanet.bodies = [0,0]

#Masses of stellar system
q = OrbitingSystem();
#q.arg_periastron = np.pi;
q.bodies = [0,0];

q.bodies[0] = OrbitingSystem()
q.bodies[0].mass = 332946. * 0.6897;
q.bodies[0].temperature = 0.985;
q.bodies[0].radius = 0.6489 * 0.00464913034;
q.bodies[0].ld_coeffs = [2.2,0,1.0,0];
q.bodies[0].ld_powers = [0.0, 0.5, 1.0, 1.5];

q.bodies[1] = OrbitingSystem();
q.bodies[1].mass = 332946. * 0.20255;
q.bodies[1].semimajor = 1. * 0.22431;
q.bodies[1].temperature = 0.015;
q.bodies[1].inclination = np.pi * 90./180.
q.bodies[1].radius = 0.22623 * 0.00464913034;
q.bodies[1].phase = -np.pi/8 + np.pi/30.;
q.bodies[1].eccentricity = 0.15944;
q.bodies[1].arg_periastron = -np.pi/30;

q.setTotalMass();

#Masses of planetary system
r = OrbitingSystem();
r.semimajor = 0.7048;
r.eccentricity = 0.0069;
r.inclination = np.pi * 90./180.;
r.phase = np.pi * -20./180.

r.bodies = [0, 0];
r.bodies[0] = OrbitingSystem();
r.bodies[0].mass = 317.83 * 0.333;
r.bodies[0].radius = 0.000477894503 * 0.7538;

#Nominal moon
r.bodies[1] = OrbitingSystem();
r.bodies[1].mass = 1.;
r.bodies[1].semimajor = 0.0075;
r.bodies[1].radius = 0.000477894503 * 0.0892141778;
#r.bodies[1].inclination = np.pi * 0.5;
r.bodies[1].phase = -np.pi/3.;

r.setTotalMass();

################################################
########### END PLANETARY DATA HERE  ###########
################################################

#Normalizing times of each body
starPlanet.bodies[0], starPlanet.bodies[1] = q, r;
starPlanet.setBaseTime();
year = starPlanet.bt
starPlanet.bodies[0].bt, starPlanet.bodies[1].bt = year, year;
starPlanet.setSystemOrbits(s = steps, r = revs);
times = starPlanet.times;

starPlanet.bodies[0].setSystemOrbits(s = steps, r = revs);
starPlanet.bodies[1].setSystemOrbits(s = steps, r = revs);

#Getting positions relative to center of mass
final_x, final_y, final_z, transit_array = gt.traverse_tree(starPlanet, times);	
gt.plot_distances(final_x, final_y, times);

#Sorting by Z-position	
fx, fy, fz, zpos = gt.sort_keys(final_x, final_y, final_z);
fx2,fy2,fz2,zpos2= gt.sort_keys(final_x[:-1], final_y[:-1], final_z[:-1]);

#getting points of intersection to check for transit
lb, intersects = gt.arrange_combinations(fx, fy, transit_array);
lb2, intersects2 = gt.arrange_combinations(fx2, fy2, transit_array[:-1]);

l = time.time()
n = 31;
ta = transit_array;
cadence, light_blocked = gt.generate_lightcurve(fx, fy, lb, intersects, n, ta, times, zpos);
cadence2, light_blocked_2 = gt.generate_lightcurve(fx2, fy2, lb2, intersects2,n, ta[:-1], times, zpos2);

print "Time to compute transits ",
print time.time() - l;

plt.clf();

plt.plot(cadence, light_blocked-light_blocked_2,   'b')
#plt.plot(cadence, -light_blocked_2, 'r');
plt.xlim(cadence[0], cadence[-1])

plt.xlabel("Time (days)")
plt.ylabel("Fraction of light blocked")
plt.title("Predicted transits of Kepler-16 versus actual detrended data")

plt.show();

chisquare = (light_blocked_2 - light_blocked)**2
print "Chi-squared: " + str(np.sum(chisquare));
