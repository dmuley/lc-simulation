if __name__ == '__main__':
    if __package__ is None:
		import sys
		from os import path
		sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) );
		print sys.path[-1];
		from gravsims.basicKepler import *;
		from integrators.polarIntegrator3 import *;
    else:
		from ..gravsims.basicKepler import *;
                from integrators.polarIntegrator3 import *;
	
import numpy as np;
from itertools import chain;
import matplotlib.pyplot as plt;
import general_transits as gt;
import time;

#OPERATING EXAMPLE TO RUN THIS CODE
#FOR RETROGRADE ORBITS: Add pi radians to inclination, add pi radians to arg_periastron

steps = 7501;
revs = 0.125;

starPlanet = OrbitingSystem();
starPlanet.bodies = [0,0]

#Masses of stellar system
q = OrbitingSystem();
q.bodies = [0,0];

q.bodies[0] = OrbitingSystem()
q.bodies[0].mass = 332946. * 0.6897;
q.bodies[0].temperature = 0.985;
q.bodies[0].radius = 0.6489 * 0.00464913034;

q.bodies[1] = OrbitingSystem();
q.bodies[1].mass = 332946. * 0.20255;
q.bodies[1].semimajor = 1. * 0.22431;
q.bodies[1].temperature = 0.015;
q.bodies[1].inclination = np.pi * 90./180.
q.bodies[1].radius = 0.22623 * 0.00464913034;
q.bodies[1].phase = -np.pi/8.;

q.setTotalMass();

#Masses of planetary system
r = OrbitingSystem();
r.semimajor = 0.7048;
r.eccentricity = 0.0069;
r.inclination = np.pi * 90./180.;
r.phase = np.pi * -30./180.

r.bodies = [0, 0];
r.bodies[0] = OrbitingSystem();
r.bodies[0].mass = 317.83 * 0.333;
r.bodies[0].radius = 0.000477894503 * 0.7538;

#Nominal moon
r.bodies[1] = OrbitingSystem();
r.bodies[1].mass = 0.00001;
r.bodies[1].semimajor = 1.;
r.bodies[1].radius = 0.0000000001;
r.bodies[1].inclination = np.pi * 0.5;

r.setTotalMass();

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
#Sorting by Z-position	
fx, fy, fz, zpos = gt.sort_keys(final_x, final_y, final_z);
#getting points of intersection to check for transit
lb, intersects = gt.arrange_combinations(fx, fy, transit_array);

l = time.time()
n = 31;
ta = transit_array;
cadence, light_blocked = gt.generate_lightcurve(fx, fy, lb, intersects, n, ta, times, zpos);

print time.time() - l;

#Plotting final light curve.
plt.clf();
plt.plot(cadence, -light_blocked, '-');
plt.xlim(cadence[0], cadence[len(times) - 1]);
plt.show();
