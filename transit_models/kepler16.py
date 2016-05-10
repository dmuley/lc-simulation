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

#OPERATING EXAMPLE TO RUN THIS CODE
#FOR RETROGRADE ORBITS: Add pi radians to inclination, add pi radians to arg_periastron

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
#r.mass = 317.83 * 0.333;
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
r.bodies[1].radius = 0.0000000001
r.bodies[1].inclination = np.pi * 0.5

r.setTotalMass();

#Normalizing times of each body
starPlanet.bodies[0], starPlanet.bodies[1] = q, r;
starPlanet.setBaseTime();
year = starPlanet.bt
starPlanet.bodies[0].bt, starPlanet.bodies[1].bt = year, year;
starPlanet.setSystemOrbits();
times = starPlanet.times;

starPlanet.bodies[0].setSystemOrbits();
starPlanet.bodies[1].setSystemOrbits();

	
print " "
print len(starPlanet.bodies[0].bodies);
print " "

final_x = [];
final_y = [];
final_z = [];
final_array = []
transit_array = [];

#print starPlanet.bodies;

#Flattening list of bodies (probably want a more general interpretation later for jagged array of bodies).

for a in range(0, len(starPlanet.bodies)):
	for b in range(0, len(starPlanet.bodies[a].bodies)):
		print (a, b)
		x = starPlanet.bodies[a].orbits[b][0] + starPlanet.orbits[a][0] - starPlanet.bodies[a].orbits[0][0] * 0 - starPlanet.orbits[0][0] * 0;
		y = starPlanet.bodies[a].orbits[b][1] + starPlanet.orbits[a][1] - starPlanet.bodies[a].orbits[0][1] * 0 - starPlanet.orbits[0][1] * 0;
		z = starPlanet.bodies[a].orbits[b][2] + starPlanet.orbits[a][2] - starPlanet.bodies[a].orbits[0][2] * 0 - starPlanet.orbits[0][2] * 0;
		
		final_array.append(np.array([x, y, z]));
		final_x.append(x);
		final_y.append(y);
		final_z.append(z);

		starPlanet.bodies[a].orbits[b][0] = x;
		starPlanet.bodies[a].orbits[b][1] = y;
		starPlanet.bodies[a].orbits[b][2] = z;
		
		transit_array.append(starPlanet.bodies[a].bodies[b]);
		
		plt.plot(times, y/np.absolute(y) * np.sqrt(x**2 + y**2), '-');
		#plt.plot(x, y);
plt.ylabel("Difference from mean (AU)");
plt.xlabel("Time (days)");
plt.title("Sky-projected deviation of bodies in planetary system from COM");
plt.xlim(times[0], times[len(times) - 1]);
plt.show();		

fx, fy, fz, zpos = gt.sort_keys(final_x, final_y, final_z);
lb, intersects = gt.arrange_combinations(fx, fy, transit_array);

n = 31;
ta = transit_array;
cadence, light_blocked = gt.generate_lightcurve(fx, fy, lb, intersects, n, ta, times, zpos);

#Plotting final light curve.
plt.clf();
plt.plot(cadence, -light_blocked, '-');
plt.xlim(cadence[0], cadence[len(times) - 1]);
plt.show();
