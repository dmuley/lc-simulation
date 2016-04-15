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

#OPERATING EXAMPLE TO RUN THIS CODE
starPlanet = OrbitingSystem();
starPlanet.bodies = [0,0]

q = OrbitingSystem();
q.bodies = [0,0];

q.bodies[0] = OrbitingSystem()
q.bodies[0].mass = 332946. * 0.6897;
q.bodies[0].temperature = 0.985;
q.bodies[0].radius = 0.6489 * 0.00464913034;

q.bodies[1] = OrbitingSystem();
q.bodies[1].mass = 332946. * 0.20255;
q.bodies[1].semimajor = 0.22431;
q.bodies[1].temperature = 0.015;
q.bodies[1].inclination = np.pi * 90./180.
q.bodies[1].radius = 0.22623 * 0.00464913034;

q.setTotalMass();

r = OrbitingSystem();
#r.mass = 317.83 * 0.333;
r.semimajor = 0.7048;
r.eccentricity = 0.0069;
r.inclination = np.pi * 90./180.;

r.bodies = [0, 0];
r.bodies[0] = OrbitingSystem();
r.bodies[0].mass = 317.83 * 0.333;
r.bodies[0].radius = 0.000477894503 * 0.7538;

r.bodies[1] = OrbitingSystem();
r.bodies[1].mass = 0.00001;
r.bodies[1].semimajor = 0.000001;
r.bodies[1].radius = 0.00001

r.setTotalMass();

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

print starPlanet.bodies;
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
		
		plt.plot(times, y/np.absolute(y) * np.sqrt(x**2 + y**2));
		#plt.plot(x, y);
plt.ylabel("Difference from mean (AU)");
plt.xlabel("Time (days)");
plt.title("Sky-projected deviation of bodies in planetary system from COM");
plt.xlim(times[0], times[len(times) - 1]);

#plt.xlim(np.array(final_x).min(), np.array(final_x).max());
#plt.xlim(np.array(final_y).min(), np.array(final_y).max());
#[plt.plot(final_x[a], final_y[a]) for a in range(0,len(final_x))];

plt.show();		

#plot_orbit(a=1.5, args = final_array);

###### LOOM ALGORITHM ######
# Using the Z-values we will index the orders of X and Y values, and then get transits for each body.

fx = np.array(final_x).T;
fy = np.array(final_y).T;
fz = np.array(final_z).T;

zpos = np.argsort(fz);
#fx2 = fx[np.arange(np.shape(fx)[0])[:,np.newaxis], np.argsort(fz)];
#fy2 = fy[np.arange(np.shape(fy)[0])[:,np.newaxis], np.argsort(fz)];

print len(fx), len(times);

#### ACTUAL TRANSIT MODELING BELOW ####
light_blocked = np.array([]);

n = 31;
ta = transit_array;

for m in range(0,len(zpos)):
	c = np.array([[ta[r].radius, fx[m][r], fy[m][r]] for r in zpos[m]]);
	L_perarea = np.array([ta[r].temperature for r in zpos[m]])
	
	if False in (np.array(zpos[m]) == np.arange(0,len(zpos[m]))):
                print zpos[m];
		print L_perarea;

	s = 0.;

	for ind in range(0,len(ta)):
		if L_perarea[ind] != 0:
			c_trunc = c[ind:len(ta)];
			ti = getTangentsIntersections(c_trunc);
			rt = generateRadiiThetas(n, ti[0], ti[1]);

			f = rd2(rt, c_trunc, opt = 0);
			sr = c_trunc[0][0];

			oh = groupAndIntegrate(bounds = f, num = n, star_rad = sr, ld_coeff = np.array([1., 1., 1., 1.]), ld_power = np.array([0., 0.5, 1., 1.5]));
			oh *= L_perarea[ind];
	
			s += oh;
	
	print times[m], s;
	light_blocked = np.append(light_blocked, s);

	'''

	for ind in range(0,len(ta)):
		if L_perarea[ind] != 0.:
			c_trunc = c[ind:len(ta)];
			sr = c_trunc[0][0];
			yi = getTangentsIntersections(c_trunc);
			m2 = generateRadiiThetas(n,yi[0], yi[1], yi[2]);
			
			f = rd2(m2, c_trunc, opt = 0);
			oh= groupAndIntegrate(bounds = f, num = n, star_rad = sr, ld_coeff = np.array([1., 0., 0., 0.]), ld_power = np.array([0., 0.5, 1., 1.5]));
			oh *= L_perarea[ind];
		
			s += oh;
	print s, times[m];
	light_blocked = np.append(light_blocked, s);

	'''
plt.clf();
plt.plot(times, -light_blocked, '.');
plt.xlim(times[0], times[len(times) - 1]);
plt.show();
