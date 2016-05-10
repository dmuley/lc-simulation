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

#STEPS = 1001;
#REVOLUTIONS = 0.125;

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


###### LOOM ALGORITHM ######
# Using the Z-values we will index the orders of X and Y values, and then get transits for each body.

fx = np.array(final_x).T;
fy = np.array(final_y).T;
fz = np.array(final_z).T;

zpos = np.argsort(fz);

print len(fx), len(times);

### measuring intersections

ta_elements = np.arange(0,len(transit_array));
ta_combos = itertools.combinations(list(ta_elements), 2);
intersects = (fx * 0).T;

print len(ta_elements), len(intersects);

for a in ta_combos:
	if (transit_array[a[0]].temperature != 0 or transit_array[a[1]].temperature != 0):
       		intersects[a[0]][np.where((fx.T[a[0]] - fx.T[a[1]])**2 + (fy.T[a[0]] - fy.T[a[1]])**2 < (transit_array[a[0]].radius + transit_array[a[1]].radius)**2)] += 1;
       		intersects[a[1]][np.where((fx.T[a[0]] - fx.T[a[1]])**2 + (fy.T[a[0]] - fy.T[a[1]])**2 < (transit_array[a[0]].radius + transit_array[a[1]].radius)**2)] += 1;

#	print np.where((fx.T[a[0]] - fx.T[a[1]])**2 + (fy.T[a[0]] - fy.T[a[1]])**2 < (transit_array[a[0]].radius + transit_array[a[1]].radius)**2);


light_blocked = np.zeros(len(fx));

#### ACTUAL TRANSIT MODELING BELOW ####

#Number of transit-modeling integration steps per timestep
n = 31;
ta = transit_array;

for m in range(0,len(light_blocked)):
	#Excluding bodies that do not intersect at this frame
	#Need >0.1 because some bodies without any intersections are noted to have -0 intersections (off by some epsilon).
	c = np.array([[ta[r].radius * 1000., fx[m][r] * 1000., fy[m][r] * 1000.] for r in zpos[m]])[intersects.T[m][zpos[m]] > 0.1];
	lum = np.array([ta[s].temperature for s in zpos[m]])[intersects.T[m][zpos[m]] > 0.1];
	
	t = 0.;
	
	#Excluding as stars bodies with zero luminosity
	for i in np.where(lum != 0.)[0]:
		ct = c[i:len(c)];
		ct.T[1] -= ct[0][1];
		ct.T[2] -= ct[0][2];

		ti = getTangentsIntersections(ct);

		if not ((ti[1] % np.pi) < np.zeros(len(ti[1])) + 0.01).all():
			rt = generateRadiiThetas(n, ti[0], ti[1], ti[2]);
			f = rd2(rt, ct, opt = 0);
			oh = groupAndIntegrate(bounds = f, num = n, star_rad = ct[0][0], ld_coeff = [1.9, 0., 1., 0.], ld_power = [0.0, 0.5, 1., 1.5]);
			
			del rt, f;
			oh *= lum[i]
			t += oh;
				
	print times[m], t
	
	light_blocked[m] = t;

#Plotting final light curve.
plt.clf();
plt.plot(times, -light_blocked, '-');
plt.xlim(times[0], times[len(times) - 1]);
plt.show();
