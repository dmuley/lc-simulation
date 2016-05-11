if __name__ == '__main__':
    if __package__ is None:
		import sys;
		from os import path;
		sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) );
		print sys.path[-1];
		from integrators.polarIntegrator3 import *;
    else:
        from integrators.polarIntegrator3 import *;
	
import numpy as np;
from itertools import chain;
import matplotlib.pyplot as plt;
import itertools;
from integrators.polarIntegrator3 import *;

def traverse_tree(starPlanet, times):
	""" Starting off with a base case, this algorithm traces the entire tree of orbits
		to an arbitrary depth in order to obtain final x, y, and z-positions of every
		body in the system. A more elegant and natural interpretation than the previous
		forced model will be added later. """
	
	final_x, final_y, final_z = [], [], [];
	transit_array = [];
		
	"""

	UNTESTED JAGGED ALGORITHM. DO NOT USE.

	final_x = [np.zeros(len(times))];
	final_y = [np.zeros(len(times))];
	final_z = [np.zeros(len(times))];
	transit_array = [starPlanet];

	indices_array = np.array([0]);
	#local_indices_array = np.array([0]);
	
	position = 0;
	for u in range(0,5):
		#Have to do a few test runs to ensure that all layers processed.
		length = len(transit_array);
		for v in range(position, length):
			if (transit_array[v].bodies != []):
				transit_array.extend(transit_array[v].bodies);
				final_x.extend(np.array([transit_array[v].orbits[a][0] for a in range(0, len(transit_array[v].orbits))]));
				final_y.extend(np.array([transit_array[v].orbits[b][1] for b in range(0, len(transit_array[v].orbits))]));
				final_z.extend(np.array([transit_array[v].orbits[c][2] for c in range(0, len(transit_array[v].orbits))]));

				indices_array = np.append(indices_array, np.ones(len(transit_array[v].orbits)) * v);
				#local_indices_array = np.append(local_indices_array, np.arange(0,len(transit_array[v].orbits));

		position += length - position;

	m = 0:
	
	valid_search_indices = np.array([]);

	for item in range(0,len(transit_array)):
		if (transit_array[item].bodies != []):
			valid_search_indices = np.append(valid_search_indices, item);

	for ind in valid_search_indices:
		newpos = ind;
		while (newpos != 0):
			newpos = indices_array[newpos];

			final_x[ind] += final_x[newpos];
			final_y[ind] += final_y[newpos];
			final_z[ind] += final_z[newpos];
	
	final_x = final_x[valid_search_indices];
	final_y = final_y[valid_search_indices];
	final_z = final_z[valid_search_indices];
	transit_array = transit_array[valid_search_indices]; """

	
	for a in range(0, len(starPlanet.bodies)):
		for b in range(0, len(starPlanet.bodies[a].bodies)):
			print (a, b)
			x = starPlanet.bodies[a].orbits[b][0] + starPlanet.orbits[a][0];
			y = starPlanet.bodies[a].orbits[b][1] + starPlanet.orbits[a][1];
			z = starPlanet.bodies[a].orbits[b][2] + starPlanet.orbits[a][2];
		
			final_x.append(x);
			final_y.append(y);
			final_z.append(z);

			starPlanet.bodies[a].orbits[b][0] = x;
			starPlanet.bodies[a].orbits[b][1] = y;
			starPlanet.bodies[a].orbits[b][2] = z;
		
			transit_array.append(starPlanet.bodies[a].bodies[b]);
		
			plt.plot(times, y/np.absolute(y) * np.sqrt(x**2 + y**2), '-');

	plt.ylabel("Difference from mean (AU)");
	plt.xlabel("Time (days)");
	plt.title("Sky-projected deviation of bodies in planetary system from COM");
	plt.xlim(times[0], times[len(times) - 1]);
	plt.show();
	
	return final_x, final_y, final_z, transit_array;
	
				
	
###### LOOM ALGORITHM ######
# Using the Z-values we will index the orders of X and Y values, and then get transits for each body.

def sort_keys(final_x, final_y, final_z):
	fx = np.array(final_x).T;
	fy = np.array(final_y).T;
	fz = np.array(final_z).T;

	zpos = np.argsort(fz);
	print len(zpos)
	
	return fx, fy, fz, zpos;

	### measuring intersections

def arrange_combinations(fx, fy, transit_array):
	ta_elements = np.arange(0,len(transit_array));
	ta_combos = itertools.combinations(list(ta_elements), 2);
	intersects = (fx * 0).T;

	print len(ta_elements), len(intersects);

	for a in ta_combos:
		if (transit_array[a[0]].temperature != 0) or (transit_array[a[1]].temperature != 0):
				intersects[a[0]][np.where((fx.T[a[0]] - fx.T[a[1]])**2 + (fy.T[a[0]] - fy.T[a[1]])**2 < (transit_array[a[0]].radius + transit_array[a[1]].radius)**2)] += 1;
				intersects[a[1]][np.where((fx.T[a[0]] - fx.T[a[1]])**2 + (fy.T[a[0]] - fy.T[a[1]])**2 < (transit_array[a[0]].radius + transit_array[a[1]].radius)**2)] += 1;

	#	print np.where((fx.T[a[0]] - fx.T[a[1]])**2 + (fy.T[a[0]] - fy.T[a[1]])**2 < (transit_array[a[0]].radius + transit_array[a[1]].radius)**2);


	light_blocked = np.zeros(len(fx));
	
	return light_blocked, intersects;

#### ACTUAL TRANSIT MODELING BELOW ####

#Number of transit-modeling integration steps per timestep

def generate_lightcurve(fx, fy, light_blocked, intersects, n, ta, times, zpos):	
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
				
		print str(times[m]) + "\t" + str(t);
	
		light_blocked[m] = t;
	
	lb_2 = light_blocked
	return times, lb_2;

