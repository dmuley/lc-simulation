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

###### LOOM ALGORITHM ######
# Using the Z-values we will index the orders of X and Y values, and then get transits for each body.

def sort_keys(final_x, final_y, final_z):
	fx = np.array(final_x).T;
	fy = np.array(final_y).T;
	fz = np.array(final_z).T;

	zpos = np.argsort(fz);

	print len(fx), len(times);
	
	return fx, fy, fz, zpos;

	### measuring intersections

def arrange_combinations(fx, fy, transit_array):
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
	
	return light_blocked, intersects;

def generate_lightcurve(fx, fy, light_blocked, intersects, n, ta, times):	
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
	
	lb_2 = light_blocked
	return times, lb_2;

