import numpy as np;
import itertools;
import matplotlib.pyplot as plt;
import time;
import copy;

def realign(circ):
	circles = np.array(circ);
	circles.T[1] -= circles[0][1];
	circles.T[2] -= circles[0][2];
	return circles;

def solveCircle(c, a, b, theta):
	'''theta_constant = np.arctan2(a, b);
	theta_variable = np.pi;
	if c < np.sqrt(a**2 + b**2):
		theta_variable = np.arcsin(c/np.sqrt(a**2 + b**2));
	
	f = np.zeros(len(theta));
	g = np.zeros(len(theta));
	
	k = ((theta > (theta_constant - theta_variable) % (2 * np.pi)) | (theta < (theta_constant + theta_variable) % (2 * np.pi))) & (theta_constant > (theta_constant + theta_variable) % (2 * np.pi));
	k = (k | ((theta > (theta_constant - theta_variable) % (2 * np.pi)) & (theta < (theta_constant + theta_variable) % (2 * np.pi))) & (theta_constant < (theta_constant + theta_variable) % (2 * np.pi)));
	print len(k) - len(theta);'''
	
	f = a * np.sin(theta) + b * np.cos(theta);
	g = np.sqrt(f**2 - (a**2 + b**2 - c**2))
	
	r = f + g;
	s = f - g;
	r[np.isnan(r) | (r < 0)] = 0.;
	s[np.isnan(s) | (s < 0)] = 0.;
	
	return np.array([r, s]);
    
def solveIntersections(plots, star_rad):
	plots2 = np.copy(plots);
	pairs = itertools.combinations(np.arange(1,len(plots)), 2);
	for m in pairs:
		#if a point is within another body, remove it from consideration
		plots[m[0]][0][(plots2[m[0]][0] < plots2[m[1]][0]) & (plots2[m[0]][0] > plots2[m[1]][1])] = 0.;
		plots[m[0]][1][(plots2[m[0]][1] < plots2[m[1]][0]) & (plots2[m[0]][1] > plots2[m[1]][1])] = 0.;
		plots[m[1]][0][(plots2[m[1]][0] < plots2[m[0]][0]) & (plots2[m[1]][0] > plots2[m[0]][1])] = 0.;
		plots[m[1]][1][(plots2[m[1]][1] < plots2[m[0]][0]) & (plots2[m[1]][1] > plots2[m[0]][1])] = 0.;
		
	for u in plots:
		u[0][u[0] > star_rad] = star_rad;
		u[1][u[1] > star_rad] = star_rad;
		
	upper = np.sum([a[0] for a in plots[1:]], axis=0);
	lower = np.sum([b[1] for b in plots[1:]], axis=0);
	
	upper[upper > star_rad] = star_rad;
	lower[lower > star_rad] = star_rad

	#print (plots == plots2);	
		
	return plots[1:];
	
def integrate(plots, ld_coeffs, ld_powers, star_rad, theta):
	ld_powers, ld_coeffs = np.array(ld_powers), np.array(ld_coeffs);
	d_theta = (theta[-1] - theta[0])/(len(theta) - 1)
	answer = 0.;

	for m in range(0,len(ld_powers)):
		for p in plots:
			integral = -np.trapz((star_rad**2 - p[0][p[0] != star_rad]**2)**(ld_powers[m]/2. + 1.)) + np.trapz((star_rad**2 - p[1][p[1] != star_rad]**2)**(ld_powers[m]/2. + 1.))
			answer += integral * ld_coeffs[m] / star_rad**(ld_powers[m] + 2.)
				
	answer /= np.pi * np.sum(ld_coeffs) * 2.;
	answer *= d_theta;
	
	return answer;
				
'''

import fastIntegrator as fi
import time
import numpy as np
import matplotlib.pyplot as plt

lc = [];
e = time.time();
n = 100
q = 2000;
for t in range(0, q):
	cir = [[5., 0., 0.],[1., -6. + 12 * t/1000., 1.],[2.,-7. + 14. * t/1000,2. * 2. * t/1000.]]
	g = fi.realign(cir)
	theta = np.linspace(0., 2. * np.pi, 2 * np.pi * int(g.T[0][0]/g.T[0].min() * n));
	m = [fi.solveCircle(a[0], a[1], a[2], theta) for a in g]
	k = fi.solveIntersections(m, 5.)
	dip = fi.integrate(k, np.array([1., 1., 1.]), np.array([0., 1., 2.]), 5., theta);
	lc.append(dip);
lc = np.array(lc);
print (time.time() - e)/q;

'''
