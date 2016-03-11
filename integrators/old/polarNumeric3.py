#An even faster version of the polar numeric integrator.
#Dhruv Muley, 2016

import numpy as np;
from scipy.optimize import fsolve;
import matplotlib.pyplot as plt;
import time;
import itertools;

def getTangentsIntersections(circ):
	cir = np.array(circ);
	circles = np.array(circ);
	
	x_base = circles[0][1];
	y_base = circles[0][2];
	
	circles.T[1] -= x_base;
	circles.T[2] -= y_base;
	
	#get angles that define lines extending from (0,0) tangent to the circles
	
	constant_angle = np.zeros(len(circles)) + np.pi;
	variable_angle = np.zeros(len(circles)) + np.pi;
	
	valid_range = np.where((circles.T[1] != 0) & np.isfinite(circles.T[1]) & np.isfinite(circles.T[2]));
	v2 = np.where((circles.T[1] != 0) & np.isfinite(circles.T[1]) & np.isfinite(circles.T[2]) & (circles.T[0] <= np.sqrt(circles.T[1]**2 + circles.T[2]**2)));
	
	constant_angle[valid_range] = np.arctan(circles.T[2]/circles.T[1]);
	variable_angle[v2] = np.arcsin(circles.T[0]/np.sqrt(circles.T[1]**2 + circles.T[2]**2));
	
	tangents = np.append(constant_angle + variable_angle - 0.00001, constant_angle - variable_angle + 0.00001);

	#now getting the intersections between each circle
	inter = np.array([a for a in itertools.combinations(np.arange(0,len(circles)), 2)]);
	
	xa = circles[inter.T].T[1].T[0]
	ya = circles[inter.T].T[2].T[0]
	
	xpos = (xa - circles[inter.T].T[1].T[1]);
	ypos = (xa - circles[inter.T].T[2].T[1]);
	
	dist_1 = np.sqrt(xpos**2 + ypos**2);
	dist_2 = circles[inter.T].T[0].T[0];
	dist_3 = circles[inter.T].T[0].T[1];
	
	theta_constant = np.zeros(len(dist_1)) + np.pi;
	theta_variable = np.zeros(len(dist_1)) + np.pi;
	
	valid_theta = np.where((xpos != 0) & (ypos != 0));
	theta_constant[valid_theta] = np.arctan(ypos/xpos);
	theta_variable[valid_theta] = np.arccos((dist_1**2 + dist_2**2 - dist_3**2)/(2 * dist_1 * dist_2));
	
	intersect_1 = np.arctan((dist_2 * np.sin(theta_constant - theta_variable) + xa)/(dist_2 * np.cos(theta_constant - theta_variable) + ya)) + 0.00001;
	intersect_2 = np.arctan((dist_2 * np.sin(theta_constant + theta_variable) + xa)/(dist_2 * np.cos(theta_constant + theta_variable) + ya)) - 0.00001;	
	
	return np.concatenate((intersect_1, intersect_2, tangents), axis=0);			
	

def generateRadiiThetas(n, circles, r):
	tt = time.time();
	r = np.unique(np.nan_to_num(r));
	r[(r > (np.pi * 2)) | (r < 0)] %= (2 * np.pi);

	angles = np.array([]);
	for a in range(0,len(r) - 1):
		angles = np.concatenate((angles, np.linspace(r[a], r[a + 1], n + 1), np.linspace(r[a], r[a + 1], n + 1) + np.pi));

	angles = np.unique(angles);
	print angles;
	
	radii = np.zeros((len(circles) * 4, len(angles)));
	
	for c in range(1, len(circles)):
	
		a = 1.;
		b = -2. * (circles[c][1] *np.cos(angles) + circles[c][2] * np.sin(angles))
		ci = circles[c][1]**2 + circles[c][2]**2 - circles[c][0]**2
			
		radius_1 = (-b + np.sqrt(b**2 - 4.*a*ci))/(2. * a);
		radius_2 = (-b - np.sqrt(b**2 - 4.*a*ci))/(2. * a);
		
		radius_1[radius_1 >= circles[0][0]] = circles[0][0];
		radius_1[radius_1 <= -circles[0][0]] = -circles[0][0];
		radius_2[radius_2 >= circles[0][0]] = circles[0][0];
		radius_2[radius_2 <= -circles[0][0]] = -circles[0][0];
		
		
		radii[4 * c - 4][np.where(radius_1 > 0.)] = radius_1[radius_1 > 0.];
		radii[4 * c - 3][np.where(radius_2 > 0.)] = radius_2[radius_2 > 0.];
		
		offset_1 = angles[radius_1 < 0.] + np.pi;
		offset_2 = angles[radius_2 < 0.] + np.pi;
		
		for o in offset_1:
			q = np.where(np.absolute(angles - o) < 0.00005);
			radii[4 * c - 2][q] = -radius_1[radius_1 < 0.];
			
		for p in offset_2:
			q = np.where(np.absolute(angles - p) < 0.00005)
			radii[4 * c - 1][q] = -radius_2[radius_2 < 0.];			
		
	
	print time.time() - tt;	
		
	return np.transpose(radii), angles;

def rd2(coords, circles,opt = 0):
		h = time.time();

		ind = np.zeros((len(coords[0]), len(coords[0][0])));
		xs = (coords[0].T * np.cos(coords[1]).T).T;
		ys = (coords[0].T * np.sin(coords[1]).T).T;

		
		indicator = np.zeros((len(coords[0]), len(coords[0][0])));
		
		for m in circles[1:]:
				indicator[(((xs - m[1])**2 + (ys - m[2])**2 < (m[0] * 0.9999)**2) & (coords[0] - circles[0][0] < 0.0001)) & ((xs != 0) | (ys != 0))] += 1;
				indicator[(((xs - m[1])**2 + (ys - m[2])**2 < (m[0] * 0.9999)**2) & (coords[0] - circles[0][0] < 0.0001)) & ((xs == 0) & (ys == 0))] -= 1000;
				indicator[(((xs - m[1])**2 + (ys - m[2])**2 < (m[0] * 0.9999)**2) & (coords[0] >= circles[0][0])) & ((xs != 0) | (ys != 0))] -= 1000;
				indicator[np.logical_not((((xs - m[1])**2 + (ys - m[2])**2 < (m[0] * 0.9999)**2) & (abs(coords[0] - circles[0][0]) > 0.0001))) & ((xs == 0.) & (ys == 0.))] += 1;
				
				indicator[((xs - m[1])**2 + (ys - m[2])**2 > (m[0] * 1.0001)**2)] += 1./(len(circles) - 1);	   
		
		#print len(indicator), len(xs);
		l = time.time() - h;
		print l
		
		#print "peeb peeb";
		
		#print coords[0][indicator < 1].shape == np.vstack([coords[1]] * len(coords[0][0])).T[indicator < 1].shape;
		if opt == 1:
			return l;
		else:
			return coords[0][indicator < 1], np.vstack([coords[1]] * len(coords[0][0])).T[indicator < 1];
        
def restrictDomain(coords, circles):
	time_time = time.time()
	
	restricted_radii = [];
	ind = np.zeros((len(coords[0]), len(coords[0][0])));
	xss = (coords[0].T * np.cos(coords[1]).T).T;
	yss = (coords[0].T * np.sin(coords[1]).T).T;
	
	
	print xss.shape == ind.shape;
	print xss.shape == coords[0].shape;
	
	
	for a in range(0,len(coords[0])):
		xs = xss[a];
		ys = yss[a];
		indicator = ind[a];
		
		for m in circles[1:]:
			u = np.where((xs - m[1])**2 + (ys - m[2])**2 < (m[0] * 0.9999)**2);
			v = np.where(abs(coords[0][a] - circles[0][0]) > 0.00001);
			f = np.where((xs - m[1])**2 + (ys - m[2])**2 > (m[0] * 1.0001)**2)
			
			w = np.where(xs == 0);
			w2= np.where(ys == 0);
			
			intersect1= np.intersect1d(u[0], v[0]);
			intersect2= np.intersect1d(w[0], w2[0]);
			
			intersect_comp = np.intersect1d(intersect1, intersect2);
			
			indicator[intersect1] += 1;
			indicator[intersect2] += 1;
			indicator[intersect_comp] -= 100;
			indicator[f] += 1./(len(circles) - 1);
			
		#print indicator;
		#print coords[0][a][indicator < 1]
		restricted_radii.append(np.unique(coords[0][a][indicator < 1]));
				
		#print len(xs);
	print time.time() - time_time;
	return coords[1], restricted_radii;

def plot_circles(circ):
	for q in circ:
		o = np.linspace(-q[0], q[0], 101);
		p = np.sqrt(q[0]**2 - o**2);
		
		plt.plot(o + q[1], p + q[2]);
		plt.plot(o + q[1], q[2] - p);
		
def plot_tangent_lines(tangents, circ):
	for q in tangents:
		o = np.linspace(-10, 10,201);
		plt.plot(o + circ[0][1], o * np.tan(q) + circ[0][2]);
				
		
u = time.time();
c = [[5., 0., 0.],[1,2,0],[2,1,1],[3,2,3],[0,1,1]]; #needs a dummy at the end for some reason
y = getTangentsIntersections(c);
m = generateRadiiThetas(5,c, y);
h = restrictDomain(m, c);
f = rd2(m, c);
print time.time() - u;

time_f = 0.;
#for q in range(0,100):
#	oh = rd2(m, c,opt=1);
#	time_f += oh;

print "Final time (100 iterations): " + str(time_f/100.);

plot_circles(c);
plot_tangent_lines(m[1], c);

#for a in range(0, len(h[0])):
#	plt.plot(h[1][a] * np.cos(h[0][a]), h[1][a] * np.sin(h[0][a]), 'x');

#for a in range(0, len(m[0])):
#	plt.plot(m[0][a] * np.cos(m[1][a]), m[0][a] * np.sin(m[1][a]), '.');

#for a in range(0, len(m[0])):
#	plt.plot(m[0][a] * np.cos(m[1][a]), m[0][a] * np.sin(m[1][a]), '.');

for ee in range(0,len(f[0])):
    plt.plot(f[0][ee] * np.cos(f[1][ee]), f[0][ee] * np.sin(f[1][ee]), '.');
	
plt.ylim(-10,10);
plt.xlim(-10,10);
plt.show();