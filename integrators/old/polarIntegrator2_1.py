#A cleaner version of the polar numeric integrator.
#Dhruv Muley, 2016

import numpy as np;
from scipy.optimize import fsolve;
import matplotlib.pyplot as plt;
import time;
import itertools;

def getTangentsIntersections(circ):
	circles = np.array(circ);
	circles.T[1] -= circles[0][1];
	circles.T[2] -= circles[0][2];
	
	t_circles = circles.T;
	
	tans = np.array([]);
	intersections = [];
	radii = [];
	circles_associated = [];
	
	tangents = np.zeros(len(circles)) + np.pi;
	tangents_var = np.zeros(len(circles)) + np.pi
	tangents[t_circles[1] != 0.] = np.arctan2(t_circles[2][t_circles[1] != 0.],t_circles[1][t_circles[1] != 0.]);
	tangents_var[t_circles[2]**2 + t_circles[1]**2 != 0] = np.arcsin(t_circles[0][t_circles[2]**2 + t_circles[1]**2 != 0]/np.sqrt(t_circles[2]**2 + t_circles[1]**2)[t_circles[2]**2 + t_circles[1]**2 != 0]);
	
	tangents = np.concatenate((tangents - tangents_var + 5 * np.finfo(np.float32).eps, tangents, tangents + tangents_var - 5 * np.finfo(np.float32).eps));
		
		#maybe will be used later, but not now.
		#radii.append([constant_radius - variable_radius, constant_radius + variable_radius]);
		

	
		#now to get intersections
	for m in range(0,len(circles)):
		for s in range(m, len(circles)):
			if (s != m):
				#using law of cosines
				dist_1 = np.sqrt((circles[m][1] - circles[s][1])**2 + (circles[m][2] - circles[s][2])**2);
				dist_2 = circles[m][0];
				dist_3 = circles[s][0];
				
				#dist_3**2 = dist_1**2 + dist_2**2 - 2(dist_1)(dist_2)cos(theta)
				
				try:
					theta_constant = np.arctan((circles[s][2] - circles[m][2])/(circles[s][1] - circles[m][1]));
				except ZeroDivisionError:
					theta_constant = np.pi;
				
				theta_variable = np.arccos((dist_1**2 + dist_2**2 - dist_3**2)/(2 * dist_1 * dist_2))
				
				x0 = dist_2 * np.cos(theta_constant - theta_variable) + circles[m][1];
				y0 = dist_2 * np.sin(theta_constant - theta_variable) + circles[m][2];
				
				x1 = dist_2 * np.cos(theta_constant + theta_variable) + circles[m][1];
				y1 = dist_2 * np.sin(theta_constant + theta_variable) + circles[m][2];
				
				intersect_1 = np.arctan2(y0,x0) + 0 * np.finfo(np.float32).eps;
				intersect_2 = np.arctan2(y1,x1) - 0 * np.finfo(np.float32).eps;
				
				if ((not np.isnan(intersect_1)) or (not np.isnan(intersect_2))):
					intersections.append([intersect_1, intersect_2]);
				else:
					intersections.append([0.,np.pi * 2])
				#circles_associated.append([m, s]);
					
				
	tangents = np.array(tangents);
	radii = np.array(radii);	
	intersections = np.array(intersections);
	circles_associated = np.array(circles_associated);
		
	return (circles, tangents, intersections, radii, circles_associated);
	
def generateRadiiThetas(n, circles, *args):
	#tt = time.time();
	r = np.array([]);
	for item in args:
		r = np.append(r, item.ravel());
	
	r = np.unique(np.concatenate((r, r + np.pi/2, r + np.pi, r + 3 * np.pi/2, 2 * np.pi - r, np.pi - r)));
	
	r[(r > np.pi * 2) | (r < 0)] %= (2 * np.pi);

	
	#print len(r);
	bounded_region = np.array([]);
	for o in range(0,len(r) - 1):
		ang = np.linspace(r[o], r[o + 1], n);
		bounded_region = np.append(bounded_region, ang);
	#print "Linspace Time: " + str(time.time() - tt);
		
	#print len(bounded_region);
		
	angles = np.unique(bounded_region);

	radii = np.zeros((len(circles) * 4, len(angles)));
	#print len(radii[0]), len(angles);
	
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
		
	
	#print time.time() - tt;	
		
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
		
		l = time.time() - h;
		#print l
			
		if opt == 1:
			return l;
		if opt == 0:
			return coords[0][indicator < 1], np.vstack([coords[1]] * len(coords[0][0])).T[indicator < 1];
		if opt == 2:
			coords[0][indicator >= 1] -= 1000;
			oh = coords[0];
			u = [np.concatenate((np.unique(oh[a][oh[a] == 0.]), oh[a][(oh[a] > 0.) & (oh[a] < oh[a].max())], np.unique(oh[a][oh[a] == oh[a].max()]))) for a in range(0,len(coords[1]))];
		
			return u, coords[1];
				
def groupAndIntegrate(bounds, num):
	rad = [];
	theta = [];

	for a in range(0,len(np.unique(bounds[1]))):
		q = bounds[0][bounds[1] == np.unique(bounds[1])[a]];
		if (len(np.unique(q)) % 2 == 0 and len(q) > 0):
			rad.append(np.unique(q));
			theta.append(np.unique(bounds[1])[a]);

	theta = np.array(theta);	
	d_theta = np.diff(theta);

	#need to tighten up loose bounds by grouping them
	#provably (except in tangent line degenerate case or floating point error) there are always
	#an even number of intersections
	
	area = 0;
	
	for i in range(0,len(d_theta))[::-1]:
		h = 0;
		if (d_theta[i] < 2 * np.pi / (num - 1.)):
			h += np.sum(rad[i][1::2]**2 - rad[i][0::2]**2);
			h += np.sum(rad[i - 1][1::2]**2 - rad[i - 1][0::2]**2);

			h *= d_theta[i]**2/np.sin(d_theta[i]);
		
		area += h;
	
	area /= 4.;
			
	
	return area;

def plot_circles(circ):
	for q in circ:
		o = np.linspace(-q[0], q[0], 101);
		p = np.sqrt(q[0]**2 - o**2);
		
		plt.plot(o + q[1], p + q[2]);
		plt.plot(o + q[1], q[2] - p);
		
def plot_tangent_lines(tangents, circ):
	for q in tangents:
		o = np.linspace(-circ[0][0], circ[0][0],201);
		plt.plot(o * np.cos(q) + circ[0][1], o * np.sin(q) + circ[0][2]);
				
"""		
u = time.time();
c = [[100., 0., 0.],[1,98.,-2.]]; #no need for dummy here
y = getTangentsIntersections(c);
print time.time() - u;
u = time.time();
m = generateRadiiThetas(16,y[0], y[1], y[2]);
print time.time() - u;
f = rd2(m, c, opt = 0);
oh= groupAndIntegrate(f, 15);
print "Light output (normalized) = " + str(oh);
print time.time() - u;
time_f = 0.;
#for q in range(0,100):
#	oh = rd2(m, c,opt=1);
#	time_f += oh;
#print "Final time (100 iterations): " + str(time_f/100.);
plot_circles(c);
#plot_tangent_lines(m[1], c);
plt.title("Bodies (position of star = (0,0); radius of star = " + str(c[0][0]) + "\n normalized light output = " + str(1. - oh/(c[0][0]**2 * np.pi)) + ")");
plt.xlabel("sky-projected x-position (arbitrary units)");
plt.ylabel("sky-projected y-position (arbitrary units)");
for ee in range(0,len(f[0])):
    plt.plot(f[0][ee] * np.cos(f[1][ee]), f[0][ee] * np.sin(f[1][ee]), '.');
	
plt.ylim(-2,-4);
plt.xlim(96,100);
plt.show();
"""