#Later on we can add more checks in here, like creating a feasible region for integration.
#For now, it's just a proof-of-concept, so there is no need.

import numpy as np;
from scipy.optimize import fsolve;
import matplotlib.pyplot as plt;
import time;

def getTangentsIntersections(circ):
	cir = np.array(circ);
	circles = circ;

	x_base = circles[0][1];
	y_base = circles[0][2];
	
	for q in circles:
		q[1] -= x_base;
		q[2] -= y_base;
		
	tangents = [];
	radii = [];
	intersections = [];
	circles_associated = [];
	
	#get tangent lines and minimal/maximal radii of the circle with respect to the basis
	
	for m in range(0, len(circles)):
		try:
			constant_angle = np.arctan(circles[m][2]/circles[m][1]);
		except ZeroDivisionError:
			constant_angle = np.pi;

		try:
			variable_angle = np.arcsin(circles[m][0]/np.sqrt(circles[m][1]**2 + circles[m][2]**2));
		except ZeroDivisionError:
			variable_angle = np.pi;
			
		if (np.isnan(constant_angle) or np.isnan(variable_angle)):
			constant_angle = np.pi;
			variable_angle = np.pi;

		constant_radius = np.sqrt(circles[m][1]**2 + circles[m][2]**2);
		variable_radius = circles[m][0];
		
		if (np.abs(constant_radius) < np.abs(variable_radius)):
			constant_angle = np.pi;
			variable_angle = np.pi;
			
		#correction for floating-point error
		tangents.append([constant_angle - variable_angle+0.0000001, constant_angle + variable_angle-0.0000001]);
		
		#maybe will be used later, but not now.
		#radii.append([constant_radius - variable_radius, constant_radius + variable_radius]);
		

	
		#now to get intersections
	
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
				
				intersect_1 = np.arctan(y0/x0);
				intersect_2 = np.arctan(y1/x1);
				
				if ((not np.isnan(intersect_1)) or (not np.isnan(intersect_2))):
					intersections.append([intersect_1, intersect_2]);
				else:
					intersections.append([0.,np.pi * 2])
				circles_associated.append([m, s]);
					
				
	tangents = np.array(tangents);
	radii = np.array(radii);	
	intersections = np.array(intersections);
	circles_associated = np.array(circles_associated);
	
	circ = cir;
	
	return (cir, tangents, intersections, circles_associated);
	
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
				
def get_integration_thetas(n, *args):
	r = np.array([])
	for item in args:
		r = np.append(r, item.ravel());
		
	print len(r)
	print r;

	print len(r);
	r = np.array(list(set(r)));
	region = np.array([]);
	print len(r);
	#for item in r:
	#	region = np.append(region, np.array([normalize(item), normalize(item + np.pi)]));
	region = r;
	print len(region);

	region = np.sort(region);
		
	bounded_region = np.array([]);
	for o in range(0,len(region) - 1):
		angles = np.linspace(region[o], region[o + 1], n);
		bounded_region = np.append(bounded_region, angles);
		
	bounded_region = np.sort(np.array(list(set(bounded_region))));
	print len(bounded_region);
	return bounded_region;
	
def normalize(angle):
	if ((angle/(2 * np.pi)) % 1 == 0 and angle != 0):
		angle = 2 * np.pi;
	else:
		angle = angle % (2 * np.pi);
		
	return angle;
	
def get_integration_radii(integration_thetas, c):
	x_base = c[0][1];
	y_base = c[0][2];
	
	for q in c:
		q[1] -= x_base;
		q[2] -= y_base;
	
	points = [];
	leaves = [];
	thetas = [];
	for theta in integration_thetas:
		points.append([0., theta]);
		leaves.append(theta);
		
		append_r1 = False;
		append_r2 = False;
		
		for circle in c[1:]:

			a = 1.;
			b = -2. * (circle[1] *np.cos(theta) + circle[2] * np.sin(theta))
			ci = circle[1]**2 + circle[2]**2 - circle[0]**2
			
			radius_1 = (-b + np.sqrt(b**2 - 4.*a*ci))/(2. * a);
			radius_2 = (-b - np.sqrt(b**2 - 4.*a*ci))/(2. * a);
			
			"""Up to a maximum of 3 points per iteration."""
			
			if 0 <= radius_1 < c[0][0]:
				points.append([radius_1, normalize(theta)]);
				leaves.append(normalize(theta));
			elif 0 > radius_1 > -c[0][0]:
				points.append([-radius_1, normalize(theta + np.pi)]);
				leaves.append(normalize(theta + np.pi));
					
			if 0 <= radius_2 < c[0][0]:
				points.append([radius_2, normalize(theta)]);
				leaves.append(normalize(theta));
			elif 0 > radius_2 > -c[0][0]:
				points.append([-radius_2, normalize(theta + np.pi)]);
				leaves.append(normalize(theta + np.pi));
				
			if abs(radius_1) >= c[0][0] or abs(radius_2) >= c[0][0]:
				points.append([c[0][0], normalize(theta)]);
				leaves.append(normalize(theta));
				
				if (radius_1/radius_2)/(abs(radius_1/radius_2)):			
					points.append([c[0][0], normalize(theta + np.pi)]);
					leaves.append(normalize(theta + np.pi));
					
	
	print len(points);						
		
	return np.array([points,leaves]);
	
def isInside(points, leaves, circles_integrated):
	circles = circles_integrated;
	insides = np.array([])
	h = time.time();
	for pt in points:
		x = pt[0] * np.cos(pt[1]);
		y = pt[0] * np.sin(pt[1]);
		
		indicator = 0;
		
		for m in circles[1:]:
			if ((x - m[1])**2 + (y - m[2])**2 < (m[0] * 0.9999)**2) and abs(abs(pt[0]) - circles[0][0]) > 0.00001:
			#if inside the circles, within some tolerance
				if (x != 0 or y != 0):
					indicator += 1;
				elif (x == 0 and y == 0):
					indicator -= 1000;
			elif (x == 0 and y == 0):
					indicator += 1;
			
			elif abs(abs(pt[0]) - circles[0][0]) < 0.00001 and (x - m[1])**2 + (y - m[2])**2 > (m[0]*1.0001)**2:
				indicator += 1./(len(circles) - 1);
				

		insides = np.append(insides, indicator);
	
	o = np.where(insides < 1);
	pts = points[o];
	lvs = leaves[o];
	
	print time.time() - h;
	
	f = time.time();
	pts_final = [];
	lvs_final = [];
	
	print len(pts), len(lvs);
	
	for n in range(0, len(pts)):
		inside = False;
		if len(pts_final) > 0:
			for i in pts_final:
				if (abs(i[0] - pts[n][0]) < 0.000001 and abs(i[1] - pts[n][1]) < 0.000001):
					inside = True;
					continue;
					
		if inside == False:
			pts_final.append(list(pts[n]));
			lvs_final.append(lvs[n]);
	
	print len(pts_final), len(lvs_final);
	pts_final = np.array(pts_final);
	lvs_final = np.array(lvs_final);
	
	h = np.argsort(lvs_final)
	
	print (time.time() - f);
	
	return [pts_final[h], lvs_final[h]];

#groups functions by lines and eliminates extraneous lines
def group(points, leaves, circles):
	base_time = time.time();
	leaves = np.array([normalize(a) for a in leaves]);
	leaves_set = np.unique(leaves);
	grouped_radii = [];
	grouped_radii_2 = [];
	
	for unique_leaf in leaves_set:
		u = np.where(leaves == unique_leaf);
		to_append = np.array([a[0] for a in points[u]])
		if 0. not in to_append:
			to_append = np.append(0., to_append);
		
		grouped_radii.append(list(to_append));
		
	for rad in range(0,len(grouped_radii)):
		radii = grouped_radii[rad];
		append_array = [];
		for r in range(0,len(radii) - 1):
			midpoint = (radii[r] + radii[r + 1])/2.;
			to_append = False;
			
			x = midpoint * np.cos(leaves_set[rad]);
			y = midpoint * np.sin(leaves_set[rad]);
			
			to_append = False;
			for m in circles[1:]:
				if (x - m[1])**2 + (y - m[2])**2 < (m[0] * 1.001)**2:
					to_append = True;
					
			if to_append == True:
				append_array.append([radii[r], radii[r + 1]]);	
				to_append = False;
	
		grouped_radii_2.append(append_array);
	
	print time.time() - base_time
	return grouped_radii_2, leaves_set;

	
#test code to get one instance; will test later with many

u = time.time();
c = [[5., 0., 0.],[1,-1,1],[3,3,3]];
y = getTangentsIntersections(c);
m = get_integration_thetas(9, y[1], y[2]);
intrad = get_integration_radii(m, c);
p_inside = isInside(intrad[0], intrad[1], c);
p_grouped = group(p_inside[0], p_inside[1], c);


u -= time.time();
print (-u);

plot_circles(c)
#plot_tangent_lines(m, c)

h = np.transpose(p_inside[0]);

plt.plot(h[0] * np.cos(h[1]), h[0] * np.sin(h[1]), 'b.');


for o in range(0, len(p_grouped[0])):
	for p in range(0, len(p_grouped[0][o])):
		#print p_grouped[0][o][p];
		ls = np.linspace(p_grouped[0][o][p][0], p_grouped[0][o][p][1], 101);
		#print p_grouped[1][o];
		plt.plot(ls * np.cos(p_grouped[1][o]), ls * np.sin(p_grouped[1][o]));
		plt.plot(p_grouped[0][o][p][0] * np.cos(p_grouped[1][o]), p_grouped[0][o][p][0] * np.sin(p_grouped[1][o]), '.');
		plt.plot(p_grouped[0][o][p][1] * np.cos(p_grouped[1][o]), p_grouped[0][o][p][1] * np.sin(p_grouped[1][o]), '.');

    
plt.plot(0,0); 

plt.ylim(-10,10);
plt.xlim(-10,10);
plt.title("Bounds of integration for n-body intersection");
plt.xlabel("x (R_earth)");
plt.ylabel("y (R_earth)");
plt.show()