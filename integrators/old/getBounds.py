import numpy as np;
import scipy;
import matplotlib.pyplot as plt;

'''This program is designed to work with ***any*** number of bodies, for any possible limb-darkening profile. The only condition (which ensures that integration can take place with analytic pieces, not numerical squares) is that the light profile have a circular footprint and be radially symmetrical).'''

class MovingBody:
	radius = 1.;
	x = 0.;
	y = 0.;
	z = 0.;

	'''Can get data from N-body simulator and plug it into here. To transform from N-body into telescope view,
	simply set the telescope position to be arbitrarily far away (say, 10000 AU) from the system.'''

	#functions of radial distance from the center.
	#Possible to get an analytic solution over here.
	coefficient_series = np.array([]);
	power_series = np.array([]);

	def setRadius(rad):
		this.radius = rad;

	def getRadius():
		return this.radius;

	def setXYZ(x0, y0, z0=0):
		this.x = x0;
		this.y = y0;
		this.z = z0;

	def moveXY(Zdx, dy, dz=0;):
		this.x += dx;
		this.y += dy;
		this.z += dz;

	def getXYZ():
		return (this.x, this.y, this.z);

	#Gets tangent lines from other body.
	def getTangentLines(moving):
		x1 = moving.getXYZ()[0];
		y1 = moving.getXYZ()[1];

		r = moving.getRadius();

		theta1 = np.arctan((this.y - y1)/(this.y - y1));
		theta2 = np.arcsin(this.radius/np.sqrt((this.y - y1)**2 + (this.x - x1)**2));

		return np.array(theta1 - theta2, theta1 + theta2);

	def getRadialRestrictions(moving):
 		x1 = moving.getXYZ()[0];
		y1 = moving.getXYZ()[1];

		r = moving.getRadius();

		radius_inner = np.sqrt((this.x - x1)**2-(this.y - y1)**2) - this.radius;
		radius_outer = np.sqrt((this.x - x1)**2-(this.y - y1)**2) + this.radius;

		if radius_inner < 0.:
			radius_inner = 0.;

		if radius_outer > r:
			radius_outer = r;

		return np.array(radius_inner, radius_outer);

class AnnularSector:
	#base case is the unit circle

	radius_inner = 0.;
	radius_outer = 1.;
	theta_lower = 0.;
	theta_upper = 2 * np.pi;

	a = 0.;
	b = 0.;

	def setRadialRestriction(inner, outer);
		this.radius_inner = inner;
		this.radius_outer = outer;

	def setAngularRestriction(lower, upper):
		this.theta_lower = lower;
		this.theta_upper = upper;

	def setPosition(x, y);
		this.a = x;
		this.b = y;

	def getParameters():
		return np.array([this.radius_inner, this.radius_outer, this.theta_lower, this.theta_upper, this.a, this.b]);

	#intersection needs to be checked with respect to EVERY obscuring body. This will be done in later functions.
	def intersectsCircle(moving_body):
		intersects = false;
		
		#checks if each part of the restriction (two curvilinear regions, two straight lines) intersect the circle
		#Getting all variables out of the way here
		x_mb = moving_body.getXYZ()[0];
		y_mb = moving_body.getXYZ()[1];

		pos_x = x_mb - this.a;
		pos_y = x_mb - this.b;

		ang_0 = this.theta_lower;
		ang_1 = this.theta_upper;

		p = moving_body.getRadius();

		r_inner = this.radius_inner;
		r_outer = this.radius_outer;	
		
		#First: is it entirely outside of the circle? If so, the *surface within* needs to be marked for splitting
		#Avoid dealing with negative angles, so add 2pi radians to everything
		
		tangentLines = moving_body.getTangentLines();
		radialRestrictions = moving_body.getRadialRestrictions();

		if ((ang_0 <= tangentLines[0] and ang_1 >= tangentLines[1]) and (radialRestrictions[0] >= r_inner and radialRestrictions[1] <= r_outer)):
			intersects = true;
			return intersects;

		#If one of the lines is intersected, then the overall figure is intersected
		def r(theta, x1, y1, rad):
			center = (x1 * np.cos(theta) + y1 * np.sin(theta));
			branch = np.sqrt(center**2 - (x1**2 + y1**2 - rad**2));

			return np.array([center - branch, center + branch]);

		for angle_restriction in [ang_0, ang_1]:
			if (tangentLines[0] <= angle_restriction <= tangentLines[1]):
				if ((r_inner <= r(angle_restriction, pos_x, pos_y, p)[0] <= r_outer) or (r_inner <= r(angle_restriction, pos_x, pos_y, p)[1] <= r_outer)):
					except TypeError:
						pass;
					intersects = true;
					return intersects;

		#if one of the bounding curves is intersected, then the overall figure is intersected
		#use law of cosines here, very simple solution

		def intersection_angles(r_large, r_small, x1, y1):
			center = np.arctan(y1/x1);
			branch = np.arccos((x1**2 + y1**2 + r_large**2 - r_small**2)/(2 * r_large * np.sqrt(x1**2 + y1**2)));
			
			return np.array([center - branch, center + branch]);

		for radius_restriction in [r_inner, r_outer]:
			q = intersection_angles(r_outer, r_inner, pos_x, pos_y):
			for item in q:
				if (ang_0 <= q <= ang_1):
					except TypeError:
						pass;
					intersects = true;
					return intersects;
		
		#Are there any other ways of gauging intersection of the circle?

		#if it doesn't intersect, the function will return false
		return intersects;

	def checkWithinOutside(listOfMovingBodies):
		'''Check to see if there's an intersection. If the restriction region intersects even ONE circle, play it cautious and split it up anyway. If there are no intersections, measure one point and see if it is within at least ONE of the circles; if that is true, then all of the points are within the circle, and no splitting is necessary and the area can be added to the running total (in the helper function). If that one point is outside all circles, then (since it doesn't intersect) no points are within the circle, and that restriction can be discarded.'''

	def splitRestriction():
		#Simple function, just splits the restriction into 4 additional restrictions. These can be passed back through the helper function to see if *they* ned splitting, and so on.

def sortListOfBodies(listOfMovingBodies):
	q = np.array([]);
	obj = [];
	for p in listOfMovingBodies:
		q = np.append(q, p.getXYZ()[2]);

	ind = np.argsort(q);
	a = 0;
	for n in ind:
		obj[n] = listOfMovingBodies[na;
		a++;

	return obj;

def getExcludedRegion(listOfMovingBodies): #has to be sorted, or else won't work correctly. currently only works with 1st moving body.
	#grouped according to which body they're from. body index = restriction index.
	radialRestrictions = [];
	angleRestrictions = [];

	#not grouped, to be used for ordering
	rr = np.array([]);
	ang= np.array([]);
	for i in listOfMovingBodies[1:]:
		radialRestrictions.append(i.getRadialRestrictions(listOfMovingBodies[0]));
		angleRestrictions.append(i.getTangentLines(listOfMovingBodies[0]));

		rr = np.append(rr, i.getRadialRestrictions(listOfMovingBodies[0]));
		ang =np.array(ang, i.getTangentLines(listOfMovingBodies[0]));

	#come up with a way of consistently ordering angles (don't want excessively long or short domains of integration!)
	#need to pick min and max for range

	a = np.sort(ang + 2 * np.pi)

	a0 = a[0];
	a1 = a[len(a) - 1];

	#Also, get minimum and maximum radii.
	#this is a trivial task for the other function.

	r = np.sort(rr);
	
	r0 = r[0];
	r1 = r[len(r) - 1];

	x = listOfMovingBodies[0].getXYZ[0];
	y = listOfMovingBodies[0].getXYZ[1];

	return np.array([[r0, r1], [a0, a1], [x, y]);

#Need to create a function here to integrate
