import numpy as np;
from scipy import stats, constants;
from scipy.optimize import fsolve;
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

STEPS = 10000; #should always be an int
REVOLUTIONS = 1.; #should always be a float

class OrbitingSystem:
	"""	OrbitingSystem has a wide range of attributes that correspond
		to orbital parameters of stellar systems such as inclination,
		argument of periastron, inclination, eccentricity, et cetera.
		Particularly importantly, its bodies attribute is an array
		which can contain more OrbitingSystems, and its orbits attribute
		stores the orbits of these child OrbitingSystems around their
		mutual center of mass.

		The first body in the bodies array of each OrbitingSystem is
		the "star". One can set a mass or radius for the body, but it
		cannot be directly given any of the other orbital parameters. The
		final shape of its orbit is entirely dependent on the orbital
		parameters of the succeeding bodies in the system.

		To achieve a retrograde orbit, add 180 degrees to the inclination
		of the body in question, and add 180 degrees to its argument of
		periastron. """
	#add more orbiting systems in, with mass of 0
	#list the most massive body in the system first
	mass = 0;
	"""Earth masses"""
	semimajor = 0;
	"""AU"""
	arg_periastron = 0;
	inclination = 0;
	phase = 0;
	"""Radians"""
	eccentricity = 0;
	radius = 0;
	"""AU"""
	temperature = 0;
	"""Arbitrary units"""
	
	bt = 1;
	
	ld_coeffs = [1];
	ld_powers = [0];
	bodies = [];
	orbits = [np.array([])];
	times = [];
	
	def setTotalMass(self):
		""" Based on the masses of each system component, can calculate the mass of the
		whole OrbitingSystem. Please use this to avoid guesswork. Now updated to work
		recursively along the lines of the traverse_tree code in general_transits.py."""
		m = 0;

		transit_array = [self];
		
		indices_array = np.zeros(1).astype('int');

		position = 0;
		for u in range(0, 5):
			length = len(transit_array);
			for v in range(position, length):
				if (transit_array[v].bodies != []):
					transit_array.extend(transit_array[v].bodies);
					indices_array = np.append(indices_array, np.ones(len(transit_array[v].bodies)) * v);

			position += length - position;

		#print len(transit_array);
		masses = np.array([body.mass if body.bodies == [] else 0. for body in transit_array])
		valid_search_indices = np.array([a for a in range(0,len(transit_array)) if (transit_array[a].bodies == [])]).astype('int');

		for ind in valid_search_indices:
			ind = int(ind);
			newpos = ind;

			while (newpos != 0):
				#traces UP the tree, rather than down.
				newpos = int(indices_array[newpos]);
				masses[newpos] += masses[ind];

		for m in range(0,len(transit_array)):
			transit_array[m].mass = masses[m];
		
	def setBaseTime(self):
		"""Sets the size of each timestep to ensure that the orbits of all bodies
		and compared to one another."""
		base_time = getMeanAnomaly(self.bodies[0].mass, self.bodies[1].mass, self.bodies[1].semimajor, self.bodies[1].eccentricity)[1];
		self.bt = base_time;
				
	def setSystemOrbits(self, s = STEPS, r = REVOLUTIONS, verbose = False):
		"""Sets the orbit of each body in the system based on the base timestep
		and the orbital parameters given. The central body responds to forces by
		the orbiting bodies, but other orbiting bodies do not."""

		base_time = self.bt;
		if verbose == True:
			print len(self.bodies);
			print base_time/86400;
		basemass_x, basemass_y, basemass_z = 0, 0, 0;
		self.orbits = list(np.zeros(len(self.bodies)));
		for satellite in range(1, len(self.bodies)):
			a = getMeanAnomaly(self.bodies[0].mass, self.bodies[satellite].mass, self.bodies[satellite].semimajor, self.bodies[satellite].eccentricity);
			anom = a[0]; 
			time_taken = a[1];
			
			if verbose == True:
				print a[1]/(60 * 60 * 24);
			
			ot = getOrbitTimes(anom, t=base_time, phase = self.bodies[satellite].phase, scale_factor = base_time/time_taken, STEPS=s, REVOLUTIONS = r); 
			self.times = ot[0];
			#gets a time- and phase- dependent set of radii and angles, needed eventually for transits
			q = computeOrbit(self.bodies[0].mass, self.bodies[satellite].mass, self.bodies[satellite].semimajor, self.bodies[satellite].eccentricity, ot[1], ot[2]);
			
			base_xyz = transformOrbit(q[0], ot[1], self.bodies[satellite].inclination, self.bodies[satellite].arg_periastron);
			
			#find out if += works for NumPy arrays, highly doubt so
			#mathematically: seems legit?
			basemass_x += base_xyz[0];
			basemass_y += base_xyz[1];
			basemass_z += base_xyz[2];
			
			planet_xyz = transformOrbit(q[1], ot[2], self.bodies[satellite].inclination, self.bodies[satellite].arg_periastron);
			self.orbits[satellite] = (np.array(planet_xyz)/149597870700.);
			
		self.orbits[0] = np.array([basemass_x, basemass_y, basemass_z])/149597870700.;
		satellite = 0;

	def insertBody(self, position = [],
        
		mass = 0,
		semimajor = 0,
		arg_periastron = 0,
		inclination = 0,
		phase = 0,
		eccentricity = 0,
		radius = 0,
		temperature = 0,
		bt = 1,
		ld_coeffs = [1],
		ld_powers = [0]):

		'''This function will allow one to dynamically insert a body at a given tree depth
		 and position by appending it to the bodies list of an existing orbiting system. 
		 It neither computes orbits for the overall system or adjusts existing orbits, 
		 nor generates a light curve.'''

	    	newBody = OrbitingSystem();
	       	newBody.mass = mass;
	      	newBody.semimajor = semimajor;
	      	newBody.arg_periastron = arg_periastron;
	       	newBody.inclination = inclination;
	       	newBody.phase = phase;
	       	newBody.eccentricity = eccentricity;
	       	newBody.radius = radius;
	       	newBody.temperature = temperature;
	       	newBody.bt = bt;
	       	newBody.ld_coeffs = ld_coeffs;
	       	newBody.ld_powers = ld_powers;

		if (position == []):
			self.bodies.append(newBody);
		else:
			m = self;
			for i in position:
				m = m.bodies[i];
            
        		m.bodies.append(newBody);

	def modifyBody(self, position = [],

		mass = 1j,
		semimajor = 1j,
		arg_periastron = 1j,
		inclination = 1j,
		phase = 1j,
		eccentricity = 1j,
		radius = 1j,
		temperature = 1j,
		bt = 1j,
		ld_coeffs = 1j,
		ld_powers = 1j):

		'''This function enables a body's properties, excluding its total mass and sub-bodies, which can
		be set later, to be altered dynamically while a run is being performed. One needs to get the right
		position, however, to determine which body to affect.'''

		newBody = self;

		if (position != []): 
			for i in position:
				newBody = newBody.bodies[i];

		if (mass != 1j):
	        	newBody.mass = mass;
	        if (semimajor != 1j):
			newBody.semimajor = semimajor;
	        if (arg_periastron != 1j):
			newBody.arg_periastron = arg_periastron;
		if (inclination != 1j):
	        	newBody.inclination = inclination;
	        if (phase != 1j):
			newBody.phase = phase;
	        if (eccentricity != 1j):
			newBody.eccentricity = eccentricity;
	        if (radius != 1j):
			newBody.radius = radius;
	        if (temperature != 1j):
			newBody.temperature = temperature;
	        if (bt != 1j):
			newBody.bt = bt;
	        if (ld_coeffs != 1j):
			newBody.ld_coeffs = ld_coeffs;
	        if (ld_powers != 1j):
			newBody.ld_powers = ld_powers;

	def removeBody(self, position = []):
		'''Removes a body from the OrbitingSystem given its position. It does not, however
		alter the total mass of the body or any other attribute.'''

		if (position == []):
			#Do nothing, should not delete the main body.
			print "No action performed."

		else:
			m = self;
			for i in position[:-1]:
				m = m.bodies[i];
			m.bodies.pop(position[-1]);

def getMeanAnomaly(m1, m2, a, e):
	"""Obtains the mean anomaly (circularized change in angle) as a function of the true
	angular anomaly. Requires solving Kepler's equation, which is seen in mean_anomaly
	child function. Often times, fsolve is the slowest link in the chain, but problems are
	negligible even on a large scale."""

	#Masses of each body are in Earth masses
	
	G = constants.G;
	m1 *= 5.9726 * 10**24;
	m2 *= 5.9726 * 10**24;
	a *= 149597870700.;
	
	t = 2 * np.pi * np.sqrt((a)**3/(G * (m1 + m2)));
	apoapsis_distance = a/(1 - e);
	periapsis_distance = a/(1 + e);
	
	def mean_anomaly(theta, a):
		ma = 2 * np.arctan(np.sqrt((1 - e)/(1 + e)) * np.tan(theta/2));
		ta = ma - e * np.sin(ma);
		
		return ta - a;
		
	return (mean_anomaly, t);
				
def getOrbitTimes(mean_anomaly, t, phase = 0, scale_factor = 1, STEPS = STEPS, REVOLUTIONS = REVOLUTIONS):
	"""Uses Kepler's Second Law of equal areas in equal times (which is true of the mean anomaly) in order to
	calculate the time at which a body has a given true anomaly. This is then used, along with the orbital radius
	at that true anomaly, to determine the shape of the orbit.

	Performs rounding to ensure that the right number of timesteps is used."""
	#equal areas in equal times --> dA ~ dt
	#making sure that times do not become slightly asynchronous due to floating point error, need our dAs right.

	#Note that phase is in terms of mean anomaly (should it be?) rather than true
	dA = 0;
	dA = np.linspace(phase, (2 * REVOLUTIONS * np.pi) * int(scale_factor * STEPS)/STEPS + phase, STEPS)-np.pi;

	theta_m1 = np.array([]);
	for p in range(0, len(dA)):
		angle = dA[p]
		val = fsolve(mean_anomaly, x0=angle, args=((angle + np.pi) % (2 * np.pi) - np.pi));
		theta_m1 = np.append(theta_m1, val);
		
	theta_m2 = np.pi + theta_m1;
	
	dt = (dA+np.pi -phase)/np.max(dA+np.pi-phase) * REVOLUTIONS * t/86400;
	
	tm1_f, tm2_f, dt_f = theta_m1, theta_m2, dt;
	
	return np.array([dt_f, tm1_f, tm2_f]);

	
def computeOrbit(m1, m2, a, e, theta_m1, theta_m2):
	"""Uses the angles of the true anomaly in order to compute the actual distances
	involved in the orbit of each body around the center of mass."""

	G = constants.G;
	m1 *= 5.9726 * 10**24;
	m2 *= 5.9726 * 10**24;
	a *= 149597870700.

	r1 = a * m2/(m1 + m2) * (1-e**2)/(1 + e * np.cos(theta_m1));
	r2 = a * m1/(m1 + m2) * (1-e**2)/(1 - e * np.cos(theta_m2));
	
	return np.array([r1, r2]);
	
def transformOrbit(r, theta, inclination, arg_periastron):
	"""Transforms the orbit according to inclination and argument of periastron."""
	x0 = r * np.cos(theta)
	y0 = r * np.sin(theta)
	
	'''xa = xa0 * np.cos(inclination);
	ya = ya0
	
	z = xa0 * np.sin(inclination);	
	'''
	x = x0 * np.cos(arg_periastron) - y0 * np.sin(arg_periastron);
	y = x0 * np.sin(arg_periastron) + y0 * np.cos(arg_periastron);
	
	x *= np.cos(inclination);
        z = x * np.sin(inclination);


	return x, y, z;

def whole_orbit(m1, m2, a, e, arg_periastron, inclination=0):
	"""Combines all other functions for orbits into one. Used in the
	getSystemOrbits class function for each body in the system."""

	anom, times = getMeanAnomaly(m1, m2, a, e);
	got = getOrbitTimes(anom, times, time_selection=0);
	ol = computeOrbit(m1, m2, a, e, got[1], got[2]);

	r = np.linspace(-2 * a * 149597870700, 2 * a * 149597870700);
	l = a * 149597870700
	
	xa, ya, za = transformOrbit(ol[0], got[1], inclination, arg_periastron);
	xb, yb, zb = transformOrbit(ol[1], got[2], inclination, arg_periastron);
	
	plot_orbit(a, np.array([xa, ya, za]), np.array([xb, yb, zb]));
	
	return np.array([xa, ya, za]), np.array([xb, yb, zb]);



def plot_orbit(a, args):
	"""Plots every orbit three-dimensionally. Useful for visualization purposes
	and sanity checks, but without much further scientific use. For scientific
	purposes a 2-D visual (distance vs. time) should be generated."""
	l = a
	fig = plt.figure();
	ax = fig.add_subplot(111, projection='3d');
	for obj in args:
		ax.plot(xs = obj[0], ys = obj[1], zs = obj[2]);
	
	ax.plot([0.], [0.], [0.], markerfacecolor='k', markeredgecolor='k', marker='o', markersize=5, alpha=1)
	ax.set_xlim(-l, l);
	ax.set_ylim(-l, l);
	ax.set_zlim(-l, l);
	
	fig.show();
