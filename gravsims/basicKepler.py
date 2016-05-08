import numpy as np;
from scipy import stats, constants;
from scipy.optimize import fsolve;
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

STEPS = 7501;
REVOLUTIONS = 0.125;

class OrbitingSystem:
	#add more orbiting systems in, with mass of 0
	#list the most massive body in the system first
	mass = 0;
	semimajor = 0;
	arg_periastron = 0;
	inclination = 0;
	eccentricity = 0;
	phase = 0;
	radius = 0;
	temperature = 0;
	
	bt = 1;
	
	bodies = [];
	orbits = [np.array([])];
	times = [];
	
	"""semimajor_to_main = []; #first element should always be 0 or self makes no sense
	arguments_to_main = []; #in every one of these lists
	inclinations_to_main = [];
	eccentricities_to_main = [];"""
	
	def setTotalMass(self):
		m = 0;
		for p in self.bodies:
			m += p.mass;
		
		self.mass = m;
		
	def setBaseTime(self):
		base_time = getTrueAnomaly(self.bodies[0].mass, self.bodies[1].mass, self.bodies[1].semimajor, self.bodies[1].eccentricity)[1];
		self.bt = base_time;
				
	def setSystemOrbits(self):
		base_time = self.bt;
		print len(self.bodies);
		print base_time/86400;
		basemass_x, basemass_y, basemass_z = 0, 0, 0;
		self.orbits = list(np.zeros(len(self.bodies)));
		for satellite in range(1, len(self.bodies)):
			a = getTrueAnomaly(self.bodies[0].mass, self.bodies[satellite].mass, self.bodies[satellite].semimajor, self.bodies[satellite].eccentricity);
			anom = a[0]; 
			time_taken = a[1];
			
			print a[1]/(60 * 60 * 24);
			
			ot = getOrbitTimes(anom, t=base_time, phase = self.bodies[satellite].phase, scale_factor = base_time/time_taken); 
			self.times = ot[0];
			#gets a time- and phase- dependent set of radii and angles, needed eventually for transits
			q = computeOrbit(self.bodies[0].mass, self.bodies[satellite].mass, self.bodies[satellite].semimajor, self.bodies[satellite].eccentricity, ot[1], ot[2]);
			
			base_xyz = transformOrbit(q[0], ot[1], self.bodies[satellite].inclination, self.bodies[satellite].arg_periastron);
			
			#find out if += works for NumPy arrays, highly doubt so
			#mathematically: seems legit?
			basemass_x = basemass_x + base_xyz[0];
			basemass_y = basemass_y + base_xyz[1];
			basemass_z = basemass_z + base_xyz[2];
			
			planet_xyz = transformOrbit(q[1], ot[2], self.bodies[satellite].inclination, self.bodies[satellite].arg_periastron);
			self.orbits[satellite] = (np.array(planet_xyz)/149597870700.);
			
		self.orbits[0] = np.array([basemass_x, basemass_y, basemass_z])/149597870700.;
		satellite = 0;

def getTrueAnomaly(m1, m2, a, e):
	#Masses of each body are in Earth masses
	
	G = constants.G;
	m1 *= 5.9726 * 10**24;
	m2 *= 5.9726 * 10**24;
	a *= 149597870700.;
	
	t = 2 * np.pi * np.sqrt((a)**3/(G * (m1 + m2)));
	apoapsis_distance = a/(1 - e);
	periapsis_distance = a/(1 + e);
	
	#NOTE that true anomaly and mean anomaly are flipped in the terminology here
	#But this does not impact the veracity of the orbits generated, only a name change is needed
	def true_anomaly(theta, a):
		ma = 2 * np.arctan(np.sqrt((1 - e)/(1 + e)) * np.tan(theta/2));
		ta = ma - e * np.sin(ma);
		
		return ta - a;
		
	return (true_anomaly, t);
				
def getOrbitTimes(true_anomaly, t, phase = 0, scale_factor = 1):
	#equal areas in equal times --> dA ~ dt
	#making sure that times do not become slightly asynchronous due to floating point error, need our dAs right.

	#Note that phase is in terms of mean anomaly (should it be?) rather than true
	dA = 0;
	dA = np.linspace(phase, (2 * REVOLUTIONS * np.pi) * int(scale_factor * STEPS)/STEPS + phase, STEPS)-np.pi;

	theta_m1 = np.array([]);
	for p in range(0, len(dA)):
		angle = dA[p]
		val = fsolve(true_anomaly, x0=angle, args=((angle + np.pi) % (2 * np.pi) - np.pi));
		theta_m1 = np.append(theta_m1, val);
	#theta_m1 = np.array([fsolve(true_anomaly, angle, (angle + np.pi) % (2 * np.pi) - np.pi) for angle in dA]);
		
	theta_m2 = np.pi + theta_m1;
	
	dt = (dA+np.pi -phase)/np.max(dA+np.pi-phase) * REVOLUTIONS * t/86400;
	
	tm1_f, tm2_f, dt_f = theta_m1, theta_m2, dt;
	
	return np.array([dt_f, tm1_f, tm2_f]);

	
def computeOrbit(m1, m2, a, e, theta_m1, theta_m2):
	G = constants.G;
	m1 *= 5.9726 * 10**24;
	m2 *= 5.9726 * 10**24;
	a *= 149597870700.

	r1 = a * m2/(m1 + m2) * (1-e**2)/(1 + e * np.cos(theta_m1));
	r2 = a * m1/(m1 + m2) * (1-e**2)/(1 - e * np.cos(theta_m2));
	
	return np.array([r1, r2]);
	
def transformOrbit(r, theta, inclination, arg_periastron):
	xa0 = r * np.cos(theta)
	ya0 = r * np.sin(theta)
	
	xa = xa0 * np.cos(inclination);
	ya = ya0
	
	z = xa0 * np.sin(inclination);	

	x = xa * np.cos(arg_periastron) - ya * np.sin(arg_periastron);
	y = xa * np.sin(arg_periastron) + ya * np.cos(arg_periastron);
	
	return x, y, z;

def whole_orbit(m1, m2, a, e, arg_periastron, inclination=0):
	anom, times = getTrueAnomaly(m1, m2, a, e);
	got = getOrbitTimes(anom, times, time_selection=0);
	ol = computeOrbit(m1, m2, a, e, got[1], got[2]);

	r = np.linspace(-2 * a * 149597870700, 2 * a * 149597870700);
	l = a * 149597870700
	
	xa, ya, za = transformOrbit(ol[0], got[1], inclination, arg_periastron);
	xb, yb, zb = transformOrbit(ol[1], got[2], inclination, arg_periastron);
	
	plot_orbit(a, np.array([xa, ya, za]), np.array([xb, yb, zb]));
	
	return np.array([xa, ya, za]), np.array([xb, yb, zb]);



def plot_orbit(a, args):

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


"""JUST A SAMPLE DEMONSTRATION, SHOULD NOT BE EXECUTED
	
q = OrbitingSystem();
q.bodies = [0,0];

q.bodies[0] = OrbitingSystem();
q.bodies[0].mass = 332946.0487 #0.8 includes inner planets
q.bodies[0].radius = 695700./149597870.700;

q.bodies[1] = OrbitingSystem();
q.bodies[1].inclination = np.pi/2;
q.bodies[1].semimajor =  149598023./149597870.700;
q.bodies[1].bodies = [0,0];

q.bodies[1].bodies[0] = OrbitingSystem();
q.bodies[1].bodies[0].mass = 1.
q.bodies[1].bodies[1] = OrbitingSystem();
q.bodies[1].bodies[1].mass = 0.0123031469;
q.bodies[1].bodies[1].eccentricity = 0.0549;
q.bodies[1].bodies[1].inclination = np.pi/2;
q.bodies[1].bodies[1].semimajor = 384748./149597870.700;
q.bodies[1].setTotalMass();
q.bodies[1].eccentricity = 0.0167086;

'''q.bodies[2] = OrbitingSystem();
q.bodies[2].mass = 317.83;
q.bodies[2].semimajor = 778570000./149597870.700;
q.bodies[2].eccentricity = 0.0489;
q.bodies[2].inclination = 1.304 * np.pi/180.;

q.bodies[2].bodies = [OrbitingSystem()]
q.bodies[2].mass = 317.83;'''



q.setBaseTime();
q.bodies[0].bt, q.bodies[1].bt, q.bodies[1].bodies[0].bt, q.bodies[1].bodies[1].bt = q.bt, q.bt, q.bt, q.bt;
#q.bodies[2].bt, q.bodies[2].bodies[0].bt = q.bt, q.bt;

q.setSystemOrbits();
q.bodies[1].setSystemOrbits();
#q.bodies[2].setSystemOrbits();
	
starPlanet = q;
times = starPlanet.times;

final_array = []

for a in range(0, len(starPlanet.bodies)):
	for b in range(0, len(starPlanet.bodies[a].bodies)):
		print (a, b)
		x = starPlanet.bodies[a].orbits[b][0] + starPlanet.orbits[a][0] - starPlanet.bodies[a].orbits[0][0] * 0 - starPlanet.orbits[0][0] * 0;
		y = starPlanet.bodies[a].orbits[b][1] + starPlanet.orbits[a][1] - starPlanet.bodies[a].orbits[0][1] * 0 - starPlanet.orbits[0][1] * 0;
		z = starPlanet.bodies[a].orbits[b][2] + starPlanet.orbits[a][2] - starPlanet.bodies[a].orbits[0][2] * 0 - starPlanet.orbits[0][2] * 0;
		
		final_array.append(np.array([x, y, z]));
		
		plt.plot(times,np.sqrt(x**2 + y**2));
		
plt.ylabel("Difference from mean (AU)");
plt.xlabel("Time (days)");
plt.title("Sky-projected deviation of bodies in planetary system from COM");
#plt.xlim(times[0], times[len(times) - 1]);
plt.xlim(q.bt/(2 * 60 * 60 * 24 * 365) - 1, q.bt/(2 * 60 * 60 * 24 * 365) + 1);
plt.ylim(0, q.bodies[0].radius);
		
plt.show();

plot_orbit(a=5, args = final_array);


	

#OPERATING EXAMPLE TO RUN THIS CODE
starPlanet = OrbitingSystem();
starPlanet.bodies = [0,0,0,0]

q = OrbitingSystem();
q.bodies = [0,0];

q.bodies[0] = OrbitingSystem()
q.bodies[0].mass = 332946.-125000.;

q.bodies[1] = OrbitingSystem();
q.bodies[1].mass = 125000.;
q.bodies[1].semimajor = 0.05;

q.setTotalMass();

q.mass = 332946.;

r = OrbitingSystem();
r.bodies = [0,0,0,0,0]
for n in range(0,len(r.bodies)):
	o = OrbitingSystem();
	o.mass = 4**(-n);
	o.semimajor = 0.001 * n;
	o.eccentricity = n/(n + 3);
	o.arg_periastron = n * np.pi/(n + 1);
	o.inclination = n * np.pi/(n + 5) * 0.5;
	r.bodies[n] = o;
	
r.semimajor = 1.0;
r.eccentricity = 0.5;
r.inclination = np.pi/6
r.setTotalMass();
	
s = OrbitingSystem();
s.bodies = [0,0,0,0,0]
for m in range(0,len(s.bodies)):
	o = OrbitingSystem();
	o.mass = 2**(-m) * 10;
	o.semimajor = 0.01 * m;
	o.eccentricity = m/(m + 3);
	o.arg_periastron = m * np.pi/(m + 1);
	o.inclination = m * np.pi/(m + 5) * 0.5;
	s.bodies[m] = o;

s.semimajor = 1.5;
s.eccentricity = 0.7;
s.inclination = np.pi/5
s.setTotalMass();

t = OrbitingSystem();
t.bodies = [0,0,0,0,0]
for m in range(0,len(t.bodies)):
	o = OrbitingSystem();
	o.mass = 2**(-m) * 100;
	o.semimajor = 0.01 * m;
	o.eccentricity = m/(m + 3);
	o.arg_periastron = m * np.pi/(m + 1);
	o.inclination = m * np.pi/(m + 5) * 0.5;
	t.bodies[m] = o;

t.semimajor = 1.5;
t.eccentricity = 0.7;
t.phase = np.pi;
t.arg_periastron = 2 * np.pi/3.;
t.setTotalMass();


 
starPlanet.bodies[0], starPlanet.bodies[1], starPlanet.bodies[2], starPlanet.bodies[3] = q, r, s, t;
starPlanet.setBaseTime();
year = starPlanet.bt
starPlanet.bodies[0].bt, starPlanet.bodies[1].bt, starPlanet.bodies[2].bt, starPlanet.bodies[3].bt = year, year, year, year;
starPlanet.setSystemOrbits();
times = starPlanet.times;

for p in starPlanet.bodies[1:]:
	p.setSystemOrbits();

final_array = []

for a in range(0, len(starPlanet.bodies)):
	for b in range(0, len(starPlanet.bodies[a].bodies)):
		print (a, b)
		x = starPlanet.bodies[a].orbits[b][0] + starPlanet.orbits[a][0] - starPlanet.bodies[a].orbits[0][0] * 0 - starPlanet.orbits[0][0] * 0;
		y = starPlanet.bodies[a].orbits[b][1] + starPlanet.orbits[a][1] - starPlanet.bodies[a].orbits[0][1] * 0 - starPlanet.orbits[0][1] * 0;
		z = starPlanet.bodies[a].orbits[b][2] + starPlanet.orbits[a][2] - starPlanet.bodies[a].orbits[0][2] * 0 - starPlanet.orbits[0][2] * 0;
		
		final_array.append(np.array([x, y, z]));
		
		plt.plot(times,np.sqrt(x**2 + y**2));
		
plt.ylabel("Difference from mean (AU)");
plt.xlabel("Time (days)");
plt.title("Sky-projected deviation of bodies in planetary system from COM");
plt.xlim(times[0], times[len(times) - 1]);
		
plt.show();

plot_orbit(a=1.5, args = final_array);

"""
