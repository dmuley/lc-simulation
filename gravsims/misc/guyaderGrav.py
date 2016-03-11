#uses the method of Le Guyader (1993), which I independently
#derived

import numpy as np;
from scipy import constants, misc;
import matplotlib.pyplot as plt;
from itertools import combinations;

STEPS = 10;
LEN_TIMESTEP = 1000.;
NUM_BODIES = 25.;
DIST_MINIMUM = 0.01;

list_of_bodies = [];

class OrbitingBody:
	mass = 10000.;

	uid = np.random.rand();
	basic = np.zeros(10) * 1.;
	

	x_variables = np.array([basic]);
	y_variables = np.array([basic]);
	z_variables = np.array([basic]);

	#x_variables[0] is the 0th time derivative of x, i. e. x itself
	#ergo with the others

	def setUID(self):
		self.uid = np.random.rand();
	def setTimeDerivative(self, x, y, z, order):
		self.x_variables[0][order] = x;
		self.y_variables[0][order] = y;
		self.z_variables[0][order] = z;
		
		#print (self.x_variables[0][order], self.y_variables[0][order], self.z_variables[0][order]);
		
	def stepTimeForward(self):
		self.x_variables = np.vstack([self.basic, self.x_variables]);
		self.y_variables = np.vstack([self.basic, self.y_variables]);
		self.z_variables = np.vstack([self.basic, self.z_variables]);
		
	def updateLowerDerivatives(self):
		self.stepTimeForward();
		#lowest 2 derivatives, get them directly since they have
		#real, physical antecedents
		
		for n in range(0,2):
			for m in range(0, len(self.basic) - n):
				recurrence = 1./misc.factorial(m) * 1. * LEN_TIMESTEP**m;
				self.x_variables[0][n] += recurrence * self.x_variables[1][m+n];
				self.y_variables[0][n] += recurrence * self.y_variables[1][m+n];
				self.z_variables[0][n] += recurrence * self.z_variables[1][m+n];
	
	'''def updateGravitation(self, g_accels):
		#set gravitational accelerations		
		gr = g_accel;
		#print gr;
		self.x_variables[0][2] = gr[0];
		self.y_variables[0][2] = gr[1];
		self.z_variables[0][2] = gr[2];
		'''
		
	def updateHigherDerivatives(self):
		
		#update higher derivatives
		for n in range(2, len(self.basic)):
			self.x_variables[0][n] = (self.x_variables[0][n-1] - self.x_variables[1][n-1])/LEN_TIMESTEP;
			self.y_variables[0][n] = (self.y_variables[0][n-1] - self.y_variables[1][n-1])/LEN_TIMESTEP;
			self.z_variables[0][n] = (self.z_variables[0][n-1] - self.z_variables[1][n-1])/LEN_TIMESTEP;
			
def getGravAcceleration(listOfBodies):
	c = range(0,len(listOfBodies));
	cx = np.zeros(len(listOfBodies));
	cy = np.zeros(len(listOfBodies));
	cz = np.zeros(len(listOfBodies));
	
	masses = np.array([m.mass for m in listOfBodies]);
	
	combos = combinations(c, 2); #gets all pairwise interactions
	for i in combos:
		a, b = listOfBodies[i[0]], listOfBodies[i[1]];
		lx = a.x_variables[0][0] - b.x_variables[0][0];
		ly = a.y_variables[0][0] - b.y_variables[0][0];
		lz = a.z_variables[0][0] - b.z_variables[0][0];
		print lx, ly, lz
		
		dist = np.sqrt(lx**2 + ly**2 + lz**2);
		
		if dist<DIST_MINIMUM:
			dist = DIST_MINIMUM;
		
		A = -(constants.G * a.mass * b.mass)/(dist**3);
		
		cx[i[0]] += A * lx/a.mass;
		cx[i[1]] -= A * lx/b.mass;
		
		cy[i[0]] += A * ly/a.mass;
		cy[i[1]] -= A * ly/b.mass;
	
		cz[i[0]] += A * lz/a.mass;
		cz[i[1]] -= A * lz/b.mass;
		
		
		
	return np.transpose(np.vstack([cx, cy, cz]));
		
list_of_bodies = []

for o in np.arange(0, NUM_BODIES):
	body = OrbitingBody();
	body.setUID();
	#print body.uid;
	for n in range(0, 2):
		#print n;
		fx = np.random.rand() * 1. * 10**(-n);
		fy = np.random.rand() * 1. * 10**(-n);
		fz = np.random.rand() * 1. * 10**(-n);
		
		print (fx, fy, fz)
		body.setTimeDerivative(order = n, x = fx, y = fy, z = fz);
		
		print (body.x_variables[0][n], body.y_variables[0][n], body.z_variables[0][n]);
		#print ""
	print  " " 
	list_of_bodies.append(body);
	
n = 0;
for q in range(0,len(list_of_bodies)):
	o = list_of_bodies[q];
	print (o.x_variables[0][n], o.y_variables[0][n], o.z_variables[0][n]);
'''


accels = getGravAcceleration(list_of_bodies);
for r in range(0, len(list_of_bodies)):
		list_of_bodies[r].setTimeDerivative(order = 2, x = accels[r][0], y = accels[r][1], z = accels[r][2]);
		
		#print list_of_bodies[r].x_variables[0][2]
		#print list_of_bodies[r].y_variables[0][2]
		#print list_of_bodies[r].z_variables[0][2];

Xs = [[] for a in np.arange(0,len(list_of_bodies))];
Ys = [[] for a in np.arange(0,len(list_of_bodies))];
Zs = [[] for a in np.arange(0,len(list_of_bodies))];

for q in range(1, STEPS):	
	#print q;
	for r in range(0, len(list_of_bodies)):
		list_of_bodies[r].updateLowerDerivatives();
		
	accels = getGravAcceleration(list_of_bodies);
	#print accels;
	
	for s in range(0, len(list_of_bodies)):
		#print s;
		list_of_bodies[s].setTimeDerivative(order = 2, x = accels[s][0], y = accels[s][1], z = accels[s][2]);
		
		list_of_bodies[s].updateHigherDerivatives();
		
		print (list_of_bodies[s].x_variables[0], list_of_bodies[s].y_variables[0], list_of_bodies[s].z_variables[0]);
	
		Xs[s].append(list_of_bodies[s].x_variables[0][0]);
		Ys[s].append(list_of_bodies[s].y_variables[0][0]);	
		Zs[s].append(list_of_bodies[s].z_variables[0][0]);'''