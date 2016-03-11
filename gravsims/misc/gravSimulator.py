import numpy as np
from scipy import stats, constants
import matplotlib.pyplot as plt
import scipy.integrate as integrate

NUMOFPLANETS = 20;
RBOUNDS = (-5. *150000000000, 5. *150000000000);
VBOUNDS = (-1 * 1000, 1 * 1000);
MAX_MASS = 1.989e27;
GRAVCONST = 6.67428e-11;
DIST_MINIMUM = 0.01 * 150000000000;

TIMESTEP = 400000;

class Planet:
	def __init__(self):
		pass;
	def setMassRadius(self, mass=1., radius=1.):
		self.mass = mass;
		self.radius = radius;
		
	#Set motions of everything in here
	def setXY(self, x, y, z):
		self.x = x;
		self.y = y;
		self.z = z;
	def setVelocity(self, vx=0., vy=0., vz=0.):
		self.vx = vx;
		self.vy = vy;
		self.vz = vz;
	def setAcceleration(self, ax=0., ay=0., az=0.):
		self.ax = ax;
		self.ay = ay;
		self.az = az;
	
#FUNCTIONS THAT HELP OUR PROGRAM ARE LISTED BELOW#
#												 #
##################################################

def createPlanetsList():
	lPlanets = []
	for a in range(0, NUMOFPLANETS):
		t = Planet();
		t.setMassRadius(np.random.uniform(0, MAX_MASS));
		t.setXY(x=np.random.uniform(RBOUNDS[0], RBOUNDS[1]), y=np.random.uniform(RBOUNDS[0], RBOUNDS[1]), z=np.random.uniform(RBOUNDS[0], RBOUNDS[1]));
		t.setVelocity(np.random.uniform(VBOUNDS[0], VBOUNDS[1]),np.random.uniform(VBOUNDS[0], VBOUNDS[1]),np.random.uniform(VBOUNDS[0], VBOUNDS[1]));
		t.setAcceleration();
		
		lPlanets.append(t);
	return lPlanets;

def getGravitationalAccel(listOfPlanets):
	xAccels = np.array([])
	yAccels = np.array([])
	zAccels = np.array([])

	for planet in listOfPlanets:
		xComps = 0.;
		yComps = 0.;
		zComps = 0.;
		for o in listOfPlanets:
			if (planet != o):
				distSquared = ((planet.x - o.x)**2 + (planet.y - o.y)**2 + (planet.x - o.z)**2);
				if distSquared < DIST_MINIMUM:
					distSquared = DIST_MINIMUM;
					
					
				fG = GRAVCONST * (planet.mass * o.mass)/distSquared;
				gX = -fG * (planet.x - o.x)/np.sqrt(distSquared);
				gY = -fG * (planet.y - o.y)/np.sqrt(distSquared);
				gZ = -fG * (planet.z - o.z)/np.sqrt(distSquared);
				
				
				xComps += gX;
				yComps += gY;
				zComps += gZ;
			else:
				xComps += 0;
				yComps += 0;
				zComps += 0;
		xAccels = np.append(xAccels, xComps/planet.mass);
		yAccels = np.append(yAccels, yComps/planet.mass);
		zAccels = np.append(zAccels, zComps/planet.mass);
		
		#print xAccels, yAccels
		
	return np.array([xAccels, yAccels, zAccels]);
	
def setGravitationalAccel(listOfPlanets, listOfAccels):
	for a in range(0, len(listOfAccels[0])):
		listOfPlanets[a].ax = listOfAccels[0][a];
		listOfPlanets[a].ay = listOfAccels[1][a];
		listOfPlanets[a].az = listOfAccels[2][a];
	
def movePlanets(listOfPlanets):
	u = 0;
	for body in listOfPlanets:
		body.x += body.vx * TIMESTEP + body.ax * TIMESTEP**2/2.;
		body.vx += body.ax * TIMESTEP;
		body.y += body.vy * TIMESTEP + body.ay * TIMESTEP**2/2.;
		body.vy += body.ay * TIMESTEP;
		body.z += body.vz* TIMESTEP + body.az * TIMESTEP**2/2.;
		body.vz += body.az * TIMESTEP;
		
		#print (body.vx * TIMESTEP + body.ax * TIMESTEP**2/2., body.vy * TIMESTEP + body.ay * TIMESTEP**2/2.);
	#print "";

l = createPlanetsList();
x = [];
y = [];
z = [];
for f in np.arange(0, TIMESTEP * 2000, TIMESTEP):
	acc = getGravitationalAccel(l);
	setGravitationalAccel(l, acc);
	movePlanets(l);
	x.append(np.array([p.x for p in l]));
	y.append(np.array([t.y for t in l]));
	z.append(np.array([q.y for q in l]));
	print f;
	
x = np.array(x);
y = np.array(y);
z = np.array(z);

x = np.transpose(x);
y = np.transpose(y);
z = np.transpose(z);

for a in range(0, NUMOFPLANETS):
    plt.plot(x[a], y[a]);
    plt.plot(x[a][len(x[a]) - 1], y[a][len(y[a]) - 1], '.', markersize = np.sqrt(l[a].mass/5.972e25) * 5);
    
plt.show()