if __package__ is None:
	import sys;
	from os import path;
	sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) );
	print sys.path[-1];
	from gravsims.basicKepler import *;
else:
	from ...gravsims.basicKepler import *;

import numpy as np;
from itertools import chain;
import matplotlib.pyplot as plt;
import time;


################################################
######## PLACE ALL PLANETARY DATA HERE  ########
################################################
#FOR RETROGRADE ORBITS: Add pi radians to inclination, add pi radians to arg_periastron

#steps = 7501;
#revs = 0.25;

s = OrbitingSystem();
s.bodies = [0,0]
#Masses of stellar system
q = OrbitingSystem();
#q.arg_periastron = np.pi;
q.bodies = [0,0];

q.bodies[0] = OrbitingSystem()
q.bodies[0].mass = 332946. * 0.6897;
q.bodies[0].temperature = 0.985;
q.bodies[0].radius = 0.6489 * 0.00464913034;
q.bodies[0].ld_coeffs = [2.2,0,1.0,0];
q.bodies[0].ld_powers = [0.0, 0.5, 1.0, 1.5];

q.bodies[1] = OrbitingSystem();
q.bodies[1].mass = 1.;
q.bodies[1].semimajor = 1.;
q.bodies[1].temperature = 0.015;
q.bodies[1].inclination = np.pi * 90./180.
q.bodies[1].radius = .22623 * 0.00464913034;
q.bodies[1].phase = -np.pi/8 + np.pi/30.;
q.bodies[1].eccentricity = 0.15944;
q.bodies[1].arg_periastron = -np.pi/30;
q.bodies[1].radius = 0.000477894503 * 0.0892141778;

q.setTotalMass();
"""
#Masses of planetary system
r = OrbitingSystem();
r.semimajor = 0.7048;
r.eccentricity = 0.0069;
r.inclination = np.pi * 90./180.;
r.phase = np.pi * -20./180.

r.bodies = [0, 0];
r.bodies[0] = OrbitingSystem();
r.bodies[0].mass = 317.83 * 0.333;
r.bodies[0].radius = 0.000477894503 * 0.7538;

#Nominal moon
r.bodies[1] = OrbitingSystem();
r.bodies[1].mass = 0.00001;
r.bodies[1].semimajor = 1.;
r.bodies[1].radius = 0.0000000001;
#r.bodies[1].inclination = np.pi * 0.5;

r.setTotalMass();
"""
s = q;
s.setTotalMass();
################################################
########### END PLANETARY DATA HERE  ###########
################################################
