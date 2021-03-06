if __package__ is None:
	import sys;
	from os import path;
	sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) );
	print sys.path[-1];
	from gravsims.basicKepler import *;
	#from processing_pipeline.pipeline import import_idl;	
else:
	from ..gravsims.basicKepler import *;
	#from ..processing_pipeline.pipeline import import_idl;

import numpy as np;
from itertools import chain;
import matplotlib.pyplot as plt;
import general_transits as gt;
import time;
import copy;

def getBodies(datafiles = ["K16"]):
	'''This function is used to divide the ability to import existing datafiles from the capacity
	to create new ones. In the future, there will simply be a "base tree" which is imported and 
	other bodies to be apended to the various .bodies attributes of the tree.'''

	finalBodies = [];

	for u in range(0, len(datafiles)):
		#Ensuring that import is always possible in any case
		if __package__ is None:
			datafile = __import__("transit_models.systems." + datafiles[u], fromlist=['']);
		else:
			datafile = __import__("systems." + datafiles[u], fromlist = ['']);
		
		system = datafile.s;
		finalBodies.append(system);

	return finalBodies;

def genTransits(data, steps = 7501., revolutions = 1., n = 81):
	'''This function can be used to loop through any number of datafiles containing an OrbitingSystem
	class instance with all orbital parameters. Although it is not implemented here, the same frame- 
	work could be used theoretically to tweak orbital parameters in order to perform MCMC. Unlike pre-
	vious versions of the code, this means that orbital parameters need no longer be static.

	By default, the program uses the first system's orbital period as the "base time".'''	

	lightcurves = [];

	revmax = np.max(np.array(revolutions));
	print revmax;

	year = 1.;
	timecoords = np.zeros(steps) + .0;

	for u in range(0, len(data)):
	#	l = time.time();
		system = data[u];
		if u == 0:
			#always make sure there is a consistent time basis to compare transits
			system.setBaseTime();
			year = system.bt;
		else:
			system.bt = year;

		#Generating coordinates relative to COM of each system
		system.setSystemOrbits(s = steps, r = revmax);
		times = system.times;

		for m in system.bodies:
			m.bt = year;
			#if len(m.bodies) > 1:
			m.setSystemOrbits(s = steps, r = revmax, verbose = True);

		#Getting final absolute X, Y, and Z coordinates in order to generate lightcurve
		l = time.time();
		final_x, final_y, final_z, transit_array = gt.traverse_tree(system, times);	
		fx, fy, fz, zpos = gt.sort_keys(final_x, final_y, final_z);
		lb, intersects = gt.arrange_combinations(fx, fy, transit_array);

		if type(revolutions) is not (float or int):
			lccr = np.sum([((times/max(times) * revmax >= a[0]) & (times/max(times) * revmax <= a[1])) for a in revolutions], axis=0).astype('bool');
		else:
			lccr = np.ones(len(times)).astype('bool');

		#print lccr;

		#print len(fx), len(fy), len(lb), len(intersects), len([n]), len(transit_array), len(times), len(zpos);

		cadence, light_blocked = gt.generate_lightcurve(fx[lccr], fy[lccr], lb[lccr], intersects[lccr], n, transit_array, times[lccr], zpos[lccr]);
		lightcurves.append(light_blocked);
		if (u == 0):
			timecoords = cadence
		print time.time() - l;
		print "";
	return timecoords, lightcurves;

def kepler16AlterPlanet():
	s = getBodies()[0];

	#Actual with planet only
	s.setTotalMass();
	bodiesList = [];
	bodiesList.append(s);

	#Actual with moon
	m = copy.deepcopy(s);
	m.modifyBody(position=[1,1], mass=1., semimajor = 0.005, radius = 0.000477894503 * 0.089214177);
	m.setTotalMass();
	bodiesList.append(m);

	#Predicted with moon
	for u in np.linspace(-50,50,101)/2.5:
		t = copy.deepcopy(s);
		t.modifyBody(position=[1], semimajor = 0.7048 + 0.001 * u);
		t.setTotalMass();
		bodiesList.append(t);

	finals = genTransits(bodiesList, revolutions = 11./220., steps = 5001, n = 51);

	return finals;

def kepler16AlterMoon():
	s = getBodies()[0];

	bodiesList = [];
	s.setTotalMass();
	bodiesList.append(s);

	#Predicted with moon
	for u in np.linspace(0.00000000001, 1, 8):
		t = copy.deepcopy(s);
		t.modifyBody(position=[1, 1], mass = 1. * u, semimajor = 0.005, radius = 0.000477894503 * 0.0892141778 * u**(1./3.));
		t.setTotalMass();
		bodiesList.append(t);

	finals = genTransits(bodiesList, revolutions = 11./220., steps = 15001, n = 81);

	return finals;

