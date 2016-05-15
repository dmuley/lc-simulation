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

def genTransits(datafiles=["K16"], steps = 7501., revolutions = 1., n = 31):
	'''This function can be used to loop through any number of datafiles containing an OrbitingSystem
	class instance with all orbital parameters. Although it is not implemented here, the same frame- 
	work could be used theoretically to tweak orbital parameters in order to perform MCMC. Unlike pre-
	vious versions of the code, this means that orbital parameters need no longer be static.

	By default, the program uses the first system's orbital period as the "base time".'''	

	year = 1.;
	timecoords = np.zeros(steps) + .0;

	lightcurves = [];

	for u in range(0, len(datafiles)):
		l = time.time()
		#Ensuring that import is always possible in any case
		if __package__ is None:
			datafile = __import__("transit_models.systems." + datafiles[u], fromlist=['']);
		else:
			datafile = __import__("systems." + datafiles[u], fromlist = ['']);
		print datafile;
		system = datafile.s;
		if u == 0:
			#always make sure there is a consistent time basis to compare transits
			system.setBaseTime();
			year = system.bt;
		else:
			system.bt = year;

		#Generating coordinates relative to COM of each system
		system.setSystemOrbits(s = steps, r = revolutions);
		times = system.times;

		for m in system.bodies:
			m.bt = year;
			#if len(m.bodies) > 1:
			m.setSystemOrbits(s = steps, r = revolutions);

		#Getting final absolute X, Y, and Z coordinates in order to generate lightcurve
		final_x, final_y, final_z, transit_array = gt.traverse_tree(system, times);	
		fx, fy, fz, zpos = gt.sort_keys(final_x, final_y, final_z);
		lb, intersects = gt.arrange_combinations(fx, fy, transit_array);

		cadence, light_blocked = gt.generate_lightcurve(fx, fy, lb, intersects, n, transit_array, times, zpos);
		lightcurves.append(light_blocked);
		if (u == 0):
			timecoords = cadence
		print time.time() - l;
	return timecoords, lightcurves;